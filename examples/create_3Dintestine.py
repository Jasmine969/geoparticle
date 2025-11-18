import geoparticle as gp
from lammps import lammps
import numpy as np

rho_fluid = 993
rho_wall = 1040
l_pipe_hrz = 0.002
l_pipe_vert = 0.004
r_si = 0.002
r_torus = 2 * r_si
dl = 2e-4
l_si = 0.02 - dl

# geometry creation ======================
# create solids
pipe_in_vert = gp.CylinderSide(
    r=r_si,  # radius of the cylinder
    l_axis=l_pipe_vert,  # length along the cylinder axis
    dl=dl,  # particle spacing
    axis='z'  # cylinder axis direction
).shift(x=-r_torus - l_pipe_hrz, z=r_torus)
l_pipe_vert = pipe_in_vert.l_axis
n_ring = pipe_in_vert.n_ring

torus_in = gp.TorusSurface(
    r_minor=r_si,  # minor radius
    r_major=r_torus,  # major radius
    dl=dl,  # particle spacing
    n_ring=n_ring,  # number of particle rings along minor circle
    phi_range='(180,270)',  # angle range along major circle
    plane='XOZ'  # plane where the major circle lies
).shift(z=r_torus, x=-l_pipe_hrz)
pipe_in_hrz = gp.CylinderSide(
    r=r_si, l_axis=l_pipe_hrz - dl, dl=dl, axis='x', name='pipe_in_hrz'
).shift(x=-l_pipe_hrz)
inlet = gp.Union((pipe_in_vert, torus_in, pipe_in_hrz), name='inlet')
# inlet = pipe_in_vert + torus_in + pipe_in_hrz

si = gp.CylinderSide(r=r_si, l_axis=l_si, dl=dl, axis='x')
l_si = si.l_axis

outlet = inlet.mirror(plane_name='YOZ', plane_pos=l_si / 2)
# create fluid
fluid = gp.FilledCylinder(r_si - dl, l_si, dl, axis='x', name='fluid')
# create atoms in lammps ===========================
lmp = lammps(cmdargs=['-screen', 'none', '-log', 'none'])
# This buffer is between the geometry bound and the bound of the simulation box
buffer_region = 0.1 * l_si
xlo = -(l_pipe_hrz + r_torus + r_si) - buffer_region
xhi = l_si + l_pipe_hrz + r_torus + r_si + buffer_region
ylo = -r_si - buffer_region
yhi = r_si + buffer_region
zlo = -r_si - buffer_region
zhi = l_pipe_vert + r_torus + buffer_region

lmp.commands_string(f"""
dimension	3
atom_style  rheo
units		si
newton	 	on
boundary	f f f
comm_modify vel yes
region      simulation_box block {xlo} {xhi} {ylo} {yhi} {zlo} {zhi} 
create_box  3 simulation_box
""")
n_atoms_inlet = inlet.size
lmp.create_atoms(
    n_atoms_inlet, np.arange(n_atoms_inlet) + 1 + lmp.get_natoms(),
    np.ones(n_atoms_inlet, dtype=int), inlet.flatten_coords)
n_atoms_si = si.size
lmp.create_atoms(
    n_atoms_si, np.arange(n_atoms_si) + 1 + lmp.get_natoms(),
    np.full(n_atoms_si, 2, dtype=int), si.flatten_coords)
n_atoms_outlet = outlet.size
lmp.create_atoms(
    n_atoms_outlet, np.arange(n_atoms_outlet) + 1 + lmp.get_natoms(),
    np.ones(n_atoms_outlet, dtype=int), outlet.flatten_coords)
n_atoms_fluid = fluid.size
lmp.create_atoms(
    n_atoms_fluid, np.arange(n_atoms_fluid) + 1 + lmp.get_natoms(),
    np.full(n_atoms_fluid, 3, dtype=int), fluid.flatten_coords)
n_atoms_all = lmp.get_natoms()
log_n_atom = f'n_atoms_inlet: {n_atoms_inlet},' \
             f' n_atoms_si: {n_atoms_si},' \
             f' n_atoms_outlet: {n_atoms_outlet},' \
             f' n_atoms_fluid: {n_atoms_fluid},' \
             f' n_atoms_all: {n_atoms_all}.'
print(log_n_atom)

lmp.commands_string(f"""
mass            * 1
set             group all rheo/rho {rho_fluid}
""")

filename = 'intestine3D'
lmp.commands_string(f"""
pair_style      zero {dl * 2}
pair_coeff      * *
neighbor        {dl * 0.1} bin
run             0
write_data      {filename}.data
""")

# check atom overlapping
lmp.commands_string(f"""
delete_atoms    overlap {dl * 0.8} all all
""")
n_atoms_all_now = lmp.get_natoms()
if n_atoms_all == n_atoms_all_now:
    print('Congrats, no atoms are overlapped!')
else:
    raise RuntimeError(f'{n_atoms_all - n_atoms_all_now} atoms are overlapped!')
