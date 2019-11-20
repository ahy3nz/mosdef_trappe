import unyt as u
import mbuild as mb
import foyer
import pdb

import molecules.propane.propane as propane


# This is the compound we want to study
cmpd = propane.Propane()

# These are some basic physical properties for packing purposes
density = 0.5 * u.gram/(u.cm**3)
density.convert_to_units(u.kilogram/u.m**3)

# Pack a box
box = mb.fill_box(cmpd, n_compounds=1000, density=density.value)
from mbuild.utils.geometry import wrap_coords
new_xyz = wrap_coords(box.xyz, box.periodicity)
box.xyz = new_xyz

# Apply non-atomistic, custom element naming convention
for part in box.particles():
    part.name = "_" + part.name

# Utilize foyer to parametrize our box
ff = foyer.Forcefield(forcefield_files=cmpd.xml)
parametrized_structure = ff.apply(box, combining_rule='lorentz')


# Implementing model for gromacs
#import gromacs_util.gromacs_functions as gmx_functions
#parametrized_structure.save('out.gro', overwrite=True)
#parametrized_structure.save('out.top', overwrite=True)
#gmx_functions.write_npt_mdp('npt.mdp', temperature=300*u.Kelvin, pressure=1*u.atm,
#                        n_steps=500, random_seed=42)
#
## Implementing model for openmm
#import openmm_util.openmm_functions as omm_functions
#sim = omm_functions.build_openmm_simulation(parametrized_structure, 
#        temperature=300*u.Kelvin,
#        pressure=1*u.atm, random_seed=42)
#omm_functions.run_openmm_simulation(n_steps=500)

# Implementing model for hoomd
from mbuild.formats.hoomd_simulation import create_hoomd_simulation
import hoomd_util.hoomd_functions as hoomd_functions
create_hoomd_simulation(parametrized_structure,
        ref_distance=10, ref_energy=1/4.184, r_cut=1.4)
hoomd_functions.run_hoomd_simulation(n_steps=500)
