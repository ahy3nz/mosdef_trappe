
import unyt as u
import mbuild as mb
import foyer

from mosdef_trappe.driver import build
from mosdef_trappe.molecules.ethane.ethane import Ethane
cmpd = Ethane()
structure = build(cmpd, n_compounds=5000, density=19.999999999999996*u.kg/u.m**3)
second = build(cmpd, n_compounds=5000, density=468.99999999999994*u.kg/u.m**3)

import mosdef_trappe.gomc_util.gomc_functions as gomc_funcs
p = gomc_funcs.simulate(structure, second, temperature=236.0*u.Kelvin,
         n_steps=5000000)
while p.returncode != 0:
    p = gomc_func.simulate(structure, second, temperature=236.0*u.Kelvin,
        n_steps=5000000)
             