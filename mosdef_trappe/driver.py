import unyt as u

import build
import molecules.propane.propane as propane


cmpd = propane.Propane()
temperature = 300*u.Kelvin
pressure = None
n_compounds = 1000
build.build_simulate(cmpd, temperature=temperature, pressure=pressure,
        n_compounds=n_compounds, engine='gromacs')
