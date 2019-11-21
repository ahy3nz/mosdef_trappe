import simtk.unit as unit
import simtk.openmm as openmm
import simtk.openmm.app as app
import unyt

def simulate(parametrized_structure, n_steps=500000, **kwargs):
    """ Umbrella routine for running openmm simulation 
    
    Parameters
    ----------
    **kwargs : kwargs for build_openmm_simulation """
    sim = build_openmm_simulation(parametrized_structure, **kwargs)

    run_openmm_simulation(sim, n_steps=n_steps)


def build_openmm_simulation(structure, temperature=300*unyt.Kelvin, 
        pressure=1*unyt.atm, random_seed=42, **kwargs):
    """ Build OpenMM simulation from a parmed.Structure 
    
    
    Notes
    -----
    OpenMM does not compute a virial, which prevents us from
    computing and reporting the pressure of a system. 
    However, the montecarlobarostat does allow for a robust method to sample
    under various pressures.
    """

    # First convert unyt units into something consistent for OpenMM
    temperature.convert_to_units(unyt.Kelvin)
    if pressure is not None:
        pressure.convert_to_units(unyt.bar)

    # hardcoded timestep - probably do not want to expose
    timestep = 2.0 * unyt.fs
    timestep.convert_to_units(unyt.ps)

    integrator = openmm.LangevinIntegrator(float(temperature.value), 
                                    float((1.0/unyt.ps).value), 
                                    float(timestep.value))
    integrator.setRandomNumberSeed(random_seed)


    system = structure.createSystem(nonbondedMethod=app.PME,
                                constraints=app.AllBonds,
                                nonbondedCutoff=14.0*unit.angstroms)
    if pressure is not None:
        barostat = openmm.MonteCarloBarostat(float(pressure.value), 
                float(temperature.value), 25)
        system.addForce(barostat)

    sim = app.Simulation(structure.topology, system, integrator)
    sim.context.setPositions(structure.positions)

    sim.reporters.append(app.StateDataReporter(open('thermo.log','w'), 5000, 
                                                step=True, time=True,
                                                potentialEnergy=True,
                                                temperature=True,
                                                volume=True, speed=True))
    sim.reporters.append(app.DCDReporter('trajectory.dcd', 5000))
    sim.reporters.append(app.CheckpointReporter('trajectory.chk', 5000))

    return sim


def run_openmm_simulation(sim, n_steps=500000):
    sim.minimizeEnergy(maxIterations=50000)
    sim.step(n_steps)

    pdbreporter = app.PDBReporter('final_frame.pdb', 1)
    pdbreporter.report(sim, sim.context.getState(-1))
