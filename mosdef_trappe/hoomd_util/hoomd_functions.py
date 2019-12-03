import unyt

import hoomd
import hoomd.md
kb = unyt.kb.in_units(unyt.amu*unyt.nm**2/(unyt.Kelvin*unyt.ps**2))

def simulate(parametrized_structure, **kwargs):
    """ Simulate in Hoomd 

    Parameters
    ----------
    parametrized_structure : parmed.Structure
    **kwargs : kwargs for run_hoomd_simulation
    
    """
    from mbuild.formats.hoomd_simulation import create_hoomd_simulation
    create_hoomd_simulation(parametrized_structure,
            ref_distance=10, ref_energy=1/4.184, r_cut=1.4)
    run_hoomd_simulation(**kwargs)


def run_hoomd_simulation(temperature=300*unyt.Kelvin,
        pressure=1*unyt.atm, n_steps=500000, random_seed=42):
    all_group = hoomd.group.all()
    hoomd.context.current.neighbor_lists[0].reset_exclusions(
            ['1-2', '1-3', '1-4'])
    fire = hoomd.md.integrate.mode_minimize_fire(0.002)
    nve_integrator = hoomd.md.integrate.nve(group=all_group)
    while not (fire.has_converged()):
        hoomd.run(1000)
    nve_integrator.disable()

    hoomd.md.integrate.mode_standard(dt=0.002)

    temperature.convert_to_units(unyt.Kelvin)
    if pressure is not None:
        pressure.convert_to_units(unyt.amu/(unyt.nm*unyt.ps**2))
        integrator = hoomd.md.integrate.npt(all_group, 
                kT=kb.value * temperature.value, tau=0.4, 
                P=pressure.value, tauP=2.0, couple='xyz')
    else:
        integrator = hoomd.md.integrate.nvt(all_group,
                kT=kb.value * temperature.value, tau=0.4)

    hoomd.analyze.log('thermo.log', 
            ['potential_energy', 'kinetic_energy', 'temperature', 'pressure'],
            5000, header_prefix='PE, KE, T, P', overwrite=True)
    hoomd.dump.dcd('md.dcd', 5000, all_group, overwrite=True)

    integrator.randomize_velocities(random_seed)
    hoomd.run(n_steps)

