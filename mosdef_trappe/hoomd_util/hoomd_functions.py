import unyt

import hoomd
import hoomd.md
kb = unyt.kb.in_units(unyt.amu*unyt.nm**2/(unyt.Kelvin*unyt.ps**2))

def run_hoomd_simulation(temperature=300*unyt.Kelvin,
        pressure=1*unyt.atm, n_steps=500000, random_seed=42):
    all_group = hoomd.group.all()
    hoomd.context.current.neighbor_lists[0].reset_exclusions(
            ['1-2', '1-3', '1-4'])
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
    integrator.randomize_velocities(random_seed)

    hoomd.run(n_steps)
