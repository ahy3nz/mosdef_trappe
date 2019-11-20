import unyt as u
import mbuild as mb
import foyer

def main():
    from mosdef_trappe.molecules.propane.propane import Propane
    import panedr
    cmpd = Propane()
    build_simulate(cmpd, temperature=200*u.Kelvin, pressure=None,
            density=615.5*u.kg/u.m**3, n_compounds=500, n_steps=10000)
    df = panedr.edr_to_df('md.edr')
    print(df['Pressure'])

def build_simulate(cmpd, temperature=300*u.Kelvin, pressure=None,
        density=0.5*u.gram/(u.cm**3), n_compounds=1000, 
        random_seed=42, engine='gromacs', n_steps=500000):
    """ Build and simulate a TraPPE compound at a given state

    This is the primary work function in this repo to 
    simulate and sample from a particular thermodynamic state.

    Parameters
    ----------
    cmpd : mb.Compound
        This compound should also have the path to its XML as 
        `cmpd.xml`
    temperature : unyt.Quantity
    pressure : unyt.Quantity
        If None, will run NVT
    density : unyt.Quantity
        Desired density for packing and siulation
    n_compounds : int
    random_seed : int
    engine : str
        'gromacs', 'hoomd', 'openmm'
    n_steps : int

    Notes
    -----
    Some run-control parameters will not be found inside this routine.
    Look inside the `X_utils` to further inspect the run-control
    parameters for engine X.
    While running this singular function once will work well for 
    demonstrations, see `script_util` for some modules that
    will manage and perform simulations of many state points
    """
    density.convert_to_units(u.kilogram/u.m**3)
    temperature.convert_to_units(u.Kelvin)

    # Pack a box
    box = mb.fill_box(cmpd, n_compounds=n_compounds, density=density.value)
    from mbuild.utils.geometry import wrap_coords
    new_xyz = wrap_coords(box.xyz, box.periodicity)
    box.xyz = new_xyz

    # Apply non-atomistic, custom element naming convention
    for part in box.particles():
        part.name = "_" + part.name

    # Utilize foyer to parametrize our box
    ff = foyer.Forcefield(forcefield_files=cmpd.xml)
    box = box.to_parmed(infer_residues=True)
    parametrized_structure = ff.apply(box, combining_rule='lorentz')
    parametrized_structure.save('compound.mol2', overwrite=True)

    # Implementing model for gromacs
    if engine.lower() == 'gromacs':
        import mosdef_trappe.gromacs_util.gromacs_functions as gmx_functions
        parametrized_structure.save('compound.gro', overwrite=True)
        parametrized_structure.save('compound.top', overwrite=True)

        gmx_functions.write_em_mdp('em.mdp')
        grompp = gmx_functions.run_grompp(mdp='em.mdp',
                gro='compound.gro', top='compound.top', out='em')
        if grompp.returncode == 0:
            mdrun = gmx_functions.run_md('em')

        gmx_functions.write_mdp('md.mdp', 
                temperature=temperature, 
                pressure=pressure, 
                n_steps=n_steps, random_seed=42)
        grompp = gmx_functions.run_grompp(mdp='md.mdp',
                gro='em.gro', top='compound.top', out='md')
        if grompp.returncode == 0:
            mdrun = gmx_functions.run_md('md')
    
    # Implementing model for openmm
    if engine.lower() == 'openmm':
        import mosdef_trappe.openmm_util.openmm_functions as omm_functions
        sim = omm_functions.build_openmm_simulation(parametrized_structure, 
                temperature=temperature,
                pressure=pressure, random_seed=random_seed)
        omm_functions.run_openmm_simulation(n_steps=n_steps)

    # Implementing model for hoomd
    if engine.lower() == 'hoomd':
        from mbuild.formats.hoomd_simulation import create_hoomd_simulation
        import mosdef_trappe.hoomd_util.hoomd_functions as hoomd_functions
        create_hoomd_simulation(parametrized_structure,
                ref_distance=10, ref_energy=1/4.184, r_cut=1.4)
        hoomd_functions.run_hoomd_simulation(temperature=temperature,
                pressure=pressure,
                random_seed=random_seed,
                n_steps=n_steps)

if __name__ == "__main__":
    main()

