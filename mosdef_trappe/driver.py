import unyt as u
import mbuild as mb
import numpy as np
import foyer

def main():
    from mosdef_trappe.molecules.propane.propane import Propane
    cmpd = Propane()
    structure = build(cmpd, n_compounds=500, density=615.5*u.kg/u.m**3)
    second = build(cmpd, n_compounds=500, density=115.5*u.kg/u.m**3)

    #import mosdef_trappe.gomc_util.gomc_functions as gomc_funcs
    #gomc_funcs.simulate(structure, second, temperature=300*u.Kelvin,
            #pressure=5*u.atm, n_steps=5000000)
    simulate(structure, 'gromacs', 
            temperature=300*u.Kelvin, pressure=5*u.atm,
            n_steps=50000000)

    #import panedr
    #df = panedr.edr_to_df('md.edr')
    #print(df['Pressure'])


def build(cmpd, density=0.5*u.gram/(u.cm**3), n_compounds=1000, 
        ff='ff/TraPPE_UA_3_fully_flexible_propane.xml'):
    """ Build and parametrize a TraPPE system at a given state

    Parameters
    ----------
    cmpd : mb.Compound
        This compound should also have the path to its XML as 
        `cmpd.xml`
    density : unyt.Quantity
        Desired density for packing and siulation
    n_compounds : int


    """
    density.convert_to_units(u.kilogram/u.m**3)

    # Pack a box
    box = mb.fill_box(cmpd, n_compounds=n_compounds, density=density.value)

    # Wrap coordinates
    new_xyz = box.xyz - 1 * np.floor_divide(box.xyz, box.periodicity) * box.periodicity
    box.xyz = new_xyz

    # Apply non-atomistic, custom element naming convention
    for part in box.particles():
        part.name = "_" + part.name

    # Utilize foyer to parametrize our box
    ff = foyer.Forcefield(forcefield_files=ff)
    box = box.to_parmed(infer_residues=True)
    parametrized_structure = ff.apply(box, combining_rule='lorentz')

    # Dump initial coordinates
    parametrized_structure.save('compound.pdb', overwrite=True)
    parametrized_structure.save('compound.mol2', overwrite=True)
    parametrized_structure.save('compound.gro', overwrite=True)

    return parametrized_structure

def simulate(structure, engine, **kwargs):
    """ Simulate a structure in a particular engine

    Parameters
    ---------
    structure : parmed.Structure
    engine : str
        'gromacs', 'hoomd', 'openmm', 'gomc'
    **kwargs : run-control kwargs for particular engine

    Notes
    -----
    Some run-control parameters will not be found inside this routine.
    Look inside the `X_utils` to further inspect the run-control
    parameters for engine X.
    """
    if engine.lower() == 'gromacs':
        import mosdef_trappe.gromacs_util.gromacs_functions as gmx_funcs
        gmx_funcs.simulate(structure, **kwargs) 
    elif engine.lower() == 'openmm':
        import mosdef_trappe.openmm_util.openmm_functions as omm_funcs
        omm_funcs.simulate(structure, **kwargs) 
    elif engine.lower() == 'hoomd':
        import mosdef_trappe.hoomd_util.hoomd_functions as hoomd_funcs
        hoomd_funcs.simulate(structure, **kwargs) 
    elif engine.lower() == 'gomc':
        import mosdef_trappe.gomc_util.gomc_functions as gomc_funcs
        gomc_funcs.simulate(structure, **kwargs)


if __name__ == "__main__":
    main()


