import numpy as np
import pandas as pd

import subprocess
import unyt as u

import parmed


def simulate(*args, **kwargs):
    """ Simulate using GOMC
    
    Parameters
    ----------
    *args : Should be 1 or 2 parmed.Structure
    **kwargs: kwargs for writing gomc input
    """
    if len(args) == 1:
        simulate_single(*args, **kwargs)
    elif len(args) ==2:
        simulate_double(*args, **kwargs)

def simulate_single(parametrized_structure, **kwargs):
    """ Simulate a single box in GOMC

    Parameters
    ---------
    parametrized_structure : parmed.Structure
    **kwargs : kwargs for writing gomc input
    """
    # Save PDB and PSF files
    parametrized_structure.save('coords.pdb', overwrite=True, use_hetatoms=False)
    parametrized_structure.save('structure.psf', overwrite=True)

    # Convert atomtypes to uppercase
    for atom in parametrized_structure.atoms:
        atom.type = atom.type.upper()
        atom.atom_type.name = atom.atom_type.name.upper()

    # Prepare parmed Structure by converting to CharmmParameterSet
    paramset = parmed.charmm.CharmmParameterSet.from_structure(
            parametrized_structure)
    paramset.write(par='parameters.par')

    # Post-process PAR file
    modify_par_file('parameters.par')


    # Write GOMC run file, which requires box information and 
    # knowledge of other input files
    write_gomc_single_input('in.conf', parametrized_structure, 
            coords='coords.pdb', structure='structure.psf', 
            parameters='parameters.par', **kwargs)

    run_single(**kwargs)

def simulate_double(parametrized_structure, other_structure, **kwargs):
    """ Simulate two boxes in GOMC

    Parameters
    ---------
    parametrized_structure : parmed.Structure
    other_structure : parmed.Structure
    **kwargs : kwargs for writing gomc input
    """
    # Save PDB and PSF files
    parametrized_structure.save('coords0.pdb', overwrite=True, use_hetatoms=False)
    parametrized_structure.save('structure0.psf', overwrite=True)
    other_structure.save('coords1.pdb', overwrite=True, use_hetatoms=False)
    other_structure.save('structure1.psf', overwrite=True)


    # Convert atomtypes to uppercase
    for atom in parametrized_structure.atoms:
        atom.type = atom.type.upper()
        atom.atom_type.name = atom.atom_type.name.upper()
    for atom in other_structure.atoms:
        atom.type = atom.type.upper()
        atom.atom_type.name = atom.atom_type.name.upper()


    # Prepare parmed Structure by converting to CharmmParameterSet
    paramset = parmed.charmm.CharmmParameterSet.from_structure(
            parametrized_structure)
    paramset.write(par='parameters.par')

    # Post-process PAR file
    modify_par_file('parameters.par')


    # Write GOMC run file, which requires box information and 
    # knowledge of other input files
    write_gomc_double_input('in.conf', parametrized_structure, other_structure,
            coords0='coords0.pdb', structure0='structure0.psf', 
            coords1='coords1.pdb', structure1='structure1.psf', 
            parameters='parameters.par', **kwargs)

    run_double(**kwargs)

def run_single(**kwargs):
    if kwargs.get('pressure', None):
        gomc_bin = 'GOMC_CPU_NPT'
    else:
        gomc_bin = 'GOMC_CPU_NVT'
    p = subprocess.Popen('{} in.conf'.format(gomc_bin), shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True)
    out, err = p.communicate()
    with open("gomc.out", 'w') as f:
        f.write(out)
    with open("gomc.err", 'w') as f:
        f.write(err)

def run_double(**kwargs):
    gomc_bin = 'GOMC_CPU_GEMC'
    p = subprocess.Popen('{} in.conf'.format(gomc_bin), shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True)
    out, err = p.communicate()
    with open("gomc.out", 'w') as f:
        f.write(out)
    with open("gomc.err", 'w') as f:
        f.write(err)
        
def dat_to_df(filename):
    """ Convert GOMC DAT file to pandas DataFrame
    
    Notes
    -----
    TOT_DENS [kg/m**3] converted to [g/cm**3]
    TOT_EN [K] converted to [kJ/mol]
    """
    data = np.loadtxt(filename)
    if len(data) == 0:
        return None
    columns = open(filename, 'r').readlines()[0]
    df = pd.DataFrame(data, columns=columns.split())
    
    mass = 30.07 * u.g
    df['TOT_DENS'] = df['TOT_DENS']*(1*u.kg/u.m**3).in_units(u.g/u.cm**3).value
    df['TOT_EN'] = df['TOT_EN']*(1*u.Kelvin*u.boltzmann_constant*6.022e23).in_units(u.kilojoule).value
    return df

def modify_par_file(par_file):
    """ GOMC parameter files do not use the atoms or impropers directive"""
    parlines = open(par_file, 'r').readlines()
    modified_parlines = open(par_file, 'r').readlines()
    found_atom_section = False
    to_delete = []
    for i, line in enumerate(parlines):
        if "ATOMS" in line:
            found_atom_section = True
            to_delete.append(i)
        elif found_atom_section:
            if "MASS" in line:
                to_delete.append(i)
            else:
                found_atom_section = False

        if 'IMPROPERS' in line:
            modified_parlines[i]='IMPROPER\n'
    to_delete.reverse()
    for delete_index in to_delete:
        modified_parlines.pop(delete_index)
    with open(par_file, 'w') as f:
        for line in modified_parlines:
            f.write(line)


def write_gomc_double_input(filename, param_structure0, param_structure1, 
                            temperature=273*u.Kelvin, pressure=None,
                            parameters='parameters.par',
                            coords0='coords0.pdb', structure0='structure0.psf', 
                            coords1='coords1.pdb', structure1='structure1.psf',
                            output='output', n_steps=500000):
    """ Write GOMC input file for a GEMC simulation 
    
    Some simulation input parameters have been hard-coded in accordance with the 
    TraPPE specification. 
    MC moves, CBMC parameters, and output control have been adapted from
    https://github.com/GOMC-WSU/GOMC_Examples/tree/master/NVT_GEMC/pure_fluid/octane_T_360_00_K
    Gibbs Ensemble Monte Carlo involves simultaneous simulations of 
    two boxes (generally corresponding to two different phases)
    
    Parameters
    ----------
    filename : str
    param_structure0: parmed.Structure
        Should be parametrized
    param_structure1: parmed.Structure
        Should be parametrized
    temperature : unyt.Quantity
        Will get converted to K
    pressure : unyt.Quantityt, default None
        Will get converted to bar if not None
    parameters : str
        CHARMM/NAMD-style PAR file for force field information
    coords0 : str
        PDB file for param_structure0 coordinates and atom names
    structure0: str
        PSF file for param_structure0 atomtypes, bonding information, and topology
    coords1 : str
        PDB file for param_structure0 coordinates and atom names
    structure1: str
        PSF file for param_structure0 atomtypes, 
        bonding information, and topology    
    output : str
        Prefix for output files from simulation
    n_steps : int
        
    Notes
    -----
    Box dimensions are assumed to be orthorhombic
    GOMC input files utilize Angstroms, which are also the units of the parmed Structure box
    TraPPE uses 14 Angstrom cutoffs, but here we use 10 Angstrom because our simulations are small
    """
    if pressure is None:
        gemc_type = "GEMC   NVT"
    else:
        gemc_type = "GEMC   NPT\nPressure {}\n".format(pressure.in_units(
            u.bar).value)

    with open(filename, 'w') as f:
        f.write(""" 
###########################################################################
#  ========-------------------- INPUT --------------------------===========
############################################################################

#########################
# enable, step
#########################
Restart     false


####################################
# kind (RESTART, RANDOM, INTSEED)
####################################
PRNG        RANDOM

####################################
# FORCE FIELD
####################################
ParaTypeCHARMM   true
ParaTypeEXOTIC   false
Parameters      {parameters}

####################################
# INPUT PDB FILES
####################################
Coordinates 0    {coords0}
Coordinates 1    {coords1}

####################################
# INPUT PSF FILES
####################################
Structure 0      {structure0}
Structure 1      {structure1}



############################################################################
#  =======--------------------- SYSTEM --------------------------===========
############################################################################

##################################
# GEMC TYPE (DEFULT IS NVT_GEMC)
##################################
{gemc_type}

#############################
# SIMULATION CONDITION
#############################
Temperature     {temperature}
Potential       VDW
LRC     true
Rcut        14
Exclude     1-4

#############################
# ELECTROSTATIC
#############################
ElectroStatic   false
Ewald           false

################################
# PRESSURE FREQ
################################
PressureCalc  true  1000

################################
# STEPS
################################
RunSteps       {n_steps}
EqSteps        500000
AdjSteps       1000

################################
# MOVE FREQUENCY
################################
DisFreq               0.49
RotFreq           0.10
VolFreq           0.01
SwapFreq          0.20
RegrowthFreq          0.10
CrankShaftFreq        0.10

################################
# BOX DIMENSION #, X, Y, Z
################################
CellBasisVector1 0  {box0_x}  0.00    0.00
CellBasisVector2 0  0.00    {box0_y}  0.00
CellBasisVector3 0  0.00    0.00    {box0_z}

CellBasisVector1 1  {box1_x}  0.00    0.00
CellBasisVector2 1  0.00    {box1_y}  0.00
CellBasisVector3 1  0.00    0.00    {box1_z}


##############################
# CBMC TRIALS
##############################
CBMC_First   10
CBMC_Nth     8
CBMC_Ang     100
CBMC_Dih     30

####################################
#          Mol. Name     Chem. Pot.
####################################
############################################################################
#  =======-------------------- OUTPUT --------------------------===========
############################################################################

##########################
# statistics filename add
##########################
OutputName  {output}

#####################################
# enable, frequency
#####################################
CoordinatesFreq    true   1000000
RestartFreq        true   1000000
ConsoleFreq        true   100000
BlockAverageFreq   true   100000
HistogramFreq      false  100000


################################
# OutHistSettings
################################


##################################
# enable: blk avg., fluct.
##################################
OutEnergy         true    true
OutPressure       true    true
OutMolNum         true    true
OutDensity        true    true
""".format(temperature=temperature.in_units(u.Kelvin).value,
            parameters=parameters, gemc_type=gemc_type,
           structure0=structure0, coords0=coords0, 
           structure1=structure1, coords1=coords1,
           box0_x=param_structure0.box[0], box0_y=param_structure0.box[1], 
           box0_z=param_structure0.box[2],
           box1_x=param_structure1.box[0], box1_y=param_structure1.box[1], 
           box1_z=param_structure1.box[2],
           output=output, n_steps=n_steps))

def write_gomc_single_input(filename, parm_structure, coords='coords.pdb',
                     structure='structure.psf', parameters='parameters.par',
                    temperature=305*u.Kelvin, pressure=None,
                    output='output', n_steps=500000):
    """ Write GOMC input file for a single-box simulation 
    
    Some simulation input parameters have been hard-coded in accordance with the 
    TraPPE specification. 
    MC moves, CBMC parameters, and output control have been adapted from
    https://github.com/GOMC-WSU/GOMC_Examples/tree/master/NVT_GEMC/pure_fluid/octane_T_360_00_K
    
    Parameters
    ----------
    filename : str
    parm_structure: parmed.Structure
        Should be parametrized
    temperature : unyt.Quantity
        Will get converted to K
    pressure : unyt.Quantityt, default None
        Will get converted to bar if not None
    parameters : str
        CHARMM/NAMD-style PAR file for force field information
    coords : str
        PDB file for param_structure coordinates and atom names
    structure: str
        PSF file for param_structure atomtypes, bonding information, and topology
    output : str
        Prefix for output files from simulation
        
    Notes
    -----
    Box dimensions are assumed to be orthorhombic
    GOMC input files utilize Angstroms, which are also the units of the parmed Structure box
    """
    if pressure is None:
        gemc_type = "GEMC   NVT"
        move_set = """DisFreq               0.60
RotFreq           0.10
RegrowthFreq          0.30\n"""

    else:
        gemc_type = "GEMC   NPT\nPressure {}\n".format(pressure.in_units(
            u.bar).value)
        move_set = """DisFreq               0.59
RotFreq           0.10
RegrowthFreq          0.30
VolFreq             0.01\n"""


    with open(filename, 'w') as f:
        f.write("""
###########################################################################
#  ========-------------------- INPUT --------------------------===========
############################################################################

#########################
# enable, step
#########################
Restart     false


####################################
# kind (RESTART, RANDOM, INTSEED)
####################################
PRNG        RANDOM

####################################
# FORCE FIELD
####################################
ParaTypeCHARMM   true
ParaTypeEXOTIC   false
Parameters      {parameters}

####################################
# INPUT PDB FILES
####################################
Coordinates 0    {coords}

####################################
# INPUT PSF FILES
####################################
Structure 0      {structure}



############################################################################
#  =======--------------------- SYSTEM --------------------------===========
############################################################################

##################################
# GEMC TYPE (DEFULT IS NVT_GEMC)
##################################
{gemc_type}

#############################
# SIMULATION CONDITION
#############################
Temperature     {temperature}
Potential       VDW
LRC     true
Rcut        14
Exclude     1-4

#############################
# ELECTROSTATIC
#############################
ElectroStatic   false
Ewald           false
PressureCalc  true  1000

################################
# STEPS
################################
RunSteps       {n_steps}
EqSteps        500000
AdjSteps       1000

################################
# MOVE FREQUENCY
################################
{move_set}

################################
# BOX DIMENSION #, X, Y, Z
################################
CellBasisVector1 0  {box_x}  0.00    0.00
CellBasisVector2 0  0.00    {box_y}  0.00
CellBasisVector3 0  0.00    0.00    {box_z}


##############################
# CBMC TRIALS
##############################
CBMC_First   10
CBMC_Nth     8
CBMC_Ang     100
CBMC_Dih     20

####################################
#          Mol. Name     Chem. Pot.
####################################
############################################################################
#  =======-------------------- OUTPUT --------------------------===========
############################################################################

##########################
# statistics filename add
##########################
OutputName  {output}

#####################################
# enable, frequency
#####################################
CoordinatesFreq    true   1000000
RestartFreq        true   1000000
ConsoleFreq        true   100000
BlockAverageFreq   true   100000
HistogramFreq      false  100000


################################
# OutHistSettings
################################


##################################
# enable: blk avg., fluct.
##################################
OutEnergy         true    true
OutPressure       true    true
OutMolNum         true    true
OutDensity        true    true
""".format(structure=structure, coords=coords, gemc_type=gemc_type,
        parameters=parameters, output=output, temperature=temperature.value,
        box_x=parm_structure.box[0], box_y=parm_structure.box[1], 
        box_z=parm_structure.box[2], n_steps=n_steps, move_set=move_set))
