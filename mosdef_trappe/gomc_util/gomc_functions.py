import parmed
import subprocess
import unyt as u


def simulate(parametrized_structure, **kwargs):
    """ Simulate using GOMC"""
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

    # Write GOMC run file, which requires box information and 
    # knowledge of other input files
    write_gomc_nvt_input('in.conf', parametrized_structure, 
            coords='coords.pdb', structure='structure.psf', 
            parameters='parameters.par', **kwargs)

    # Post-process PAR file
    modify_par_file('parameters.par')

    run_nvt()

def run_nvt():
    p = subprocess.Popen('GOMC_CPU_NVT in.conf', shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True)
    out, err = p.communicate()
    with open("gomc.out", 'w') as f:
        f.write(out)
    with open("gomc.err", 'w') as f:
        f.write(err)

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


def write_gomc_gemc_input(filename, param_structure0, param_structure1, 
                            temperature=273*u.Kelvin, parameters='parameters.par',
                            coords0='coords0.pdb', structure0='structure0.psf', 
                            coords1='coords1.pdb', structure1='structure1.psf',
                            output='output'):
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
    temperature : float
        Implicitly in Kelvin
    parameters : str
        CHARMM/NAMD-style PAR file for force field information
    coords0 : str
        PDB file for param_structure0 coordinates and atom names
    structure0: str
        PSF file for param_structure0 atomtypes, bonding information, and topology
    coords1 : str
        PDB file for param_structure0 coordinates and atom names
    structure1: str
        PSF file for param_structure0 atomtypes, bonding information, and topology    
    output : str
        Prefix for output files from simulation
        
    Notes
    -----
    Box dimensions are assumed to be orthorhombic
    GOMC input files utilize Angstroms, which are also the units of the parmed Structure box
    TraPPE uses 14 Angstrom cutoffs, but here we use 10 Angstrom because our simulations are small
    """
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
GEMC      NVT

#############################
# SIMULATION CONDITION
#############################
Temperature     {temperature}
Potential       VDW
LRC     true
Rcut        10
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
RunSteps       1000000
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
""".format(temperature=temperature.value, parameters=parameters,
           structure0=structure0, coords0=coords0, 
           structure1=structure1, coords1=coords1,
           box0_x=param_structure0.box[0], box0_y=param_structure0.box[1], 
           box0_z=param_structure0.box[2],
           box1_x=param_structure1.box[0], box1_y=param_structure1.box[1], 
           box1_z=param_structure1.box[2],
           output=output))

def write_gomc_nvt_input(filename, parm_structure, coords='coords.pdb',
                     structure='structure.psf', parameters='parameters.par',
                    output='output', temperature=305*u.Kelvin):
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


#############################
# SIMULATION CONDITION
#############################
Temperature     {temperature}
Potential       VDW
LRC     true
Rcut        10
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
RunSteps       1000000
EqSteps        500000
AdjSteps       1000

################################
# MOVE FREQUENCY
################################
DisFreq               0.60
RotFreq           0.10
RegrowthFreq          0.30


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
""".format(structure=structure, coords=coords, 
        parameters=parameters, output=output, temperature=temperature.value,
        box_x=parm_structure.box[0], box_y=parm_structure.box[1], 
        box_z=parm_structure.box[2]))
