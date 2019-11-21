import unyt as u 
import subprocess

def simulate(parametrized_structure, **kwargs):
    """ Umbrella routine to simulate in gromacs 

    Parameters
    ----------
    kwargs : keyword argumnets passed to write_mdp
    """
    parametrized_structure.save('compound.gro', overwrite=True)
    parametrized_structure.save('compound.top', overwrite=True)

    write_em_mdp('em.mdp')
    grompp = run_grompp(mdp='em.mdp',
            gro='compound.gro', top='compound.top', out='em')
    if grompp.returncode == 0:
        mdrun = run_md('em')

    write_mdp('md.mdp', **kwargs)
    grompp = run_grompp(mdp='md.mdp',
            gro='em.gro', top='compound.top', out='md')
    if grompp.returncode == 0:
        mdrun = run_md('md')


def run_grompp(mdp='md.mdp', gro='compound.gro', top='compound.top',
        out='out'):
    """ Run GMX grompp """
    to_run = 'gmx grompp -f {mdp} -c {gro} -p {top} -o {out} -maxwarn 1'.format(
            mdp=mdp, gro=gro, top=top, out=out)
    grompp = subprocess.Popen(to_run, shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
    out, err = grompp.communicate()
    with open('grompp.out', 'w') as f:
        f.write(out)
    with open('grompp.err', 'w') as f:
        f.write(err)

    return grompp


def run_md(output):
    """Run gmx mdrun """
    mdrun = subprocess.Popen('gmx mdrun -deffnm {}'.format(output),
                shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True)
    out, err = mdrun.communicate()
    with open('mdrun.out', 'w') as f:
        f.write(out)
    with open('mdrun.err', 'w') as f:
        f.write(err)

    return mdrun


def write_mdp(filename, pressure=None, **kwargs):
    if pressure is None:
        write_nvt_mdp(filename, **kwargs)
    else:
        write_npt_mdp(filename, pressure=pressure, **kwargs)


def write_npt_mdp(filename, temperature=300*u.Kelvin,
        pressure=1*u.atm, random_seed=42,
        n_steps=500000):
    """ Some args are made accessible since they should be changeable,
    others are "hardcoded" as they probably shouldn't change unless
    absolutely necessary """

    temperature.convert_to_units(u.Kelvin)
    pressure.convert_to_units(u.bar)

    with open(filename, 'w') as f:
        f.write("""
title                     = NPT
; Run parameters
integrator                = md        ; leap-frog integrator
nsteps                    = {n_steps}     
dt                        = 0.002     ; 2 fs

; Output control
nstxout                   = 5000       ; Every 10 ps
nstvout                   = 5000
nstenergy                 = 5000
nstlog                    = 5000

;Bond parameters
continuation              = no
constraint_algorithm    = lincs
constraints             =   all-bonds
lincs_iter              = 1
lincs_order             = 4

; Neighbor searching
cutoff-scheme           = Verlet
nstype                  = grid
nstlist                 = 10
rcoulomb                = 1.4
rvdw                    = 1.4

; Electrostatics
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

; Temperature coupling
tcoupl                  = nose-hoover
tc-grps                 = system
tau_t                   = 0.4  
ref_t                   = {temperature}  

; Pressure coupling
pcoupl                      = Parrinello-Rahman
pcoupltype                  = isotropic
tau_p                       = 2.0           ; ps
ref_p                       = {pressure}          ; bar 
compressibility             = 4.5e-5
refcoord_scaling            = com

;Periodic boundary conditions
pbc                     = xyz

;Dispersion correction
DispCorr                = EnerPres

;Velocity generation
gen_vel                 = yes
gen_temp                = {temperature}
gen_seed                = {random_seed}""".format(
        temperature=temperature.value,
        pressure=pressure.value,
        n_steps=n_steps,
        random_seed=random_seed))

def write_nvt_mdp(filename, temperature=300*u.Kelvin,
        random_seed=42,
        n_steps=500000):
    """ Some args are made accessible since they should be changeable,
    others are "hardcoded" as they probably shouldn't change unless
    absolutely necessary """

    temperature.convert_to_units(u.Kelvin)

    with open(filename, 'w') as f:
        f.write("""
title                     = NVT
; Run parameters
integrator                = md        ; leap-frog integrator
nsteps                    = {n_steps}     
dt                        = 0.002     ; 2 fs

; Output control
nstxout                   = 5000       ; Every 10 ps
nstvout                   = 5000
nstenergy                 = 5000
nstlog                    = 5000

;Bond parameters
continuation              = no
constraint_algorithm    = lincs
constraints             =   all-bonds
lincs_iter              = 1
lincs_order             = 4

; Neighbor searching
cutoff-scheme           = Verlet
nstype                  = grid
nstlist                 = 10
rcoulomb                = 1.4
rvdw                    = 1.4

; Electrostatics
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

; Temperature coupling
tcoupl                  = nose-hoover
tc-grps                 = system
tau_t                   = 0.4  
ref_t                   = {temperature}  

; Pressure coupling
pcoupl                      = no

;Periodic boundary conditions
pbc                     = xyz

;Dispersion correction
DispCorr                = EnerPres

;Velocity generation
gen_vel                 = yes
gen_temp                = {temperature}
gen_seed                = {random_seed}""".format(
        temperature=temperature.value,
        n_steps=n_steps,
        random_seed=random_seed))


def write_em_mdp(filename):
    with open(filename, 'w') as f:
        f.write("""
integrator  = steep     ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0    ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps      = 50000     ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions (yes/no)
""")
