import unyt as u 


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
