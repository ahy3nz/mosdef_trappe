# MoSDeF + TraPPE
Utilizing the Molecular Simulation Design Framework with the TraPPE force field
to build and parametrize systems suitable for a variety of simulation engines

## TraPPE implementation details
* 1.4 nm cutoff
* Tail corrections to LJ interactions
* Ewald summation for Coulombic interactions
* 1-4 LJ and Coulombic interactions excluded
    * 1-2 and 1-3 LJ and Coulombic interactions also excluded
* Fixed bond lengths
* 2 fs timestep
* Lorentz-Berthelot mixing rules
    * Arithmetic mean for sigma
    * Geometric mean for epsilon

## Helpful Links
* [NIST reference examples](https://www.nist.gov/mml/csd/chemical-informatics-research-group/lammps-md-equation-state-pressure-vs-density-trappe-0)
