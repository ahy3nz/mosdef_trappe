# MoSDeF + TraPPE
Utilizing the Molecular Simulation Design Framework with the TraPPE force field
to build and parametrize systems suitable for a variety of simulation engines

## Examples
To run the examples from this git repo, `cd` into `mosdef_trappe/mosdef_trappe`.
The examples should be run from the same directory as the examples themselves. 

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
* [NIST standard reference simulation website](https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website)
