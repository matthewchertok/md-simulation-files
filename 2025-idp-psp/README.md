# EXAMPLE LAMMPS SCRIPTS

This directory contains example LAMMPS simulation input scripts for the paper [When B<sub>2</sub> is Not Enough: Evaluating Simple Metrics for Predicting Phase Separation of Intrinsically Disordered Proteins](https://doi.org/10.48550/arXiv.2507.12312).
Scripts for three distinct types of simulations are contained within the [example_lammps_scripts](/2025-idp-psp/example_lammps_scripts) directory: single chain simulations for calculating the radius of gyration, adaptive biasing force simulations for calculating the second virial coefficient, and 100 chain simulations used to compute an approximate pressure-density equation-of-state.
All files are for simulations of an example IDP given by the amino acid sequence EQEFSDNELQELSTQGSRYV.
View the [list of all sequences](https://g-ef94ef.f0ad1.36fe.data.globus.org/10.34770/6tnm-7b56/390/seq_heteromeric.txt) used in the study.
Force field parameters contained within the sys.settings files are based off of work by [Regy et al.](https://doi.org/10.1002/pro.4094)


### Single chain simulations

The single chain simulation scripts can be found within the [rg_scripts](/2025-idp-psp/example_lammps_scripts/rg_scripts) directory contained within example_lammps_scripts. 
The files are:

- **start.lmp** --- This file contains general LAMMPS instructions and utilizes information from sys.data and sys.settings
- **sys.data** --- This file contains the simulation starting configuration
- **sys.settings** --- This file contains the force field settings including non-bonded and bonded interactions


### Adaptive biasing force simulations

The adaptive biasing force simulation scripts can be found within the [b2_scripts](/2025-idp-psp/example_lammps_scripts/b2_scripts) directory contained within example_lammps_scripts. 
The files are:

- **colvars.inp** -- This file contains instructions to run the adaptive biasing force portion of the simulation
- **start.lmp** --- This file contains general LAMMPS instructions and utilizes information from colvars.inp, sys.data, and sys.settings
- **sys.data** --- This file contains the simulation starting configuration
- **sys.settings** --- This file contains the force field settings including non-bonded and bonded interactions


### Equation-of-state simulations

The equation-of-state simulation scripts can be found within the [eos_scripts](/2025-idp-psp/example_lammps_scripts/eos_scripts) directory contained within example_lammps_scripts. 
The files are:

- **start.lmp** --- This file contains general LAMMPS instructions and utilizes information from sys.data and sys.settings
- **sys.data** --- This file contains the simulation starting configuration. Note that the box dimensions for this simulation are set to reflect a density of 0.05 g/mL
- **sys.settings** --- This file contains the force field settings including non-bonded and bonded interactions
