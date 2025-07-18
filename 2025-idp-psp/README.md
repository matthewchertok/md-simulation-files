### EXAMPLE LAMMPS SCRIPTS

## Summary

This directory contains example LAMMPS simulation input scripts for the paper [When B<sub>2</sub> is Not Enough: Evaluating Simple Metrics for Predicting Phase Separation of Intrinsically Disordered Proteins](https://doi.org/10.48550/arXiv.2507.12312).
Scripts for three distinct types of simulations are contained within [example_lammps_scripts](/2025-idp-psp/example_lammps_scripts): single chain simulations for calculating the radius of gyration, adaptive biasing force simulations for calculating the second virial coefficient, and 100 chain simulations used to compute an approximate pressure-density equation-of-state.
All files are for simulations of a single example IDP sequence, EQEFSDNELQELSTQGSRYV.
View the [list of all sequences](https://g-ef94ef.f0ad1.36fe.data.globus.org/10.34770/6tnm-7b56/390/seq_heteromeric.txt).

# Single chain simulations

The single chain scripts can be found within the [rg_scripts](/2025-idp-psp/example_lammps_scripts/rg_scripts) directory contained within example_lammps_scripts. 
Several files can be found within:

- **start.lmp** --- This file contains general LAMMPS instructions
- **sys.data** --- This file contains the simulation starting configuration
- **sys.settings** --- This file contains the force field settings
