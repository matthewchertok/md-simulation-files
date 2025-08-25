# 2025-mol-containers

This repository contains all files required to recreate the simulations and analysis present in "Molecular Cage–Directed Hierarchical Precipitation for Selective Recovery of Aliphatic and Perfluoro Surfactants" by Pérez-Ferreiro *et al*.

### MD Simulations

MD simulations require the LAMMPS simulation package and were executed using the stable release from June 23, 2022. The files required to recreate simulations for p-A4B4, p-A2C3, and p-A are available in directories `a4b4/`, `a2c3/`, and `a/`, respectively. The `sys.data`, `sys.settings` and `tip4p.mol` files are used for system preparation. `equil.in` prepares the system and performs minimization and equilibration. `prod.in` performs production simulations used for analysis.

### Aggregation Analysis

All files and data related to the aggregation analysis are available in the `analysis/` directory. `compute_aggregates.py` implements the aggregation calculations, which requires LAMMPS trajectory files for execution. The `analysis/data/` directory contains the outputs of this script applied to the trajectories generated using the files in `a4b4/`, `a2c3/`, and `a/`.
