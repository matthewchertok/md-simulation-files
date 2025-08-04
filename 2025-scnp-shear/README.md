# Reproducing the Simulations

**Paper:** *Effect of Shear Flow and Precursor Polymer Design on Single‑Chain Nanoparticle Formation*

This repository provides all scripts required to replicate the molecular‑dynamics (MD) simulations reported in the paper cited above.

---

## 1. Workflow Overview

### Step 1 – Generate Relaxed Precursor Chains  
Submit the job array:

```bash
sbatch job_gen_seed_conformations.slurm
```

`main_generate_seed_conformations.py` invokes `generate_seed_conformations.py` to create **10** equilibrated conformations for every unique combination of the parameters below.

| Parameter | Values |
|-----------|--------|
| Linker fraction *f* | 0.1, 0.2, 0.4 |
| Linker blockiness *β* | 0.2, 0.8 |
| Backbone stiffness *κ* | 0, 5, 10 |
| Plate speed *v<sub>x</sub>* | 0, 0.05, 0.25, 0.5 |
| Chain length *N* | 50, 100, 150 |

*Plate separation* = `max(100, N)`  
*Shear rate* = \( \dot{\gamma}= 2 v_x /(\text{plate separation}) \) (plates move at ± *v<sub>x</sub>*).

---

### Step 2 – Run Cross‑Linking Production Simulations  
Submit the second job array:

```bash
sbatch job_production_irreversible.slurm
```

`main_generate_production_data.py` reads every seed conformation (10 per parameter set) and calls `run_reactive_steps.py`.  
If your scheduler limits array size, adjust `START_IDX` and `END_IDX` inside the script so the slice of `file_names` aligns with the array‑task indices.

All trajectories are saved as `.gsd` files.

---

## 2. Script Catalogue

### 2.1 Top‑Level Drivers
| Script | Description |
|--------|-------------|
| `main_generate_seed_conformations.py` | Generates relaxed precursor chains. |
| `main_generate_production_data.py` | Launches reactive production runs. |
| `generate_seed_conformations.py` | Helper routine called by the driver above. |
| `run_reactive_steps.py` | Executes reactive MD with on‑the‑fly bond formation. |

### 2.2 `sim_setup/` Utilities
| File | Purpose |
|------|---------|
| `parse_args.py` | Centralised parameter definitions (e.g., `init_conformations_folder`, `production_save_folder`). |
| `place_md_particles.py` | Builds SCNP backbone & linker beads; defines bond, angle, and excluded‑volume potentials. |
| `place_mpcd_particles.py` | Adds MPCD solvent particles to the box. |
| `create_simulation.py` | Assembles the full system, including parallel plates and integrator. |
| `make_bonds_custom_action.py` | Detects eligible cross‑links during the run and updates the topology. |

---

## 3. Reproducing the Results

Running only the two Slurm submission scripts listed above (`job_gen_seed_conformations.slurm` and `job_production_irreversible.slurm`) will generate every dataset analysed in the manuscript.

---