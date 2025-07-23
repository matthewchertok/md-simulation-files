# README

This directory contains simulation and analysis data for single-chain nanoparticle (SCNP) conformations generated under shear and their corresponding quiescent relaxation states. Detailed descriptions of each dataset and instructions for merging auxiliary files are provided below.

---

## 1. Selected Conformations Under Shear

**Location:** `selected_conformations_under_shear`

These files are pickled pandas DataFrames capturing ten independent SCNP snapshots from a 20 million–step production run under shear. The final 10 million steps were sampled at 1 million–step intervals. Each DataFrame includes:

- **filename**  
  Name of the source GSD file, e.g.  
  `production_combo0_chain0_f=0.1,beta=0.2,kappa=0,plate_speed=constant=0.0_bath_temp=1.0_box_dim=100_n_backbone=50_n_sidechain=5_replicate0.gsd`
- **anisotropy**  
  Relative shape anisotropy.
- **rg**  
  Radius of gyration.
- **lc**  
  Contour length.
- **rg_over_lc**  
  Radius of gyration divided by contour length.
- **frame_index**  
  Snapshot index (0–9).
- **pdist_hist_normalized_density**  
  Histogram normalized so that ∑(count × bin width) = 1.
- **f**  
  Linker fraction.
- **beta**  
  Blockiness parameter of linkers.
- **kappa**  
  Backbone stiffness.
- **shear_rate**  
  Shear rate during cross-linking.
- **combo_chain_id**  
  Identifier encoding simulation parameters and replicate.
- **temperature**
  Thermal energy (k_BT), always 1.0.
- **n_back**  
  Number of backbone monomers.
- **relaxation_time**  
  Calculated as  
  `0.0513 × n_backbone × (kappa_over_kbt × n_backbone + kappa_over_kbt + n_backbone)`.
- **kappa_over_kbt**  
  Backbone stiffness divided by temperature.
- **weissenberg**  
  Weissenberg number (shear_rate × relaxation_time).

---

## 2. Selected Conformations Under Quiescent Conditions

**Location:** `selected_conformations_quiescent`

Structure and content mirror the “Under Shear” dataset, but snapshots follow a 2 million–step relaxation (cross‑linking disabled) and a subsequent 10 million–step production run at zero shear. The `filename` field uses the suffix  
`_crosslinked_equilibrium.gsd`.

---

## 3. UMAP Weights

**Location:** `umap_weights`

Contains precomputed weights for the UMAP reducer used in embedding generation. Due to file size, the weights have been split into multiple parts. To reconstruct the original file:

```bash
cd umap_weights
cat umap_weights.sav.part_* > umap_weights.sav
```

Load the merged weights in Python via Joblib:

```python
import joblib
umap_model = joblib.load('umap_weights/umap_weights.sav')
```

---

## 4. Topological Domain Data

**File:** `topological_domain_data.pkl`

This pickle stores SCNP morphological metrics extracted from the first frame of each quiescent relaxation simulation:

- **n_topological_domains**  
  Number of distinct topological domains after cross-linking.
- **topological_domain_sizes**  
  List of domain sizes.
- **n_free_backbone**  
  Number of backbone beads not assigned to any domain.
- **f**, **kappa_over_kbt**, **n_back**, **beta**  
  Linker fraction, stiffness parameter, backbone monomer count, and blockiness, respectively.
- **weissenberg**  
  Weissenberg number.
- **median_domain_size**  
  Median size of all topological domains.