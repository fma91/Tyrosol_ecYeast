# Methods — Tyrosol production model and ecFactory strain design

## Objective

Identify metabolic engineering targets in *Saccharomyces cerevisiae* to improve
tyrosol production using an enzyme-constrained genome-scale model (ecModel)
and the ecFactory pipeline (flux-scanning with enforced objective function,
FSEOF, followed by enzyme usage variability analysis).

## Software

| Step | Software |
|---|---|
| Model construction | **MATLAB**, **RAVEN Toolbox** (`addMets`, `addRxns`, `setParam`) |
| Strain design | **MATLAB**, **RAVEN**, **GECKO 2.0.3**, **Gurobi**, **ecFactory** |

Models are built with RAVEN model-editing functions only — not by direct
hand-editing of `.mat` files.

## Base model

- **Source:** `ecYeastGEM_batch.mat` from CellFactory-ecYeastGEM.
- **Background:** wild-type yeast (no chassis gene deletions in this build).

## Pathway and model extensions

Tyrosol production is represented through the **Ehrlich pathway** from
L-tyrosine:

1. L-tyrosine → 4-hydroxyphenylpyruvate (4-HPP) — native yeast reactions.
2. 4-HPP → 4-hydroxyphenylacetaldehyde (4-HPAA) — **ARO10** (YDR380W).
3. 4-HPAA → tyrosol — **ADH7** (YCR105W, NADPH-dependent).
4. Tyrosol transport (cytoplasm → extracellular) and tyrosol exchange.

`model/build_ecTyrosol_model_raven.m` adds three metabolites and four
reactions and sets `new_tyrosol_ex` as the product objective. Outputs:
`ecTyrosol.mat` and `ecTyrosol_native.mat` (alias).

### Modeling assumptions

| ID | Assumption |
|---|---|
| A1 | Tyrosol is produced via the native Ehrlich route from tyrosine, not via a heterologous tyrosine-to-tyramine bypass. |
| A2 | ARO10 catalyses decarboxylation of 4-HPP to 4-HPAA. |
| A3 | ADH7 catalyses reduction of 4-HPAA to tyrosol (NADPH-dependent). |
| A4 | Dedicated enzyme arms for the new steps use **kcat = 1000 s⁻¹** (S-coefficient = 1 / (kcat × 3600)). |
| A5 | Wild-type genetic background (no pre-applied gene deletions). |
| A6 | Product is exported through a transport + exchange pair; exchange flux is the optimization objective during model construction. |

## Strain design simulation (ecFactory)

**Medium:** minimal medium, D-glucose as sole carbon source.

**Yield scan:** `WT_yield = 0.48` g biomass / g glucose;  
`expYield = 0.49 × WT_yield`.

**Procedure:** `scripts/run_ecTyrosol_native.m` loads `ecTyrosol_native.mat`,
applies the medium via `changeMedia_batch`, switches `model.c` to the biomass
reaction for FSEOF, and runs ecFactory. Outputs in `results/`:

- `candidates_L1.txt` — initial FSEOF targets
- `candidates_L2.txt` — after removing essential genes
- `candidates_L3.txt` — after EUVA filtering
- `transporter_targets.txt` — transport reactions without gene association

## Reproducibility

```matlab
cd Tyrosol_ecYeast/model
build_ecTyrosol_model_raven

addpath('../scripts')
run_ecTyrosol_native
```

**External dependencies (not in this repo):**

| Component | Location |
|---|---|
| Base ecModel | `CellFactory-ecYeastGEM/ModelFiles/ecYeastGEM_batch.mat` |
| ecFactory + GECKO 2.0.3 | `~/Documents/ecFactory/code/GECKO` |
| RAVEN Toolbox | User MATLAB installation |
| Gurobi | LP solver for RAVEN |

`ecTyrosol.mat` is an alias of `ecTyrosol_native.mat` for envelope and dFBA
scripts in this repository.
