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
| Strain design | **MATLAB**, **RAVEN Toolbox**, **GECKO 2.0.3**, **Gurobi**, **ecFactory** |

Reviewers can reproduce model construction and target prediction from
`ecYeastGEM_batch.mat` using only the scripts in this repository and the
external dependencies listed below.

## Base model

- **Source:** `ecYeastGEM_batch.mat` from CellFactory-ecYeastGEM (enzyme-constrained
  yeast GEM).
- **Background:** wild-type yeast (no chassis gene deletions in this build).

## Pathway and model extensions

Tyrosol production is represented through the **Ehrlich pathway** from
L-tyrosine:

1. L-tyrosine → 4-hydroxyphenylpyruvate (4-HPP) — native yeast reactions.
2. 4-HPP → 4-hydroxyphenylacetaldehyde (4-HPAA) — **ARO10** (YDR380W).
3. 4-HPAA → tyrosol — **ADH7** (YCR105W, NADPH-dependent).
4. Tyrosol transport (cytoplasm → extracellular) and tyrosol exchange.

`model/build_ecTyrosol_model_raven.m` adds three metabolites and four
reactions to the batch model and sets `new_tyrosol_ex` as the product
objective.

### Modeling assumptions

| ID | Assumption |
|---|---|
| A1 | Tyrosol is produced via the native Ehrlich route from tyrosine, not via a heterologous tyrosine-to-tyramine bypass. |
| A2 | ARO10 catalyses decarboxylation of 4-HPP to 4-HPAA (reported substrate promiscuity of Aro10p). |
| A3 | ADH7 catalyses reduction of 4-HPAA to tyrosol (NADPH-dependent aldehyde reductase). |
| A4 | Dedicated enzyme arms for the new steps use **kcat = 1000 s⁻¹** in the GECKO representation (S-coefficient = 1 / (kcat × 3600)). |
| A5 | Wild-type genetic background (no pre-applied gene deletions). |
| A6 | Product is exported through a transport + exchange pair; exchange flux is the optimization objective during model construction. |

## Strain design simulation (ecFactory)

**Medium:** minimal medium, D-glucose as sole carbon source.

**Yield scan:** `WT_yield = 0.48` g biomass / g glucose;  
`expYield = 0.49 × WT_yield`.

**Procedure:** `scripts/run_tyrosol_ecFactory.m` loads `ecTyrosol.mat`, applies
the medium via GECKO (`changeMedia_batch`), and runs ecFactory. Outputs in
`results/`:

- `candidates_L1.txt` — initial FSEOF targets
- `candidates_L2.txt` — after removing essential genes
- `candidates_L3.txt` — after EUVA filtering
- `transporter_targets.txt` — transport reactions without gene association

## Reproducibility

```matlab
% 1) Build ecTyrosol.mat (RAVEN)
cd('~/Documents/Tyrosol_ecYeast/model')
build_ecTyrosol_model_raven

% 2) Strain design (RAVEN + GECKO 2.0.3 + ecFactory)
addpath('~/Documents/Tyrosol_ecYeast/scripts')
run_tyrosol_ecFactory
```

**External dependencies (not in this repo):**

- `CellFactory-ecYeastGEM/ModelFiles/ecYeastGEM_batch.mat`
- [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)
- [GECKO 2.0.3](https://github.com/SysBioChalmers/GECKO) on the MATLAB path
  (e.g. via ecFactory: `~/Documents/ecFactory/code/GECKO`)
- [ecFactory](https://github.com/SysBioChalmers/ecFactory)
- Gurobi (LP solver for RAVEN)

Set `ECYEASTGEM_BATCH` if the base model is not at the default path above.
