# Methods — Tyrosol production model and ecFactory strain design

## Objective

Identify metabolic engineering targets in *Saccharomyces cerevisiae* to improve
tyrosol production using an enzyme-constrained genome-scale model (ecModel)
and the ecFactory pipeline (flux-scanning with enforced objective function,
FSEOF, followed by enzyme usage variability analysis).

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

Three metabolites and four reactions are added to the batch model
(`build_ecTyrosol_model.py`). The product exchange reaction is set as the
model objective (`model.c`).

### Modeling assumptions

| ID | Assumption |
|---|---|
| A1 | Tyrosol is produced via the native Ehrlich route from tyrosine, not via a heterologous tyrosine-to-tyramine bypass. |
| A2 | ARO10 catalyses decarboxylation of 4-HPP to 4-HPAA (reported substrate promiscuity of Aro10p). |
| A3 | ADH7 catalyses reduction of 4-HPAA to tyrosol (NADPH-dependent aldehyde reductase). |
| A4 | Dedicated enzyme arms for the new steps use **kcat = 1000 s⁻¹** in the GECKO representation (S-coefficient = 1 / (kcat × 3600)). This parametrizes a high-capacity production route in the enzyme-constrained framework. |
| A5 | Wild-type genetic background (no pre-applied gene deletions). |
| A6 | Product is exported through a transport + exchange pair; exchange flux is the optimization objective during model construction. |

## Strain design simulation (ecFactory)

**Software:** ecFactory (GECKO 2.0.3), RAVEN Toolbox, Gurobi.

**Medium:** minimal medium, D-glucose as sole carbon source.

**Yield scan:** `WT_yield = 0.48` g biomass / g glucose;  
`expYield = 0.49 × WT_yield` (suboptimal biomass fraction used in the FSEOF scan).

**Procedure:** `run_tyrosol_ecFactory.m` loads `ecTyrosol.mat`, applies the
medium, sets the growth reaction as the active objective for FSEOF, and writes
filtered target lists to `results/`:

- `candidates_L1.txt` — initial FSEOF targets
- `candidates_L2.txt` — after removing essential genes
- `candidates_L3.txt` — after EUVA filtering and optimal target combination
- `transporter_targets.txt` — transport reactions without gene association

## Reproducibility

```bash
cd ~/Documents/Tyrosol_ecYeast/model
python build_ecTyrosol_model.py

matlab -nodisplay -batch "addpath('~/Documents/Tyrosol_ecYeast/scripts'); run_tyrosol_ecFactory"
```

Prerequisites: Python with `scipy`; MATLAB with RAVEN, Gurobi, and GECKO 2.0.3
linked from `~/Documents/ecFactory/code/GECKO`.
