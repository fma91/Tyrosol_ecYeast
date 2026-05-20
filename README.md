# Tyrosol_ecYeast

Genome-scale **strain design for tyrosol production** in *Saccharomyces
cerevisiae* using an enzyme-constrained model (`ecTyrosol.mat`) and the
[ecFactory](https://github.com/SysBioChalmers/ecFactory) pipeline.

Repository layout for model building, simulation, and target prediction
outputs (L1, L2, L3, transporters).

## Contents

```
Tyrosol_ecYeast/
├── README.md
├── model/
│   ├── ecTyrosol.mat              enzyme-constrained tyrosol production model
│   └── build_ecTyrosol_model.py   builds ecTyrosol.mat from ecYeastGEM_batch
├── scripts/
│   └── run_tyrosol_ecFactory.m    ecFactory strain-design pipeline (MATLAB)
├── results/
│   ├── candidates_L1.txt
│   ├── candidates_L2.txt
│   ├── candidates_L3.txt
│   └── transporter_targets.txt
└── docs/
    └── METHODS.md                 pathway assumptions and simulation settings
```

## Workflow

1. **Model construction** — Extend `ecYeastGEM_batch.mat` with the Ehrlich
   tyrosol pathway (ARO10, ADH7) and product export reactions. Modeling
   assumptions are listed in `build_ecTyrosol_model.py` and `docs/METHODS.md`.

2. **Strain design** — Run ecFactory on minimal glucose medium to obtain
   ranked gene targets (OE / KD / KO) and transporter candidates.

## Model assumptions (summary)

- Tyrosol via native Ehrlich pathway: 4-HPP → 4-HPAA (ARO10) → tyrosol (ADH7).
- Enzyme-constrained arms on the new steps with kcat = 1000 s⁻¹.
- Wild-type yeast background; minimal medium for the simulation.

## ecFactory settings

| Parameter | Value |
|---|---|
| Medium | Minimal, D-glucose |
| WT biomass yield | 0.48 g / g glucose |
| Target yield in scan | 0.49 × WT yield |

## Requirements

- Python 3 + `scipy` (model build)
- MATLAB + RAVEN + Gurobi + GECKO 2.0.3 (`~/Documents/ecFactory/code/GECKO`)

## Quick start

```bash
git clone <repository-url>
cd Tyrosol_ecYeast

# 1) Build the ecModel (requires ecYeastGEM_batch.mat from CellFactory-ecYeastGEM)
pip install -r requirements.txt
python model/build_ecTyrosol_model.py

# 2) Run ecFactory in MATLAB (GECKO 2.0.3 + RAVEN + Gurobi on your path)
matlab -nodisplay -batch "addpath('scripts'); run_tyrosol_ecFactory"
```

Outputs are written to `results/`. See `docs/METHODS.md` for assumptions and
parameters.

## External dependencies (not included in this repo)

| Component | Location / note |
|---|---|
| Base ecModel | `CellFactory-ecYeastGEM/ModelFiles/ecYeastGEM_batch.mat` |
| ecFactory + GECKO 2.0.3 | `~/Documents/ecFactory/code` with `GECKO` → GECKO 2.0.3 |
| RAVEN Toolbox | User MATLAB installation |
| Gurobi | LP solver for RAVEN |

## Citation

If you use this workflow, cite ecFactory, GECKO, and ecYeastGEM as appropriate
for your manuscript (add your publication reference here).
