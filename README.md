# Tyrosol_ecYeast

Genome-scale **strain design for tyrosol production** in *Saccharomyces
cerevisiae* using an enzyme-constrained model (`ecTyrosol.mat`) and the
[ecFactory](https://github.com/SysBioChalmers/ecFactory) pipeline (GECKO 2.0.3
+ RAVEN).

## Contents

```
Tyrosol_ecYeast/
├── model/
│   ├── build_ecTyrosol_model_raven.m   build ecTyrosol.mat (RAVEN)
│   └── ecTyrosol.mat                   enzyme-constrained tyrosol model
├── scripts/
│   └── run_tyrosol_ecFactory.m         ecFactory strain design
├── results/                            target lists (L1, L2, L3, transporters)
└── docs/
    └── METHODS.md                      assumptions and parameters
```

## Reproducibility (RAVEN + GECKO)

```matlab
cd Tyrosol_ecYeast/model
build_ecTyrosol_model_raven          % ecTyrosol.mat from ecYeastGEM_batch.mat

addpath(fullfile(pwd, '..', 'scripts'))
run_tyrosol_ecFactory                % results/*.txt
```

## Workflow

1. **Model** — `build_ecTyrosol_model_raven.m` extends `ecYeastGEM_batch.mat`
   with the Ehrlich pathway (ARO10, ADH7), export reactions, and product
   objective using RAVEN `addMets` / `addRxns` / `setParam`.

2. **Strain design** — `run_tyrosol_ecFactory.m` runs ecFactory on minimal
   glucose medium and writes ranked targets to `results/`.

## Model assumptions

- Native Ehrlich route: 4-HPP → 4-HPAA (ARO10) → tyrosol (ADH7).
- Enzyme-constrained new steps: kcat = 1000 s⁻¹ (GECKO protein coefficients).
- Wild-type background; minimal medium in simulation.

## ecFactory settings

| Parameter | Value |
|---|---|
| Medium | Minimal, D-glucose |
| WT biomass yield | 0.48 g / g glucose |
| Target yield in scan | 0.49 × WT yield |

## Requirements

- MATLAB
- [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)
- [GECKO 2.0.3](https://github.com/SysBioChalmers/GECKO) and [ecFactory](https://github.com/SysBioChalmers/ecFactory)
- Gurobi
- `ecYeastGEM_batch.mat` from [CellFactory-ecYeastGEM](https://github.com/SysBioChalmers/CellFactory-ecYeastGEM)

See `docs/METHODS.md` for full assumptions (A1–A6) and paths.

## Citation

Cite ecFactory, GECKO, ecYeastGEM, and this work as appropriate for your manuscript.
