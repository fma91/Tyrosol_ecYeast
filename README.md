# Tyrosol_ecYeast

Genome-scale **strain design for tyrosol production** in *Saccharomyces
cerevisiae* using an enzyme-constrained model (`ecTyrosol_native.mat`) and the
[ecFactory](https://github.com/SysBioChalmers/ecFactory) pipeline (GECKO 2.0.3
+ RAVEN).

## Contents

```
Tyrosol_ecYeast/
├── model/
│   ├── build_ecTyrosol_model_raven.m   build ecTyrosol.mat (RAVEN addMets/addRxns)
│   ├── ecTyrosol_native.mat            model used for ecFactory targets
│   └── ecTyrosol.mat                   alias (envelope / dFBA)
├── scripts/
│   ├── run_ecTyrosol_native.m      ecFactory strain design
│   ├── tyrosol_envelopes.m         measured-strain production envelopes
│   └── tyrosol_dfba.m              dynamic FBA (G5)
├── data/
│   └── ByTyrOH_strain_table.md     lab strain ↔ model gene mapping
├── results/                        target lists (L1, L2, L3, transporters)
└── docs/
    └── METHODS.md                  assumptions and parameters
```

## Reproducibility (ecFactory targets)

```matlab
cd Tyrosol_ecYeast/model
build_ecTyrosol_model_raven     % ecTyrosol.mat + ecTyrosol_native.mat

addpath('../scripts')
run_ecTyrosol_native            % writes results/*.txt
```

## Workflow

1. **Model** — `build_ecTyrosol_model_raven.m` extends `ecYeastGEM_batch.mat`
   with the Ehrlich pathway (ARO10, ADH7) using RAVEN `addMets` / `addRxns`.

2. **Strain design** — `run_ecTyrosol_native.m` runs ecFactory on minimal
   glucose medium and writes ranked targets to `results/`.

3. **Envelopes / dFBA** — `tyrosol_envelopes.m` and `tyrosol_dfba.m` load
   `model/ecTyrosol.mat` (identical to `ecTyrosol_native.mat`).

## Model assumptions

- Native Ehrlich route: 4-HPP → 4-HPAA (ARO10) → tyrosol (ADH7).
- Enzyme-constrained new steps: kcat = 1000 s⁻¹ (GECKO protein coefficients).
- Wild-type background; minimal medium in ecFactory simulation.

## ecFactory settings

| Parameter | Value |
|---|---|
| Medium | Minimal, D-glucose |
| WT biomass yield | 0.48 g / g glucose |
| Target yield in scan | 0.49 × WT yield |

## Requirements

- Python 3 + `scipy` (not required for model build)
- MATLAB + RAVEN + Gurobi + GECKO 2.0.3 + ecFactory
- `ecYeastGEM_batch.mat` from [CellFactory-ecYeastGEM](https://github.com/SysBioChalmers/CellFactory-ecYeastGEM)

See `docs/METHODS.md` for full assumptions and paths.

## Citation

Cite ecFactory, GECKO, ecYeastGEM, and this work as appropriate for your manuscript.
