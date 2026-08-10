# ByTyrOH cumulative strain table (in silico mapping)

Laboratory strains built as a **cumulative** engineering series. Each row
adds one or more edits on top of all previous rows. The in silico mapping is
implemented in `scripts/defineByTyrOHStrainTable.m` and used by
`tyrosol_envelopes.m` and `tyrosol_dfba.m`.

Envelopes and dFBA load **`model/ecTyrosol.mat`** (native Ehrlich pathway,
ARO10/ADH7 at kcat = 1000 s⁻¹ on minimal medium).

| Strain | ARO2 OE | ARO7 G141S OE | ARO4 K229L OE | ARO10 OE | LPP1 KD | ALD5 KO | MAE1 KO | ZWF1 KD | ARO1 OE | PHA2 KD | ADH6–ADH7 OE | PDH1 OE |
|--------|---------|---------------|---------------|----------|---------|---------|---------|---------|---------|---------|---------------|---------|
| By4743 wt | − | − | − | − | − | − | − | − | − | − | − | − |
| ByTyrOH 1 | + | − | − | − | − | − | − | − | − | − | − | − |
| ByTyrOH 2 | + | + | − | − | − | − | − | − | − | − | − | − |
| ByTyrOH 3 | + | + | + | + | − | − | − | − | − | − | − | − |
| ByTyrOH 5 | + | + | + | + | + | − | − | − | − | − | − | − |
| ByTyrOH 6 | + | + | + | + | + | + | − | − | − | − | − | − |
| ByTyrOH 7 | + | + | + | + | + | + | + | − | − | − | − | − |
| ByTyrOH 8 | + | + | + | + | + | + | + | + | − | − | − | − |
| ByTyrOH 9 | + | + | + | + | + | + | + | + | + | − | − | − |
| ByTyrOH 10 | + | + | + | + | + | + | + | + | + | + | − | − |
| ByTyrOH 11 | + | + | + | + | + | + | + | + | + | + | + | − |
| ByTyrOH 12 | + | + | + | + | + | + | + | + | + | + | − | + |
| ByTyrOH def | + | + | + | + | + | + | + | + | + | + | + | + |

Lab flask names use **ByTOH1–ByTOH11** (no ByTOH4); cumulative in silico names
skip **ByTyrOH 4** accordingly (sequence 1, 2, 3, 5, …, 11, 12, def).

## Gene IDs (S. cerevisiae)

| Label | Gene ID | Action in model |
|-------|---------|-----------------|
| ARO2 | YGL148W | OE (×1000) |
| ARO7 G141S | YPR060C | OE (×1000) |
| ARO4 K229L | YBR249C | OE (×1000) |
| ARO10 | YDR380W | OE (×1000) |
| LPP1 | YDR503C | KD (×0.21) |
| ALD5 | YER073W | KO |
| MAE1 | YKL029C | KO |
| ZWF1 | YNL241C | KD (×0.21) |
| ARO1 | YDR127W | OE (×13.93) |
| PHA2 | YNL316C | KD (×0.21) |
| ADH6 | YMR318C | OE (×1000) |
| ADH7 | YCR105W | OE (×1000) |
| PDH1 | YER178W (PDA1) | OE (×1000) |

Feedback-resistant alleles are represented as overexpression of the wild-type
gene entry in the ecModel (no separate G141S/K229L entries in ecYeastGEM).
