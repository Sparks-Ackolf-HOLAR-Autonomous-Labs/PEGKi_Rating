# Measured formation-enthalpy receiver: all beta, gamma and Mg2Si combinations

Script: `../../measured_enthalpy_combinations.py`. Data: `measured_enthalpy_combinations.json`.

## Measurements (kJ/mol of atoms; eV/atom)

| Compound | Measurement | kJ/mol-atoms | eV/atom | Source |
|---|---|---|---|---|
| beta | BP1970 | -2.47 | -0.0256 | Brown and Pratt 1970, doi:10.1007/bf03037810, via Mezbahul-Islam et al. 2014 Table 4 |
| beta | PH1978 | -1.22 | -0.0126 | Predel and Huelse 1978, Z. Metallkd. 69, 661, via Mezbahul-Islam et al. 2014 Table 4 |
| gamma | BP1970 | -3.80 | -0.0394 | Brown and Pratt 1970, doi:10.1007/bf03037810, via Mezbahul-Islam et al. 2014 Table 4 |
| gamma | PH1978 | -1.01 | -0.0105 | Predel and Huelse 1978, Z. Metallkd. 69, 661, via Mezbahul-Islam et al. 2014 Table 4 |
| Mg2Si | JANAF | -25.94 | -0.2689 | NIST-JANAF table Mg-034, doi:10.18434/T42S31 |
| Mg2Si | Feufel1997 | -22.20 | -0.2301 | Feufel et al. 1997, doi:10.1016/s0925-8388(96)02655-2, as quoted by Jund et al. 2013, doi:10.1088/0953-8984/25/3/035403 |

Brown and Pratt place gamma below beta; Predel and Huelse place beta below gamma. Mg2Si is the most strongly bound compound in every combination.

## Computed values (eV/atom)

| Source | gamma | Mg2Si | beta, hull at x_Mg = 0.389 | beta, epsilon proxy |
|---|---|---|---|---|
| QE PBE, 0 K | -0.0235 | -0.1600 | -0.0182 | -0.0203 |
| QE HSE06, 0 K | -0.0121 | -0.1297 | -0.0237 | -0.0150 |
| MACE-MP-0, 0 K | -0.0322 | -0.0813 | -0.0498 | -0.0461 |
| MACE-MP-0 MD, 300 K | -0.0321 | -0.0825 | -0.0483 | -0.0462 |
| MEAM MD, 300 K | +0.0544 | +0.0435 | -0.1130 | +0.1082 |
| NEP89 MD, 300 K | -0.0227 | +0.0888 | -0.0283 | -0.0317 |
| MEAM, 0 K | +0.0476 | +0.0423 | -0.1140 | +0.1284 |
| NEP89, 0 K | -0.0223 | +0.0881 | -0.0294 | -0.0328 |
| MP archive | -0.0093 | -0.1086 | -0.0113 | -0.0102 |

## Summary over the 8 combinations, beta from each source's Al-Mg hull (primary)

| Source | MAE range, eV/atom | Kendall tau values | beta-gamma order agrees | binding scale range |
|---|---|---|---|---|
| QE PBE, 0 K | 0.030 to 0.044 | +0.33, +1.00 | 4 of 8 | 0.60 to 0.70 |
| QE HSE06, 0 K | 0.035 to 0.059 | +0.33, +1.00 | 4 of 8 | 0.48 to 0.57 |
| MACE-MP-0, 0 K | 0.060 to 0.082 | +0.33, +1.00 | 4 of 8 | 0.31 to 0.39 |
| MACE-MP-0 MD, 300 K | 0.059 to 0.081 | +0.33, +1.00 | 4 of 8 | 0.32 to 0.39 |
| MEAM MD, 300 K | 0.142 to 0.169 | -0.33, +0.33 | 4 of 8 | -0.20 to -0.13 |
| NEP89 MD, 300 K | 0.111 to 0.130 | -1.00, -0.33 | 4 of 8 | -0.37 to -0.30 |
| MEAM, 0 K | 0.140 to 0.166 | -0.33, +0.33 | 4 of 8 | -0.19 to -0.12 |
| NEP89, 0 K | 0.111 to 0.130 | -1.00, -0.33 | 4 of 8 | -0.37 to -0.30 |
| MP archive | 0.041 to 0.068 | +0.33, +1.00 | 4 of 8 | 0.40 to 0.47 |

## Summary over the 8 combinations, beta mapped to the epsilon phase (alternative)

| Source | MAE range, eV/atom | Kendall tau values | beta-gamma order agrees | binding scale range |
|---|---|---|---|---|
| QE PBE, 0 K | 0.030 to 0.044 | +0.33, +1.00 | 4 of 8 | 0.60 to 0.70 |
| QE HSE06, 0 K | 0.035 to 0.059 | +0.33, +1.00 | 4 of 8 | 0.48 to 0.57 |
| MACE-MP-0, 0 K | 0.059 to 0.081 | +0.33, +1.00 | 4 of 8 | 0.31 to 0.38 |
| MACE-MP-0 MD, 300 K | 0.059 to 0.081 | +0.33, +1.00 | 4 of 8 | 0.32 to 0.39 |
| MEAM MD, 300 K | 0.153 to 0.180 | +0.33, +1.00 | 4 of 8 | -0.27 to -0.19 |
| NEP89 MD, 300 K | 0.112 to 0.131 | -1.00, -0.33 | 4 of 8 | -0.37 to -0.30 |
| MEAM, 0 K | 0.157 to 0.184 | +0.33, +1.00 | 4 of 8 | -0.27 to -0.19 |
| NEP89, 0 K | 0.112 to 0.131 | -1.00, -0.33 | 4 of 8 | -0.37 to -0.29 |
| MP archive | 0.042 to 0.069 | +0.33, +1.00 | 4 of 8 | 0.40 to 0.47 |

Binding scale: least-squares slope of computed against measured through the origin (1 = the source reproduces the measured binding; below 1 = underbinding).

## Per combination (primary rule): Kendall tau / MAE eV/atom

| Combination | QE PBE, 0 K | QE HSE06, 0 K | MACE-MP-0, 0 K | MACE-MP-0 MD, 300 K | MEAM MD, 300 K | NEP89 MD, 300 K | MEAM, 0 K | NEP89, 0 K | MP archive |
|---|---|---|---|---|---|---|---|---|---|
| beta:BP1970|gamma:BP1970|Mg2Si:JANAF | +1.00 / 0.044 | +0.33 / 0.056 | +0.33 / 0.073 | +0.33 / 0.072 | -0.33 / 0.165 | -1.00 / 0.126 | -0.33 / 0.162 | -1.00 / 0.126 | +0.33 / 0.068 |
| beta:BP1970|gamma:BP1970|Mg2Si:Feufel1997 | +1.00 / 0.031 | +0.33 / 0.043 | +0.33 / 0.060 | +0.33 / 0.059 | -0.33 / 0.152 | -1.00 / 0.113 | -0.33 / 0.149 | -1.00 / 0.113 | +0.33 / 0.055 |
| beta:BP1970|gamma:PH1978|Mg2Si:JANAF | +0.33 / 0.043 | +1.00 / 0.048 | +1.00 / 0.078 | +1.00 / 0.077 | +0.33 / 0.155 | -0.33 / 0.124 | +0.33 / 0.153 | -0.33 / 0.124 | +1.00 / 0.059 |
| beta:BP1970|gamma:PH1978|Mg2Si:Feufel1997 | +0.33 / 0.030 | +1.00 / 0.035 | +1.00 / 0.065 | +1.00 / 0.064 | +0.33 / 0.142 | -0.33 / 0.111 | +0.33 / 0.140 | -0.33 / 0.111 | +1.00 / 0.046 |
| beta:PH1978|gamma:BP1970|Mg2Si:JANAF | +1.00 / 0.043 | +0.33 / 0.059 | +0.33 / 0.077 | +0.33 / 0.076 | -0.33 / 0.169 | -1.00 / 0.130 | -0.33 / 0.166 | -1.00 / 0.130 | +0.33 / 0.064 |
| beta:PH1978|gamma:BP1970|Mg2Si:Feufel1997 | +1.00 / 0.031 | +0.33 / 0.046 | +0.33 / 0.064 | +0.33 / 0.064 | -0.33 / 0.156 | -1.00 / 0.117 | -0.33 / 0.154 | -1.00 / 0.117 | +0.33 / 0.051 |
| beta:PH1978|gamma:PH1978|Mg2Si:JANAF | +0.33 / 0.043 | +1.00 / 0.051 | +1.00 / 0.082 | +1.00 / 0.081 | +0.33 / 0.159 | -0.33 / 0.129 | +0.33 / 0.157 | -0.33 / 0.129 | +1.00 / 0.054 |
| beta:PH1978|gamma:PH1978|Mg2Si:Feufel1997 | +0.33 / 0.030 | +1.00 / 0.038 | +1.00 / 0.069 | +1.00 / 0.068 | +0.33 / 0.146 | -0.33 / 0.116 | +0.33 / 0.144 | -0.33 / 0.116 | +1.00 / 0.041 |
