# MP r2SCAN (archived) against the computed sources

26 pool compounds carry an MP r2SCAN formation energy; comparisons use compounds (elements dropped) that both sources and MP r2SCAN hold.

| Donor | Receiver | n | Kendall tau | MAD (eV/atom) |
|---|---|---|---|---|
| MP GGA/GGA+U | MP r2SCAN | 15 | +0.68 | 0.020 |
| QE PBE | MP r2SCAN | 15 | +0.83 | 0.018 |
| QE HSE06 | MP r2SCAN | 15 | +0.45 | 0.027 |
| MACE-MP-0 0 K | MP r2SCAN | 15 | +0.68 | 0.036 |
| MEAM 0 K | MP r2SCAN | 15 | +0.01 | 0.086 |
| NEP89 0 K | MP r2SCAN | 15 | +0.20 | 0.051 |
| MP r2SCAN | QE PBE | 15 | +0.83 | 0.018 |
| MP GGA/GGA+U | QE PBE | 15 | +0.85 | 0.009 |
| QE HSE06 | QE PBE | 15 | +0.47 | 0.026 |
| MACE-MP-0 0 K | QE PBE | 15 | +0.66 | 0.019 |
| MEAM 0 K | QE PBE | 15 | +0.03 | 0.077 |
| NEP89 0 K | QE PBE | 15 | +0.30 | 0.036 |
| MP r2SCAN | QE HSE06 | 15 | +0.45 | 0.027 |
| MP GGA/GGA+U | QE HSE06 | 15 | +0.47 | 0.020 |
| QE PBE | QE HSE06 | 15 | +0.47 | 0.026 |
| MACE-MP-0 0 K | QE HSE06 | 15 | +0.62 | 0.034 |
| MEAM 0 K | QE HSE06 | 15 | -0.16 | 0.094 |
| NEP89 0 K | QE HSE06 | 15 | -0.01 | 0.056 |

Mg2Si (mp-1367) formation energy, eV/atom; calorimetric -0.269:

- MP r2SCAN: -0.194
- MP GGA/GGA+U: -0.129
- QE PBE: -0.160
- QE HSE06: -0.130
- MACE-MP-0 0 K: -0.081
- MEAM 0 K: +0.042
- NEP89 0 K: +0.088
