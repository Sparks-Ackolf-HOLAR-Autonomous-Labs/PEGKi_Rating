# Solute-cluster binding energies against PBE (Al-Mg-Si)

PBE clusters: 9; attractive at PBE: 8

| Source | n | Kendall tau vs PBE (90% CI) | sign agreement | MAD (eV) | attractive | tau vs PBE on compounds |
|---|---|---|---|---|---|---|
| MACE-MP-0, 0 K | 9 | +0.78 [+0.48, +1.00] | 0.89 | 0.153 | 9/9 | +0.87 |
| MACE-MP-0 MD, 300 K | 9 | +0.78 [+0.52, +1.00] | 0.89 | 0.279 | 9/9 | +0.61 |
| MEAM MD, 300 K | 9 | -0.67 [-1.00, -0.35] | 0.11 | 0.501 | 0/9 | +0.49 |
| NEP89 MD, 300 K | 9 | -0.17 [-0.75, +0.42] | 0.22 | 0.067 | 3/9 | -0.18 |

## Binding energies (eV, positive = attractive)

| Cluster | PBE | MACE-MP-0, 0 K | MACE-MP-0 MD, 300 K | MEAM MD, 300 K | NEP89 MD, 300 K |
|---|---|---|---|---|---|
| Mg2Si2_tet | +0.127 | +0.389 | +0.540 | -0.703 | +0.018 |
| Mg2Si_tri | +0.045 | +0.151 | +0.268 | -0.345 | -0.017 |
| Mg3Si_tet | +0.044 | +0.106 | +0.365 | -0.908 | -0.032 |
| MgMg_1nn | n/a | -0.037 | +0.172 | -0.200 | -0.094 |
| MgSi2_tri | +0.054 | +0.319 | +0.500 | -0.490 | -0.019 |
| MgSi3_tet | +0.095 | +0.568 | +0.713 | -1.118 | -0.040 |
| MgSi_1nn | +0.022 | +0.124 | +0.171 | -0.136 | -0.003 |
| MgSi_2nn | +0.004 | +0.019 | +0.143 | -0.215 | -0.066 |
| SiSi_1nn | -0.002 | +0.071 | +0.121 | -0.075 | +0.025 |
| SiSi_2nn | +0.004 | +0.022 | +0.081 | -0.125 | +0.037 |
