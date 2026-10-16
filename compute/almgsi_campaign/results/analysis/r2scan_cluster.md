# Solute-cluster binding: r2SCAN against PBE

GPAW r2SCAN single points at the QE PBE geometry (108-site Al cells). Binding in eV, positive = attractive.

| Cluster | PBE | r2SCAN |
|---|---|---|
| Mg2Si2_tet | +0.127 | +0.049 |
| Mg2Si_tri | +0.045 | +0.088 |
| Mg3Si_tet | +0.044 | +0.089 |
| MgSi2_tri | +0.054 | +0.060 |
| MgSi3_tet | +0.095 | +0.054 |
| MgSi_1nn | +0.022 | +0.049 |
| MgSi_2nn | +0.004 | +0.002 |
| SiSi_1nn | -0.002 | -0.063 |
| SiSi_2nn | +0.004 | -0.050 |

PBE vs r2SCAN: n = 9, tau +0.33 [-0.24, +0.81], sign agreement 8/9, mean abs. difference 0.040 eV.

| Source | tau vs PBE | tau vs r2SCAN |
|---|---|---|
| MACE-MP-0, 0 K | +0.78 [+0.50, +1.00] | +0.22 [-0.23, +0.55] |
| MACE-MP-0 MD, 300 K | +0.78 [+0.52, +1.00] | +0.44 [-0.09, +0.86] |
| MEAM MD, 300 K | -0.67 [-1.00, -0.35] | -0.56 [-0.88, -0.09] |
| NEP89 MD, 300 K | -0.17 [-0.73, +0.43] | -0.39 [-0.87, +0.16] |
