# Al-Mg-Si contact on Si: interface energetics per source

dH in eV per atom moved; dE_rxn in eV per Si taken from the substrate; x_Si = exp(-dH_Si/kT) at 723 K (configurational entropy only).

| source | T | dH_Si | dH_Mg | Ef(Mg2Si)/f.u. | dE_rxn | B(Mg-Si) | dH_Si paired | x_Si(723 K) |
|---|---|---|---|---|---|---|---|---|
| QE PBE | 0K | +0.371 | +0.129 | -0.480 | -0.739 | +0.022 | +0.349 | 2.6e-03 |
| GPAW r2SCAN | 0K | +0.549 | +0.202 | -0.591 | -0.995 | +0.049 | +0.499 | 1.5e-04 |
| MACE-MP-0, 0 K | 0K | +0.406 | +0.017 | -0.244 | -0.278 | +0.124 | +0.282 | 1.5e-03 |
| MACE-MP-0 MD | 300K | +0.478 | +0.024 | -0.247 | -0.295 | +0.288 | +0.190 | 4.7e-04 |
| MEAM MD | 300K | +0.084 | -0.205 | +0.131 | +0.540 | -0.115 | +0.200 | 2.6e-01 |
| NEP89 MD | 300K | -0.171 | -0.014 | +0.266 | +0.294 | -0.076 | -0.094 | 1.5e+01 |
| MACE-MP-0 MD | 453K | +0.549 | +0.073 | -0.245 | -0.392 | +0.312 | +0.236 | 1.5e-04 |
| MEAM MD | 453K | +0.153 | -0.219 | +0.131 | +0.570 | +0.119 | +0.034 | 8.5e-02 |
| NEP89 MD | 453K | -0.214 | -0.159 | -0.083 | +0.234 | -0.166 | -0.048 | 3.1e+01 |
| MACE-MP-0 MD | 803K | +0.680 | -0.136 | -0.229 | +0.042 | +0.467 | +0.213 | 1.8e-05 |
| MEAM MD | 803K | +0.439 | -0.345 | +0.139 | +0.829 | +0.351 | +0.088 | 8.7e-04 |
| NEP89 MD | 803K | -0.423 | -0.405 | -0.215 | +0.594 | -0.350 | -0.073 | 8.9e+02 |

Mg2Si band gap per source (eV): qe_pbe 0.22, qe_hse06 0.63, mp_gga_expsynth 0.22, jarvis_mbj 0.62, snumat_hse 0.59, jarvis_optb88vdw 0.38, mp_gga 0.22, snumat_gga 0.22, expt_zhuo 2.27, mp_r2scan 0.59
