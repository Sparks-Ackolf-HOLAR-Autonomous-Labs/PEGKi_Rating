# Routing rules across query arity k

Archive sha256[:16] = `f2158ecdc78e0229`; seeds per run {3: 50, 4: 50, 5: 50, 6: 50}; pairwise arm identical across runs: True. Values: mean measured early transfer gain on held-out receivers (leave one receiver out).

| k | maximin linear program | uniform over donors | learned from LP, default score | learned from LP, learned score | best measured donor | learned minus LP |
|---|---|---|---|---|---|---|
| 2 | +0.231 | +0.231 | +0.340 | +0.348 | +0.424 | +0.109 |
| 3 | +0.185 | +0.184 | +0.280 | +0.286 | +0.355 | +0.095 |
| 4 | +0.147 | +0.147 | +0.234 | +0.243 | +0.307 | +0.087 |
| 5 | +0.134 | +0.135 | +0.220 | +0.222 | +0.288 | +0.085 |
| 6 | +0.130 | +0.130 | +0.212 | +0.217 | +0.281 | +0.082 |

Best measured donor per receiver:

- k=2: expt_zhuo: snumat_gga, jarvis_mbj: mp_gga_expsynth, jarvis_optb88vdw: mp_gga, mp_gga: jarvis_mbj, mp_gga_expsynth: jarvis_mbj, mp_r2scan: mp_gga_expsynth, snumat_gga: snumat_hse, snumat_hse: snumat_gga
- k=3: expt_zhuo: snumat_gga, jarvis_mbj: mp_gga_expsynth, jarvis_optb88vdw: mp_gga, mp_gga: jarvis_mbj, mp_gga_expsynth: jarvis_mbj, mp_r2scan: mp_gga_expsynth, snumat_gga: snumat_hse, snumat_hse: snumat_gga
- k=4: expt_zhuo: snumat_gga, jarvis_mbj: mp_gga_expsynth, jarvis_optb88vdw: mp_gga, mp_gga: jarvis_mbj, mp_gga_expsynth: mp_r2scan, mp_r2scan: mp_gga_expsynth, snumat_gga: snumat_hse, snumat_hse: snumat_gga
- k=5: expt_zhuo: snumat_gga, jarvis_mbj: mp_gga_expsynth, jarvis_optb88vdw: mp_gga, mp_gga: jarvis_mbj, mp_gga_expsynth: mp_r2scan, mp_r2scan: mp_gga_expsynth, snumat_gga: snumat_hse, snumat_hse: snumat_gga
- k=6: expt_zhuo: snumat_gga, jarvis_mbj: mp_gga_expsynth, jarvis_optb88vdw: mp_gga, mp_gga: jarvis_mbj, mp_gga_expsynth: mp_r2scan, mp_r2scan: mp_gga_expsynth, snumat_gga: snumat_hse, snumat_hse: snumat_gga
