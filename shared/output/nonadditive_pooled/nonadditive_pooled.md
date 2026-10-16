# Non-additive pooled model: learned departure from the additive cold start

Archive sha256[:16] = `f2158ecdc78e0229`; sources: expt_zhuo, jarvis_mbj, jarvis_optb88vdw, mp_gga, mp_gga_expsynth, mp_r2scan, snumat_gga, snumat_hse.
Events: 2400 within-source class rankings, 12000 cross-source bridge rankings; split 60/20/20 (train/validation/test).
Penalty grid (pull toward additive): (0.03, 0.3, 3.0, 30.0, 300.0, 3000.0); chosen on validation, refitted on train+validation, scored on test.

| model | chosen lambda | parameters | held-out gain per event vs additive [95% CI] | within-source | bridge |
|---|---|---|---|---|---|
| additive | None | 16 | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] |
| scale | 30.0 | 24 | +0.0122 [+0.0056, +0.0188] | +0.0656 [+0.0366, +0.0941] | +0.0016 [-0.0032, +0.0062] |
| scale+r1 | 30.0 | 40 | +0.0136 [+0.0070, +0.0206] | +0.0600 [+0.0285, +0.0924] | +0.0044 [-0.0007, +0.0092] |
| scale+r2 | 30.0 | 56 | +0.0157 [+0.0087, +0.0225] | +0.0627 [+0.0299, +0.0930] | +0.0063 [+0.0011, +0.0119] |

## Learned margin scales a_c (scale model; 1 = additive)

| source | a_c | bootstrap 95% CI |
|---|---|---|
| expt_zhuo | 1.658 | [1.529, 1.784] |
| snumat_hse | 1.235 | [1.152, 1.308] |
| snumat_gga | 1.140 | [1.074, 1.216] |
| jarvis_mbj | 1.100 | [1.024, 1.188] |
| mp_gga_expsynth | 0.922 | [0.874, 0.995] |
| mp_gga | 0.804 | [0.741, 0.886] |
| mp_r2scan | 0.743 | [0.682, 0.807] |
| jarvis_optb88vdw | 0.707 | [0.645, 0.763] |

## Class-dependent ladder spacings implied by `scale+r2`

| spacing | oxide | halide | chalcogenide | nitride |
|---|---|---|---|---|
| snumat_hse-snumat_gga | +3.771 | +3.718 | +3.584 | +3.671 |
| jarvis_mbj-jarvis_optb88vdw | +3.000 | +3.186 | +2.831 | +2.918 |
| expt_zhuo-mp_gga | +3.152 | +3.095 | +2.390 | +3.045 |

Family-resolved spacings fitted separately (scripts/family_offsets.py):

| spacing | oxide | halide | chalcogenide | nitride |
|---|---|---|---|---|
| snumat_hse-snumat_gga | +3.707 | +3.753 | +3.541 | +3.822 |
| jarvis_mbj-jarvis_optb88vdw | +2.892 | +3.450 | +2.956 | +2.964 |
| expt_zhuo-mp_gga | +2.702 | +2.643 | +2.605 | +4.260 |

Runtime 16.8 s.
