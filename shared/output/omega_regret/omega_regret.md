# Regret-based omega from the paired-campaign logs

Inputs: 400 (receiver, seed) checkpoint files from paired_campaign_k4 (matched seeds, fixed-trust full reuse); early window = first 5 of 10 checkpoints; interim omega from omega_decomposition.json. Bootstrap CIs over seeds.

omega_regret(A->B) = mean over matched seeds of the average over early-window checkpoints of max(0, tau_cold - tau_full) under fixed-trust full reuse -- the expected early harm of following donor A's advice on receiver B. Pairwise arm primary; k-way reported in omega_regret.json.

## Correlation with the interim dispersion omega (12 matched directed pairs)

Spearman rho = +0.5167; Kendall tau = +0.4187 (n = 12).

| donor | receiver | interim omega_total | omega_regret [95% CI] | omega_regret_final |
|---|---|---|---|---|
| snumat_gga | jarvis_mbj | 0.0475 | 0.0520 [0.0357, 0.0707] | 0.0829 |
| snumat_hse | jarvis_mbj | 0.0433 | 0.0520 [0.0347, 0.0696] | 0.0829 |
| jarvis_mbj | snumat_hse | 0.0445 | 0.0498 [0.0347, 0.0656] | 0.0586 |
| jarvis_mbj | snumat_gga | 0.0444 | 0.0402 [0.0283, 0.0525] | 0.0543 |
| jarvis_optb88vdw | snumat_hse | 0.0819 | 0.0317 [0.0220, 0.0426] | 0.0629 |
| jarvis_optb88vdw | snumat_gga | 0.0791 | 0.0262 [0.0167, 0.0374] | 0.0686 |
| snumat_gga | jarvis_optb88vdw | 0.0735 | 0.0211 [0.0130, 0.0303] | 0.0371 |
| snumat_hse | jarvis_optb88vdw | 0.0719 | 0.0211 [0.0129, 0.0297] | 0.0371 |
| jarvis_optb88vdw | jarvis_mbj | 0.0416 | 0.0139 [0.0069, 0.0214] | 0.0186 |
| jarvis_mbj | jarvis_optb88vdw | 0.0327 | 0.0071 [0.0034, 0.0117] | 0.0143 |
| snumat_gga | snumat_hse | 0.0200 | 0.0000 [0.0000, 0.0000] | 0.0000 |
| snumat_hse | snumat_gga | 0.0156 | 0.0000 [0.0000, 0.0000] | 0.0000 |

## Largest-regret cells (all 56 campaign cells, pair arm)

| rank | donor | receiver | omega_regret [95% CI] | omega_regret_final |
|---|---|---|---|---|
| 1 | mp_gga | expt_zhuo | 0.1506 [0.1251, 0.1758] | 0.1219 |
| 2 | mp_gga_expsynth | expt_zhuo | 0.1324 [0.1081, 0.1576] | 0.1048 |
| 3 | mp_gga | snumat_hse | 0.0986 [0.0761, 0.1211] | 0.1143 |
| 4 | mp_gga | snumat_gga | 0.0835 [0.0623, 0.1062] | 0.1143 |
| 5 | jarvis_mbj | expt_zhuo | 0.0816 [0.0626, 0.1013] | 0.0819 |
| 6 | snumat_gga | mp_gga | 0.0814 [0.0590, 0.1043] | 0.1157 |
| 7 | snumat_hse | mp_gga | 0.0814 [0.0608, 0.1051] | 0.1157 |
| 8 | jarvis_optb88vdw | expt_zhuo | 0.0782 [0.0594, 0.0995] | 0.0610 |

## Full pair-arm matrix (omega_regret, early window)

| donor \ receiver | expt_zhuo | jarvis_mbj | jarvis_optb88vdw | mp_gga | mp_gga_expsynth | mp_r2scan | snumat_gga | snumat_hse |
|---|---|---|---|---|---|---|---|---|
| expt_zhuo | - | 0.0092 | 0.0123 | 0.0275 | 0.0243 | 0.0077 | 0.0000 | 0.0011 |
| jarvis_mbj | 0.0816 | - | 0.0071 | 0.0000 | 0.0000 | 0.0007 | 0.0402 | 0.0498 |
| jarvis_optb88vdw | 0.0782 | 0.0139 | - | 0.0017 | 0.0099 | 0.0017 | 0.0262 | 0.0317 |
| mp_gga | 0.1506 | 0.0006 | 0.0000 | - | 0.0024 | 0.0039 | 0.0835 | 0.0986 |
| mp_gga_expsynth | 0.1324 | 0.0000 | 0.0026 | 0.0006 | - | 0.0000 | 0.0291 | 0.0362 |
| mp_r2scan | 0.0782 | 0.0006 | 0.0003 | 0.0043 | 0.0003 | - | 0.0135 | 0.0156 |
| snumat_gga | 0.0154 | 0.0520 | 0.0211 | 0.0814 | 0.0486 | 0.0201 | - | 0.0000 |
| snumat_hse | 0.0154 | 0.0520 | 0.0211 | 0.0814 | 0.0486 | 0.0201 | 0.0000 | - |

The regret omega exists for every cell the campaigns measured (including the experimental channel expt_zhuo); the interim dispersion omega exists only where phase 3 ran. K-way arm, CIs, and the final-checkpoint variant: omega_regret.json. Computed entirely from committed campaign logs; no new observations.
