# Graded cycles: decision bodies versus material classes

Archive sha256[:16] = `f2158ecdc78e0229`; bootstrap 200 resamples.

## Whole tournaments

| tournament | bodies | rounds | loop share (90% CI) | local curl | global harmonic | regularity (null) | triangles | point cycles | expected cycles | undecided triangles |
|---|---|---|---|---|---|---|---|---|---|---|
| ML models (JARVIS-Leaderboard) | 14 | 8 | 0.193 [0.193, 0.207] | 0.193 | 0.000 | 0.76 (0.45) | 364 | 0 | 0.5 | 3 |
| Reuse strategies, pairwise queries | 18 | 4000 | 0.080 [0.062, 0.108] | 0.080 | 0.000 | 0.85 (0.36) | 816 | 4 | 3.0 | 8 |
| Reuse strategies, 4-way events | 18 | 4000 | 0.075 [0.054, 0.101] | 0.075 | 0.000 | 0.85 (0.36) | 816 | 4 | 3.7 | 12 |
| Reuse strategies, 6-way events | 18 | 4000 | 0.064 [0.048, 0.091] | 0.064 | 0.000 | 0.85 (0.36) | 816 | 1 | 2.3 | 8 |
| Acquisition rules (exp 27) | 5 | 240 | 0.012 [0.009, 0.016] | 0.012 | 0.000 | 0.96 (0.92) | 10 | 0 | 0.0 | 0 |
| ML models, all leaderboard tasks | 14 | 43 | 0.142 [0.145, 0.186] | 0.142 | 0.000 | 0.65 (0.45) | 364 | 0 | 2.7 | 14 |
| Acquisition rules, long horizon (exp 27b) | 5 | 40000 | 0.005 [0.003, 0.009] | 0.005 | 0.000 | 0.99 (0.92) | 10 | 0 | 0.0 | 0 |
| Acquisition algorithms (colour testbed) | 9 | 416 | 0.024 [0.023, 0.027] | 0.024 | 0.000 | 0.98 (0.63) | 84 | 0 | 0.0 | 0 |
| Material classes, mp_gga | 8 | records | 0.031 [0.026, 0.037] | 0.031 | 0.000 | 0.99 (0.69) | 56 | 0 | 0.0 | n/a |
| Material classes, mp_r2scan | 8 | records | 0.043 [0.029, 0.061] | 0.043 | 0.000 | 0.98 (0.70) | 56 | 0 | 0.0 | n/a |
| Material classes, mp_gga_expsynth | 8 | records | 0.028 [0.021, 0.035] | 0.028 | 0.000 | 0.98 (0.69) | 56 | 0 | 0.0 | n/a |
| Material classes, expt_zhuo | 7 | records | 0.006 [0.004, 0.015] | 0.006 | 0.000 | 0.93 (0.76) | 35 | 0 | 0.0 | n/a |
| Material classes, jarvis_optb88vdw | 8 | records | 0.034 [0.026, 0.045] | 0.034 | 0.000 | 0.99 (0.69) | 56 | 0 | 0.0 | n/a |
| Material classes, jarvis_mbj | 8 | records | 0.042 [0.032, 0.055] | 0.042 | 0.000 | 0.99 (0.69) | 56 | 0 | 0.0 | n/a |
| Material classes, snumat_gga | 8 | records | 0.027 [0.023, 0.036] | 0.027 | 0.000 | 0.99 (0.71) | 56 | 0 | 0.0 | n/a |
| Material classes, snumat_hse | 8 | records | 0.041 [0.033, 0.052] | 0.041 | 0.000 | 0.99 (0.68) | 56 | 0 | 0.0 | n/a |

Loop share = curl + harmonic energy fraction of the pairwise log-odds flow. Regularity = share of the loop energy in the top Blade-Chest plane (1 = one regular cycle); null = the same share for structureless random flows on the same comparison graph. Expected cycles = sum over triangles of the bootstrap probability that the triangle is cyclic. Undecided = triangles cyclic in between 5% and 95% of resamples.

## Strongest near-cycles per tournament

**ML models (JARVIS-Leaderboard)**: kgcnn_cgcnn > kgcnn_megnet > kgcnn_schnet (strength +0.000, P(cycle) 0.12); kgcnn_coGN > kgcnn_coNGN > potnet (strength +0.000, P(cycle) 0.25); kgcnn_megnet > kgcnn_schnet > matminer_xgboost (strength +0.000, P(cycle) 0.12)

**Reuse strategies, pairwise queries**: full:expt_zhuo > full:mp_gga > selective:jarvis_optb88vdw (strength +0.024, P(cycle) 0.88); full:expt_zhuo > full:jarvis_optb88vdw > selective:jarvis_mbj (strength +0.005, P(cycle) 0.21); full:expt_zhuo > full:mp_gga > selective:jarvis_mbj (strength +0.005, P(cycle) 0.58)

**Reuse strategies, 4-way events**: full:expt_zhuo > full:mp_gga > selective:jarvis_optb88vdw (strength +0.028, P(cycle) 0.81); full:mp_gga > selective:expt_zhuo > selective:jarvis_optb88vdw (strength +0.020, P(cycle) 0.75); full:expt_zhuo > full:jarvis_optb88vdw > full:mp_gga (strength +0.016, P(cycle) 0.62)

**Reuse strategies, 6-way events**: selective:expt_zhuo > selective:jarvis_optb88vdw > selective:mp_gga (strength +0.012, P(cycle) 0.54); full:jarvis_mbj > full:mp_gga_expsynth > selective:mp_r2scan (strength -0.005, P(cycle) 0.11); full:expt_zhuo > full:jarvis_optb88vdw > selective:jarvis_mbj (strength -0.007, P(cycle) 0.12)

**Acquisition rules (exp 27)**: kway_random > kway_active > kway_fisher (strength -0.142, P(cycle) 0.00); kway_random > kway_active > kway_learned (strength -0.206, P(cycle) 0.00); kway_random > kway_fisher > kway_learned (strength -0.206, P(cycle) 0.00)

**ML models, all leaderboard tasks**: cfid > cgcnn_model > kgcnn_dimenetPP (strength +0.000, P(cycle) 0.25); cgcnn_model > kgcnn_dimenetPP > matminer_lgbm (strength +0.000, P(cycle) 0.20); kgcnn_megnet > kgcnn_schnet > matminer_rf (strength -0.020, P(cycle) 0.21)

**Acquisition rules, long horizon (exp 27b)**: kway_random > kway_active > kway_fisher (strength -0.023, P(cycle) 0.00); kway_random > kway_fisher > kway_learned (strength -0.061, P(cycle) 0.00); kway_active > kway_fisher > kway_learned (strength -0.061, P(cycle) 0.00)

**Acquisition algorithms (colour testbed)**: random_forest > thompson_sampling > ucb1_bandit (strength -0.156, P(cycle) 0.00); evolutionary_strategy > grid_search > thompson_sampling (strength -0.159, P(cycle) 0.00); evolutionary_strategy > thompson_sampling > ucb1_bandit (strength -0.159, P(cycle) 0.00)

## Contexts: reversal and commitment regret

| tournament | context | levels | loop share within one level (mean, max) | distinct winners | reversed pairs | max reversal | mean regret | max regret | regret unit |
|---|---|---|---|---|---|---|---|---|---|
| ML models (JARVIS-Leaderboard) | property family | 3 | 0.251, 0.280 | 2 | 10 / 91 | 0.500 | 0.017 | 0.027 | log MAE ratio (0.10 = 10.5% more error) |
| Reuse strategies, pairwise queries | receiver | 8 | 0.034, 0.048 | 5 | 144 / 153 | 0.498 | 0.172 | 0.449 | Kendall tau |
| Reuse strategies, pairwise queries | stage | 2 | 0.083, 0.096 | 2 | 7 / 153 | 0.035 | 0.011 | 0.020 | Kendall tau |
| Reuse strategies, 4-way events | receiver | 8 | 0.019, 0.038 | 8 | 148 / 153 | 0.483 | 0.143 | 0.389 | Kendall tau |
| Reuse strategies, 4-way events | stage | 2 | 0.076, 0.099 | 1 | 8 / 153 | 0.022 | 0.000 | 0.000 | Kendall tau |
| Reuse strategies, 6-way events | receiver | 8 | 0.024, 0.037 | 7 | 144 / 153 | 0.486 | 0.122 | 0.339 | Kendall tau |
| Reuse strategies, 6-way events | stage | 2 | 0.066, 0.085 | 1 | 15 / 153 | 0.049 | 0.000 | 0.000 | Kendall tau |
| Acquisition rules (exp 27) | held-out source | 8 | 0.056, 0.081 | 2 | 5 / 10 | 0.300 | 0.020 | 0.092 | Kendall tau |
| Acquisition rules (exp 27) | stage | 2 | 0.020, 0.030 | 1 | 1 / 10 | 0.009 | 0.000 | 0.000 | Kendall tau |
| ML models, all leaderboard tasks | property family | 5 | 0.223, 0.256 | 3 | 18 / 91 | 0.500 | 0.108 | 0.405 | log MAE ratio (0.10 = 10.5% more error) |
| Acquisition rules, long horizon (exp 27b) | held-out source | 8 | 0.012, 0.028 | 3 | 6 / 10 | 0.111 | 0.021 | 0.097 | Kendall tau |
| Acquisition rules, long horizon (exp 27b) | stage | 3 | 0.006, 0.010 | 1 | 1 / 10 | 0.006 | 0.000 | 0.000 | Kendall tau |
| Acquisition algorithms (colour testbed) | engine | 4 | 0.040, 0.050 | 3 | 1 / 36 | 0.067 | 0.741 | 2.069 | colour distance |

- ML models (JARVIS-Leaderboard), by property family: winners {'electronic': 'potnet', 'energetic': 'potnet', 'mechanical': 'kgcnn_coNGN'}; strongest reversals: kgcnn_coNGN vs potnet (+0.500); kgcnn_dimenetPP vs matminer_lgbm (+0.500); kgcnn_megnet vs kgcnn_schnet (+0.500)
- Reuse strategies, pairwise queries, by receiver: winners {'expt_zhuo': 'full:snumat_gga', 'jarvis_mbj': 'full:mp_gga_expsynth', 'jarvis_optb88vdw': 'full:mp_gga', 'mp_gga': 'full:jarvis_mbj', 'mp_gga_expsynth': 'full:jarvis_mbj', 'mp_r2scan': 'full:mp_gga_expsynth', 'snumat_gga': 'full:snumat_hse', 'snumat_hse': 'full:snumat_gga'}; strongest reversals: full:mp_gga_expsynth vs full:snumat_hse (+0.498); full:jarvis_optb88vdw vs full:snumat_gga (+0.497); full:mp_gga_expsynth vs full:snumat_gga (+0.497)
- Reuse strategies, pairwise queries, by stage: winners {'early': 'full:mp_r2scan', 'late': 'full:jarvis_mbj'}; strongest reversals: full:mp_gga_expsynth vs selective:mp_gga_expsynth (+0.035); full:jarvis_mbj vs full:mp_gga_expsynth (+0.029); selective:jarvis_optb88vdw vs selective:mp_gga (+0.018)
- Reuse strategies, 4-way events, by receiver: winners {'expt_zhuo': 'full:snumat_gga', 'jarvis_mbj': 'selective:mp_gga_expsynth', 'jarvis_optb88vdw': 'full:mp_gga', 'mp_gga': 'full:jarvis_mbj', 'mp_gga_expsynth': 'full:mp_r2scan', 'mp_r2scan': 'full:mp_gga_expsynth', 'snumat_gga': 'selective:snumat_hse', 'snumat_hse': 'selective:snumat_gga'}; strongest reversals: full:mp_gga_expsynth vs full:snumat_hse (+0.483); full:expt_zhuo vs full:jarvis_mbj (+0.478); full:snumat_gga vs selective:mp_gga_expsynth (+0.476)
- Reuse strategies, 4-way events, by stage: winners {'early': 'full:mp_r2scan', 'late': 'full:mp_r2scan'}; strongest reversals: selective:jarvis_optb88vdw vs selective:mp_gga (+0.022); full:jarvis_mbj vs full:mp_gga_expsynth (+0.020); combined vs full:jarvis_mbj (+0.016)
- Reuse strategies, 6-way events, by receiver: winners {'expt_zhuo': 'full:snumat_gga', 'jarvis_mbj': 'full:mp_gga_expsynth', 'jarvis_optb88vdw': 'full:mp_gga', 'mp_gga': 'full:jarvis_mbj', 'mp_gga_expsynth': 'full:mp_r2scan', 'mp_r2scan': 'full:mp_gga_expsynth', 'snumat_gga': 'full:snumat_hse', 'snumat_hse': 'selective:snumat_gga'}; strongest reversals: full:jarvis_optb88vdw vs full:snumat_hse (+0.486); full:jarvis_mbj vs full:snumat_hse (+0.483); full:mp_gga_expsynth vs full:snumat_hse (+0.483)
- Reuse strategies, 6-way events, by stage: winners {'early': 'full:mp_r2scan', 'late': 'full:mp_r2scan'}; strongest reversals: full:expt_zhuo vs selective:jarvis_mbj (+0.049); full:expt_zhuo vs selective:mp_gga_expsynth (+0.031); full:jarvis_mbj vs full:mp_gga_expsynth (+0.029)
- Acquisition rules (exp 27), by held-out source: winners {'expt_zhuo': 'kway_learned', 'jarvis_mbj': 'kway_random', 'jarvis_optb88vdw': 'kway_learned', 'mp_gga': 'kway_learned', 'mp_gga_expsynth': 'kway_random', 'mp_r2scan': 'kway_learned', 'snumat_gga': 'kway_learned', 'snumat_hse': 'kway_learned'}; strongest reversals: kway_random vs kway_learned (+0.300); kway_random vs kway_active (+0.267); kway_active vs kway_fisher (+0.250)
- Acquisition rules (exp 27), by stage: winners {'early': 'kway_learned', 'late': 'kway_learned'}; strongest reversals: kway_random vs kway_active (+0.009)
- ML models, all leaderboard tasks, by property family: winners {'calculation settings': 'kgcnn_coNGN', 'electronic': 'potnet', 'energetic': 'potnet', 'mechanical': 'kgcnn_coNGN', 'porous (MOF)': 'matminer_xgboost'}; strongest reversals: alignn_model vs kgcnn_coGN (+0.500); cfid vs kgcnn_dimenetPP (+0.500); cfid vs kgcnn_megnet (+0.500)
- Acquisition rules, long horizon (exp 27b), by held-out source: winners {'expt_zhuo': 'kway_learned', 'jarvis_mbj': 'kway_fisher', 'jarvis_optb88vdw': 'kway_learned', 'mp_gga': 'kway_learned', 'mp_gga_expsynth': 'kway_fisher', 'mp_r2scan': 'kway_learned', 'snumat_gga': 'kway_learned', 'snumat_hse': 'kway_active'}; strongest reversals: kway_fisher vs kway_learned (+0.111); kway_active vs kway_learned (+0.101); kway_random vs kway_learned (+0.045)
- Acquisition rules, long horizon (exp 27b), by stage: winners {'1: to 240': 'kway_learned', '2: 240-960': 'kway_learned', '3: 960-1920': 'kway_learned'}; strongest reversals: kway_random vs kway_active (+0.006)
- Acquisition algorithms (colour testbed), by engine: winners {'km': 'simulated_annealing', 'mixbox': 'bayesian_ei', 'ryb': 'bayesian_ei', 'spectral': 'bayesian_ucb'}; strongest reversals: random_forest vs ucb1_bandit (+0.067)

## Material classes across data sources

8 of 21 class pairs reverse between at least two sources; max reversal strength +0.116. Strongest: boride vs oxide (+0.116); boride vs nitride (+0.108); boride vs halide (+0.108)
