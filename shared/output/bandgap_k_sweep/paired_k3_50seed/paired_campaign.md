# Experiment 26: paired cold/full/selective campaign protocol

Receivers: expt_zhuo, jarvis_mbj, jarvis_optb88vdw, mp_gga, mp_gga_expsynth, mp_r2scan, snumat_gga, snumat_hse; k=3; budget=240; 50 seeds (matched across conditions); lambda_full=8.0; early window = first 5 of 10 checkpoints.

Reference-policy mapping: cold=random, full=per-channel, selective=transfer-aware, combined=pooled.

## Condition x arm summary (mean over all donor->receiver cells)

| condition | arm | mean T_early | mean T_final | cells>0 (early) |
|---|---|---|---|---|
| combined | kway | +0.1749 | +0.0708 | 100% |
| combined | pair | +0.2338 | +0.0968 | 100% |
| full | kway | +0.1843 | +0.0555 | 88% |
| full | pair | +0.2311 | +0.0783 | 88% |
| selective | kway | +0.1369 | +0.0461 | 84% |
| selective | pair | +0.1961 | +0.0635 | 95% |

## Hyperedge residual vs clique-expansion cyclic energy (cold k-way stream)

| receiver | mean hyperedge residual | clique cyclic fraction |
|---|---|---|
| expt_zhuo | 0.2158 | 0.0952 |
| jarvis_mbj | 0.3077 | 0.2986 |
| jarvis_optb88vdw | 0.3485 | 0.3712 |
| mp_gga | 0.3465 | 0.3798 |
| mp_gga_expsynth | 0.3068 | 0.2821 |
| mp_r2scan | 0.3186 | 0.3173 |
| snumat_gga | 0.3269 | 0.3167 |
| snumat_hse | 0.3043 | 0.2781 |

Full directed matrices with bootstrap CIs: paired_campaign.json.
Each evaluation = one real archived record (no invented observations); donor priors are ordinal-only.
