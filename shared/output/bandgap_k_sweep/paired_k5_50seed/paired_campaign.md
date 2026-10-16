# Experiment 26: paired cold/full/selective campaign protocol

Receivers: expt_zhuo, jarvis_mbj, jarvis_optb88vdw, mp_gga, mp_gga_expsynth, mp_r2scan, snumat_gga, snumat_hse; k=5; budget=240; 50 seeds (matched across conditions); lambda_full=8.0; early window = first 5 of 10 checkpoints.

Reference-policy mapping: cold=random, full=per-channel, selective=transfer-aware, combined=pooled.

## Condition x arm summary (mean over all donor->receiver cells)

| condition | arm | mean T_early | mean T_final | cells>0 (early) |
|---|---|---|---|---|
| combined | kway | +0.1500 | +0.0377 | 88% |
| combined | pair | +0.2338 | +0.0968 | 100% |
| full | kway | +0.1350 | +0.0195 | 80% |
| full | pair | +0.2311 | +0.0783 | 88% |
| selective | kway | +0.1095 | +0.0176 | 77% |
| selective | pair | +0.1961 | +0.0635 | 95% |

## Hyperedge residual vs clique-expansion cyclic energy (cold k-way stream)

| receiver | mean hyperedge residual | clique cyclic fraction |
|---|---|---|
| expt_zhuo | 0.2233 | 0.0459 |
| jarvis_mbj | 0.3074 | 0.1488 |
| jarvis_optb88vdw | 0.3509 | 0.1970 |
| mp_gga | 0.3515 | 0.2130 |
| mp_gga_expsynth | 0.3011 | 0.1427 |
| mp_r2scan | 0.3227 | 0.1687 |
| snumat_gga | 0.3263 | 0.1573 |
| snumat_hse | 0.3070 | 0.1413 |

Full directed matrices with bootstrap CIs: paired_campaign.json.
Each evaluation = one real archived record (no invented observations); donor priors are ordinal-only.
