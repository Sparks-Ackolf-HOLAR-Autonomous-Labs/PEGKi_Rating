# Hypergraph trial on the frozen band-gap archive universe

Records used (gap > 0.01 eV): 192257 of 261053.
One within-channel event ranks material classes by real archived
gaps (deterministic stride draw); one bridge event ranks the
same formula across channels (identifies fidelity offsets).

| channel | events | PL class ranking (best first) | edge residual | pairwise cyclic energy | tau vs mean order |
|---|---|---|---|---|---|
| expt_zhuo | 300 | halide > oxide > nitride > boride > chalcogenide > pnictide > other | 0.225 | 0.033 | 1.000 |
| jarvis_mbj | 300 | oxide > halide > boride > nitride > carbide > chalcogenide > other > pnictide | 0.326 | 0.037 | 0.643 |
| jarvis_optb88vdw | 300 | halide > carbide > oxide > nitride > boride > chalcogenide > other > pnictide | 0.386 | 0.028 | 0.571 |
| mp_gga | 300 | halide > oxide > boride > nitride > carbide > chalcogenide > pnictide > other | 0.369 | 0.034 | 0.714 |
| mp_gga_expsynth | 300 | boride > halide > oxide > nitride > carbide > chalcogenide > other > pnictide | 0.297 | 0.028 | 0.857 |
| mp_r2scan | 300 | oxide > halide > nitride > boride > carbide > chalcogenide > other > pnictide | 0.345 | 0.041 | 0.571 |
| snumat_gga | 300 | nitride > oxide > halide > boride > carbide > chalcogenide > other > pnictide | 0.348 | 0.069 | 0.571 |
| snumat_hse | 300 | oxide > nitride > halide > carbide > boride > chalcogenide > other > pnictide | 0.317 | 0.022 | 0.643 |

## Pooled cross-channel fit

2400 within-channel events + 4000 bridge events.

Combined class ranking: halide > oxide > nitride > boride > carbide > chalcogenide > other > pnictide

Channel offsets (log wide-gap propensity, descending):

- `snumat_hse`: +2.188
- `expt_zhuo`: +1.441
- `jarvis_mbj`: +1.434
- `mp_r2scan`: +0.443
- `mp_gga_expsynth`: -0.943
- `snumat_gga`: -1.481
- `jarvis_optb88vdw`: -1.505
- `mp_gga`: -1.577

