# Bridge budget: minimum shared-composition overlap for a stable cross-channel offset

Universe: `shared/output/bandgap_archive_build/bandgap_universe.jsonl` (261053 records seen, 150985 used with gap > 0.01 eV in the requested channels, 0 unparseable formulas skipped).

Per pair, ALL shared (canonical reduced) formulas define the ground-truth cross-channel offset delta_full (2-entity Bradley-Terry log-strength difference; delta_AB > 0 means channel A's gaps run higher) and full-overlap Kendall tau. The overlap is then thinned to n compositions (300 seeded subsamples per n) and the offset and tau are refitted.

Tau significance per subsample: scipy.stats.kendalltau two-sided p < 0.05 (exact permutation null at small n without ties, normal approximation otherwise). The tau CI half-width column is the analytic null normal approximation 1.96*sqrt(2(2n+5)/(9n(n-1))).

Certificate threshold (planning defaults): n* = smallest n with sign stability >= 0.95 AND tau power >= 0.8.

Estimator note: under (near-)complete separation (one channel wins almost every bridge event, e.g. hybrid vs GGA of the same database) the Bradley-Terry MLE diverges and the fitted |delta| is limited by the MAP damping of PlackettLuceModel, so the median-absolute-offset-error column can be non-monotone in n for such pairs. Sign stability and tau power are unaffected.

## jarvis_optb88vdw vs jarvis_mbj

- records used: 15020 (jarvis_optb88vdw), 6439 (jarvis_mbj); shared formulas: **5865**
- delta_full = -3.0821 (sign: B higher), tau_full = 0.7911 (p = 0.00e+00)
- n* = 8

| n | sign stability | median abs offset error | tau power | tau CI half-width (null) | tau median [2.5, 97.5 pct] |
|---|---|---|---|---|---|
| 3 | 0.993 | 2.6250 | 0.000 | 1.251 | 1.000 [-0.333, 1.000] |
| 5 | 1.000 | 3.1345 | 0.407 | 0.800 | 0.800 [0.200, 1.000] |
| 8 | 1.000 | 3.6037 | 0.880 | 0.566 | 0.786 [0.429, 1.000] |
| 12 | 1.000 | 4.0088 | 0.987 | 0.433 | 0.818 [0.545, 0.939] |
| 20 | 1.000 | 0.8893 | 1.000 | 0.318 | 0.800 [0.632, 0.905] |
| 35 | 1.000 | 0.4346 | 1.000 | 0.232 | 0.799 [0.682, 0.889] |
| 60 | 1.000 | 0.2804 | 1.000 | 0.174 | 0.795 [0.707, 0.854] |
| 100 | 1.000 | 0.3908 | 1.000 | 0.133 | 0.796 [0.730, 0.844] |
| 200 | 1.000 | 0.2333 | 1.000 | 0.093 | 0.792 [0.754, 0.825] |

## snumat_gga vs snumat_hse

- records used: 8412 (snumat_gga), 8729 (snumat_hse); shared formulas: **8412**
- delta_full = -5.8561 (sign: B higher), tau_full = 0.7908 (p = 0.00e+00)
- n* = 8

| n | sign stability | median abs offset error | tau power | tau CI half-width (null) | tau median [2.5, 97.5 pct] |
|---|---|---|---|---|---|
| 3 | 1.000 | 0.1490 | 0.000 | 1.251 | 1.000 [-0.333, 1.000] |
| 5 | 1.000 | 0.3605 | 0.380 | 0.800 | 0.800 [0.400, 1.000] |
| 8 | 1.000 | 0.8298 | 0.887 | 0.566 | 0.786 [0.429, 1.000] |
| 12 | 1.000 | 1.2348 | 1.000 | 0.433 | 0.818 [0.545, 0.939] |
| 20 | 1.000 | 1.7453 | 1.000 | 0.318 | 0.800 [0.632, 0.911] |
| 35 | 1.000 | 2.3047 | 1.000 | 0.232 | 0.798 [0.682, 0.877] |
| 60 | 1.000 | 2.8436 | 1.000 | 0.174 | 0.788 [0.710, 0.859] |
| 100 | 1.000 | 3.3544 | 1.000 | 0.133 | 0.793 [0.727, 0.845] |
| 200 | 1.000 | 4.0474 | 1.000 | 0.093 | 0.791 [0.751, 0.829] |

## mp_gga vs expt_zhuo

- records used: 67462 (mp_gga), 2468 (expt_zhuo); shared formulas: **1309**
- delta_full = -1.6493 (sign: B higher), tau_full = 0.6066 (p = 3.03e-236)
- n* = 12

| n | sign stability | median abs offset error | tau power | tau CI half-width (null) | tau median [2.5, 97.5 pct] |
|---|---|---|---|---|---|
| 3 | 0.933 | 4.0578 | 0.000 | 1.251 | 1.000 [-1.000, 1.000] |
| 5 | 0.970 | 1.2455 | 0.157 | 0.800 | 0.600 [-0.105, 1.000] |
| 8 | 0.963 | 0.5541 | 0.573 | 0.566 | 0.643 [0.143, 0.966] |
| 12 | 0.993 | 0.6462 | 0.887 | 0.433 | 0.616 [0.266, 0.879] |
| 20 | 1.000 | 0.5434 | 0.993 | 0.318 | 0.628 [0.386, 0.806] |
| 35 | 1.000 | 0.2641 | 1.000 | 0.232 | 0.612 [0.419, 0.746] |
| 60 | 1.000 | 0.2214 | 1.000 | 0.174 | 0.613 [0.485, 0.710] |
| 100 | 1.000 | 0.1653 | 1.000 | 0.133 | 0.609 [0.515, 0.693] |
| 200 | 1.000 | 0.1247 | 1.000 | 0.093 | 0.605 [0.541, 0.665] |

## Bridge-budget table (n* per pair)

| pair | overlap | delta_full | tau_full | n* | n=3 status |
|---|---|---|---|---|---|
| jarvis_optb88vdw vs jarvis_mbj | 5865 | -3.082 | 0.791 | 8 | sign stab 0.99, tau power 0.00 |
| snumat_gga vs snumat_hse | 8412 | -5.856 | 0.791 | 8 | sign stab 1.00, tau power 0.00 |
| mp_gga vs expt_zhuo | 1309 | -1.649 | 0.607 | 12 | sign stab 0.93, tau power 0.00 |

## What n=3 means for the Al-Mg-Si archive gate

The Al-Mg-Si retrospective gate (compute/almgsi_archive/results/TRACKA_GATE.md) found exactly 3 shared ternary compositions between aflow_dHf and mp_dHf, with ternary-only Kendall tau = 1.0 at n=3. The thinning curves above show what n=3 buys in channel pairs where the truth is KNOWN from a rich overlap: the exact permutation test at n=3 can never reach p < 0.05 (its smallest two-sided p-value is 1/3), so the tau power at n=3 is identically zero -- a perfect tau over 3 compositions certifies nothing -- and the offset sign estimated from 3 bridge events flips across subsamples unless the channel separation is very large. The n* column is the bridge budget a prospective Al-Mg-Si campaign must purchase (shared compositions measured/computed in BOTH channels) before a transfer certificate is possible.

Caveat, stated honestly: these are band-gap channel pairs, a proxy for the Al-Mg-Si formation-energy/empirical channel pairs. The transferable content is the SHAPE of the destabilisation curve and the ORDER OF MAGNITUDE of n*, not the exact number.

