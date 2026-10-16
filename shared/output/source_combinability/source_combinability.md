# Multi-source set operations: which data sources can be combined

Universe sha256[:16] = f2158ecdc78e0229; tolerance 0.3 eV; 200 bootstrap resamples; pairs with fewer than 30 shared compounds are not assessable.

| A | B | shared | coverage | scissor (eV) | values within tol | Kendall tau | margin scale B/A |
|---|---|---|---|---|---|---|---|
| Experiment | SNUMAT HSE06 | 790 | 0.076 | -0.00 | 0.39 [0.36, 0.42] | +0.63 [+0.61, +0.66] | 1.00 |
| Experiment | JARVIS TB-mBJ | 494 | 0.059 | -0.17 | 0.36 [0.33, 0.39] | +0.60 [+0.57, +0.63] | 0.96 |
| Experiment | MP r2SCAN | 735 | 0.056 | -0.32 | 0.38 [0.35, 0.41] | +0.62 [+0.60, +0.65] | 0.84 |
| Experiment | SNUMAT GGA | 745 | 0.073 | -0.78 | 0.35 [0.32, 0.38] | +0.60 [+0.57, +0.63] | 0.81 |
| Experiment | JARVIS OptB88vdW | 692 | 0.041 | -0.79 | 0.30 [0.28, 0.33] | +0.56 [+0.54, +0.59] | 0.76 |
| Experiment | MP GGA | 1288 | 0.019 | -0.66 | 0.34 [0.32, 0.37] | +0.61 [+0.59, +0.63] | 0.85 |
| Experiment | MP GGA (ICSD-matched) | 1102 | 0.044 | -0.70 | 0.35 [0.32, 0.37] | +0.61 [+0.59, +0.64] | 0.85 |
| SNUMAT HSE06 | JARVIS TB-mBJ | 2205 | 0.170 | -0.11 | 0.42 [0.40, 0.44] | +0.67 [+0.66, +0.69] | 1.03 |
| SNUMAT HSE06 | MP r2SCAN | 4277 | 0.269 | -0.54 | 0.54 [0.53, 0.56] | +0.70 [+0.68, +0.71] | 0.85 |
| SNUMAT HSE06 | SNUMAT GGA | 8412 | 0.964 | -1.01 | 0.49 [0.48, 0.50] | +0.79 [+0.79, +0.80] | 0.79 |
| SNUMAT HSE06 | JARVIS OptB88vdW | 4004 | 0.203 | -1.04 | 0.41 [0.39, 0.42] | +0.59 [+0.58, +0.61] | 0.83 |
| SNUMAT HSE06 | MP GGA | 7727 | 0.113 | -1.01 | 0.43 [0.42, 0.44] | +0.65 [+0.64, +0.66] | 0.84 |
| SNUMAT HSE06 | MP GGA (ICSD-matched) | 7411 | 0.297 | -1.00 | 0.44 [0.43, 0.45] | +0.66 [+0.65, +0.67] | 0.84 |
| JARVIS TB-mBJ | MP r2SCAN | 3088 | 0.209 | -0.40 | 0.38 [0.37, 0.40] | +0.79 [+0.78, +0.80] | 0.86 |
| JARVIS TB-mBJ | SNUMAT GGA | 2139 | 0.168 | -0.78 | 0.34 [0.33, 0.36] | +0.66 [+0.64, +0.67] | 0.78 |
| JARVIS TB-mBJ | JARVIS OptB88vdW | 5865 | 0.376 | -0.86 | 0.31 [0.30, 0.32] | +0.79 [+0.78, +0.80] | 0.78 |
| JARVIS TB-mBJ | MP GGA | 4864 | 0.070 | -0.84 | 0.29 [0.28, 0.31] | +0.74 [+0.73, +0.75] | 0.75 |
| JARVIS TB-mBJ | MP GGA (ICSD-matched) | 3270 | 0.122 | -0.88 | 0.30 [0.28, 0.31] | +0.75 [+0.74, +0.76] | 0.76 |
| MP r2SCAN | SNUMAT GGA | 4224 | 0.270 | -0.39 | 0.70 [0.68, 0.71] | +0.74 [+0.73, +0.75] | 0.92 |
| MP r2SCAN | JARVIS OptB88vdW | 5403 | 0.256 | -0.47 | 0.72 [0.71, 0.73] | +0.85 [+0.84, +0.86] | 0.92 |
| MP r2SCAN | MP GGA | 11409 | 0.169 | -0.43 | 0.64 [0.63, 0.64] | +0.83 [+0.82, +0.83] | 0.89 |
| MP r2SCAN | MP GGA (ICSD-matched) | 8163 | 0.304 | -0.43 | 0.67 [0.66, 0.68] | +0.84 [+0.83, +0.84] | 0.91 |
| SNUMAT GGA | JARVIS OptB88vdW | 3937 | 0.202 | -0.02 | 0.73 [0.72, 0.74] | +0.67 [+0.65, +0.68] | 1.06 |
| SNUMAT GGA | MP GGA | 7641 | 0.112 | +0.01 | 0.77 [0.76, 0.78] | +0.75 [+0.74, +0.76] | 1.08 |
| SNUMAT GGA | MP GGA (ICSD-matched) | 7344 | 0.298 | +0.01 | 0.79 [0.78, 0.80] | +0.76 [+0.75, +0.77] | 1.08 |
| JARVIS OptB88vdW | MP GGA | 10811 | 0.151 | +0.00 | 0.73 [0.72, 0.74] | +0.77 [+0.77, +0.78] | 0.93 |
| JARVIS OptB88vdW | MP GGA (ICSD-matched) | 6475 | 0.201 | -0.01 | 0.76 [0.75, 0.77] | +0.82 [+0.81, +0.82] | 0.95 |
| MP GGA | MP GGA (ICSD-matched) | 23595 | 0.350 | +0.00 | 0.97 [0.97, 0.97] | +0.98 [+0.98, +0.98] | 1.01 |

Groups at 1 - tau <= 0.2 (average linkage):

- SNUMAT HSE06: common compounds 8729, union 8729, min pairwise tau n/a
- SNUMAT GGA: common compounds 8412, union 8412, min pairwise tau n/a
- MP r2SCAN, JARVIS OptB88vdW, MP GGA, MP GGA (ICSD-matched): common compounds 4157, union 71700, min pairwise tau +0.77
- JARVIS TB-mBJ: common compounds 6439, union 6439, min pairwise tau n/a
- Experiment: common compounds 2471, union 2471, min pairwise tau n/a

Cut sweep (each row: the tau range of the cut that gives this grouping; singletons omitted):

| cut tau range | groups with two or more sources |
|---|---|
| +0.978 to +1.000 | none |
| +0.849 to +0.978 | {MP GGA, MP GGA (ICSD-matched)} |
| +0.814 to +0.849 | {MP GGA, MP GGA (ICSD-matched)}; {MP r2SCAN, JARVIS OptB88vdW} |
| +0.791 to +0.814 | {MP r2SCAN, JARVIS OptB88vdW, MP GGA, MP GGA (ICSD-matched)} |
| +0.766 to +0.791 | {SNUMAT HSE06, SNUMAT GGA}; {MP r2SCAN, JARVIS OptB88vdW, MP GGA, MP GGA (ICSD-matched)} |
| +0.684 to +0.766 | {SNUMAT HSE06, SNUMAT GGA}; {JARVIS TB-mBJ, MP r2SCAN, JARVIS OptB88vdW, MP GGA, MP GGA (ICSD-matched)} |
| +0.606 to +0.684 | {SNUMAT HSE06, JARVIS TB-mBJ, MP r2SCAN, SNUMAT GGA, JARVIS OptB88vdW, MP GGA, MP GGA (ICSD-matched)} |
| +0.000 to +0.606 | {Experiment, SNUMAT HSE06, JARVIS TB-mBJ, MP r2SCAN, SNUMAT GGA, JARVIS OptB88vdW, MP GGA, MP GGA (ICSD-matched)} |

Compounds reported by all eight sources: 330
