# Pair combinatorics behind the set-operation and hypergraph results

Archive sha256[:16] = `f2158ecdc78e0229`. P_ij = Pr(random record of class i has the larger gap).

## Base data sources

| source | classes | pairs | resolvability r | near-tie pairs | majority 3-cycles | median-vs-majority discordant pairs | order ceiling tau | median spread (eV) |
|---|---|---|---|---|---|---|---|---|
| Experiment | 7 | 21 | 0.539 | 2 | 0 | 1 | +0.905 | 2.53 |
| JARVIS TB-mBJ | 8 | 28 | 0.363 | 2 | 0 | 1 | +0.929 | 5.07 |
| JARVIS OptB88vdW | 8 | 28 | 0.266 | 8 | 0 | 3 | +0.786 | 1.52 |
| MP GGA | 8 | 28 | 0.275 | 3 | 0 | 0 | +1.000 | 1.87 |
| MP GGA (ICSD-matched) | 8 | 28 | 0.372 | 4 | 0 | 0 | +1.000 | 3.32 |
| MP r2SCAN | 8 | 28 | 0.327 | 5 | 0 | 2 | +0.857 | 2.25 |
| SNUMAT GGA | 8 | 28 | 0.312 | 5 | 0 | 2 | +0.857 | 1.55 |
| SNUMAT HSE06 | 8 | 28 | 0.354 | 4 | 0 | 1 | +0.929 | 2.12 |

## One full round: predicted versus observed residual

| source | observed edge residual | predicted (majority order) | predicted (median order) | observed cyclic fraction | isotropic-noise cyclic fraction | majority 3-cycles |
|---|---|---|---|---|---|---|
| Experiment | 0.225 | 0.231 | 0.232 | 0.033 | 0.714 | 0 |
| JARVIS TB-mBJ | 0.326 | 0.319 | 0.322 | 0.037 | 0.750 | 0 |
| JARVIS OptB88vdW | 0.386 | 0.367 | 0.377 | 0.028 | 0.750 | 0 |
| MP GGA | 0.369 | 0.363 | 0.363 | 0.034 | 0.750 | 0 |
| MP GGA (ICSD-matched) | 0.297 | 0.314 | 0.314 | 0.028 | 0.750 | 0 |
| MP r2SCAN | 0.345 | 0.336 | 0.339 | 0.041 | 0.750 | 0 |
| SNUMAT GGA | 0.348 | 0.344 | 0.348 | 0.069 | 0.750 | 0 |
| SNUMAT HSE06 | 0.317 | 0.323 | 0.323 | 0.022 | 0.750 | 0 |

Identity check: max |observed - predicted| = 0.019, mean 0.009, Pearson +0.984; sources with a directed 3-cycle: 0.

## k-sweep (50 seeds, budget 240) against the pair structure

| k | mean final tau | mean gap to order ceiling | Spearman(r, tau) over 8 sources | Spearman(-near ties, tau) |
|---|---|---|---|---|
| 2 | +0.644 | +0.264 | +0.81 | +0.29 |
| 3 | +0.689 | +0.219 | +0.88 | +0.40 |
| 4 | +0.721 | +0.187 | +0.88 | +0.40 |
| 6 | +0.722 | +0.186 | +0.81 | +0.29 |

## Covering designs (Schonheim lower bound)

n = 7 classes: k=2: 21 rounds, 42 evaluations; k=3: 7 rounds, 21 evaluations; k=4: 4 rounds, 16 evaluations; k=5: 3 rounds, 15 evaluations; k=6: 3 rounds, 18 evaluations; k=7: 1 rounds, 7 evaluations
n = 8 classes: k=2: 28 rounds, 56 evaluations; k=3: 11 rounds, 33 evaluations; k=4: 6 rounds, 24 evaluations; k=5: 4 rounds, 20 evaluations; k=6: 3 rounds, 18 evaluations; k=7: 3 rounds, 21 evaluations; k=8: 1 rounds, 8 evaluations

## Set-operation regions

| donor to receiver | Delta(agree) - Delta(disagree) | r agree | r disagree | near ties agree / disagree | classes |
|---|---|---|---|---|---|
| MP GGA to Experiment | -0.044 | 0.665 | 0.578 | 0 / 1 | 4 |
| SNUMAT GGA to SNUMAT HSE06 | -0.507 | 0.395 | 0.462 | 2 / 3 | 7 |
| JARVIS OptB88vdW to JARVIS TB-mBJ | -0.018 | 0.340 | 0.387 | 3 / 4 | 7 |

Components into agree / disagree regions:

- MP GGA to Experiment: concordance 0.568 / 0.550; sigma 0.886 / 0.902; rho 0.205 / 0.226; omega 0.261 / 0.255
- SNUMAT GGA to SNUMAT HSE06: concordance 0.218 / 0.436; sigma 0.427 / 0.855; rho 0.327 / 0.349; omega 0.179 / 0.122
- JARVIS OptB88vdW to JARVIS TB-mBJ: concordance 0.316 / 0.353; sigma 0.876 / 0.878; rho 0.402 / 0.323; omega 0.222 / 0.128

| term | Spearman vs pair concordance | Spearman vs receiver resolvability | mean | sd |
|---|---|---|---|---|
| sigma | +0.570 | +0.298 | 0.855 | 0.117 |
| rho | -0.315 | -0.400 | 0.365 | 0.073 |
| omega | +0.600 | +0.721 | 0.165 | 0.098 |
| delta | -0.338 | -0.605 | 1.055 | 0.194 |

Pair-concordance predictor versus committed Delta over 140 routes: Pearson -0.170, Spearman -0.338; receiver resolvability alone, Spearman -0.605.
