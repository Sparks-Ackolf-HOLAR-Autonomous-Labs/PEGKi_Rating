# PEGKi transfer scores into the DFT receivers

| Receiver | Donor | Delta | sigma | rho | omega | learned-trust early share | fixed-trust cost to 8 of 10 | donor pool cost / receiver pool cost |
|---|---|---|---|---|---|---|---|---|
| GPAW r2SCAN | MACE-MP-0, 0 K | +0.974 | +0.922 | +0.513 | +0.462 | 0.585 | 0.199 | 0.0012 |
| GPAW r2SCAN | NEP89 MD, 300 K | -0.629 | -0.365 | +0.192 | +0.456 | 0.501 | 0.905 | 0.0073 |
| GPAW r2SCAN | MEAM MD, 300 K | +0.690 | +0.748 | +0.258 | +0.315 | 0.194 | 0.684 | 0.0015 |
| GPAW r2SCAN | QE PBE | +0.993 | +0.909 | +0.473 | +0.388 | 0.712 | 0.089 | 7.6992 |
| GPAW r2SCAN | QE HSE06 | +1.322 | +0.916 | +0.507 | +0.101 | 0.706 | 0.059 | 75.0806 |

GPAW r2SCAN: Spearman(Delta, measured early share) = +0.80; top by Delta = QE HSE06; top measured = QE PBE

Over 20 MD-noise draws: top-Delta pick QE HSE06 100%; Spearman mean +0.80 (range +0.80 to +0.80); Delta 5-95%: MACE-MP-0, 0 K [+0.974, +0.974], NEP89 MD, 300 K [-0.671, -0.543], MEAM MD, 300 K [+0.592, +0.703], QE PBE [+0.993, +0.993], QE HSE06 [+1.322, +1.322]

