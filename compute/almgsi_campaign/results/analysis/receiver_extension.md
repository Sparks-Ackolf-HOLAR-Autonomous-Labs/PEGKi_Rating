# Added Al-Mg-Si receivers: PEGKi transfer scores and measured campaigns

| Receiver | Donor | tau (pool) | Delta | learned-trust early share | fixed-trust cost to 8 of 10 (k=4) |
|---|---|---|---|---|---|
| MACE-MP-0 MD, 300 K | MACE-MP-0, 0 K | +0.88 | +1.186 | 0.746 | 0.040 |
| MACE-MP-0 MD, 300 K | NEP89 MD, 300 K | -0.04 | -0.392 | 0.487 | 0.838 |
| MACE-MP-0 MD, 300 K | MEAM MD, 300 K | +0.42 | +0.664 | 0.449 | 0.572 |
| MACE-MP-0 MD, 300 K | MP archive | +0.73 | +1.122 | 0.602 | 0.419 |
| MACE-MP-0 MD, 300 K | QE PBE | +0.80 | +0.954 | 0.606 | 0.153 |
| MACE-MP-0 MD, 300 K | QE HSE06 | +0.73 | +1.071 | 0.572 | 0.160 |

MACE-MP-0 MD, 300 K (n = 168): no-donor early share 0.454; Spearman(Delta, measured) = +0.77; top by Delta = MACE-MP-0, 0 K; top measured = MACE-MP-0, 0 K

- Shapley (pooled_fixed): MACE-MP-0, 0 K +0.160, QE PBE +0.088, MP archive +0.086, QE HSE06 +0.072, MEAM MD, 300 K -0.028, NEP89 MD, 300 K -0.094; all donors 0.737, best coalition 0.825
- Shapley (pooled_learned): MACE-MP-0, 0 K +0.113, MP archive +0.042, QE HSE06 +0.030, QE PBE +0.028, NEP89 MD, 300 K -0.005, MEAM MD, 300 K -0.007; all donors 0.655, best coalition 0.746

| MP archive | MACE-MP-0, 0 K | +0.77 | +0.763 | 0.202 | 0.480 |
| MP archive | NEP89 MD, 300 K | +0.33 | -0.248 | 0.158 | 0.640 |
| MP archive | MEAM MD, 300 K | -0.16 | -0.340 | 0.030 | 0.960 |
| MP archive | QE PBE | +0.78 | +0.768 | 0.202 | 0.480 |
| MP archive | QE HSE06 | +0.73 | +0.951 | 0.202 | 0.320 |

MP archive (n = 25): no-donor early share 0.109; Spearman(Delta, measured) = +0.89; top by Delta = QE HSE06; top measured = MACE-MP-0, 0 K

- Shapley (pooled_fixed): QE HSE06 +0.055, QE PBE +0.044, MACE-MP-0, 0 K +0.028, NEP89 MD, 300 K +0.023, MEAM MD, 300 K -0.068; all donors 0.192, best coalition 0.204
- Shapley (pooled_learned): QE HSE06 +0.046, QE PBE +0.035, MACE-MP-0, 0 K +0.028, NEP89 MD, 300 K +0.001, MEAM MD, 300 K -0.036; all donors 0.182, best coalition 0.202

| MACE-MP-0 quasi-harmonic, 453 K | MACE-MP-0, 0 K | +0.78 | +1.117 | 0.395 | 0.231 |
| MACE-MP-0 quasi-harmonic, 453 K | QE PBE | +0.62 | +0.706 | 0.319 | 0.234 |
| MACE-MP-0 quasi-harmonic, 453 K | MACE-MP-0 MD, 453 K | +0.77 | +1.313 | 0.371 | 0.205 |
| MACE-MP-0 quasi-harmonic, 453 K | MEAM MD, 453 K | +0.19 | +0.243 | 0.319 | 0.954 |
| MACE-MP-0 quasi-harmonic, 453 K | NEP89 MD, 453 K | +0.35 | +0.769 | 0.255 | 0.625 |

MACE-MP-0 quasi-harmonic, 453 K (n = 49): no-donor early share 0.238; Spearman(Delta, measured) = +0.60; top by Delta = MACE-MP-0 MD, 453 K; top measured = MACE-MP-0, 0 K

- Shapley (pooled_fixed): MACE-MP-0 MD, 453 K +0.079, MACE-MP-0, 0 K +0.077, QE PBE +0.066, NEP89 MD, 453 K -0.023, MEAM MD, 453 K -0.066; all donors 0.371, best coalition 0.427
- Shapley (pooled_learned): MEAM MD, 453 K +0.052, MACE-MP-0, 0 K +0.050, MACE-MP-0 MD, 453 K +0.043, QE PBE +0.031, NEP89 MD, 453 K -0.046; all donors 0.368, best coalition 0.447

| MACE-MP-0 quasi-harmonic, 803 K | MACE-MP-0, 0 K | +0.58 | +1.041 | 0.330 | 0.380 |
| MACE-MP-0 quasi-harmonic, 803 K | QE PBE | +0.50 | +0.618 | 0.311 | 0.283 |
| MACE-MP-0 quasi-harmonic, 803 K | MACE-MP-0 MD, 803 K | +0.56 | +1.152 | 0.338 | 0.448 |
| MACE-MP-0 quasi-harmonic, 803 K | MEAM MD, 803 K | +0.08 | +0.144 | 0.240 | 0.943 |
| MACE-MP-0 quasi-harmonic, 803 K | NEP89 MD, 803 K | +0.35 | +0.702 | 0.233 | 0.715 |

MACE-MP-0 quasi-harmonic, 803 K (n = 49): no-donor early share 0.220; Spearman(Delta, measured) = +0.70; top by Delta = MACE-MP-0 MD, 803 K; top measured = MACE-MP-0 MD, 803 K

- Shapley (pooled_fixed): MACE-MP-0, 0 K +0.065, QE PBE +0.063, MACE-MP-0 MD, 803 K +0.062, NEP89 MD, 803 K +0.047, MEAM MD, 803 K -0.161; all donors 0.296, best coalition 0.375
- Shapley (pooled_learned): MACE-MP-0 MD, 803 K +0.041, QE PBE +0.039, MACE-MP-0, 0 K +0.039, NEP89 MD, 803 K +0.010, MEAM MD, 803 K -0.005; all donors 0.346, best coalition 0.346

