# Rock-paper-scissors cycles: global versus local

Archive sha256[:16] = `f2158ecdc78e0229`; bootstrap 200 resamples; campaign replays 50 seeds.

## Whole data sources (global)

| source | classes | 3-cycles | cyclic energy | max triangle curl | P(any cycle, bootstrap) | fragile triangles |
|---|---|---|---|---|---|---|
| MP GGA | 8 | 0 | 0.031 | 0.63 | 0.000 | 0 / 56 |
| MP r2SCAN | 8 | 0 | 0.043 | 0.95 | 0.020 | 1 / 56 |
| MP GGA (ICSD-matched) | 8 | 0 | 0.028 | 0.76 | 0.000 | 0 / 56 |
| Experiment | 7 | 0 | 0.006 | 0.68 | 0.000 | 0 / 35 |
| JARVIS OptB88vdW | 8 | 0 | 0.034 | 0.65 | 0.020 | 3 / 56 |
| JARVIS TB-mBJ | 8 | 0 | 0.042 | 1.15 | 0.015 | 1 / 56 |
| SNUMAT GGA | 8 | 0 | 0.027 | 1.01 | 0.010 | 2 / 56 |
| SNUMAT HSE06 | 8 | 0 | 0.041 | 1.43 | 0.020 | 4 / 56 |

## Sub-populations (local, L1 and L2)

| stratum | source | classes | records | 3-cycles | P(any cycle, bootstrap) | robust cycles | worst triangle (freq) |
|---|---|---|---|---|---|---|---|
| triclinic | MP GGA | 6 | 18005 | 0  | 0.000 | none | none |
| triclinic | MP r2SCAN | 4 | 996 | 0  | 0.000 | none | none |
| triclinic | MP GGA (ICSD-matched) | 6 | 2475 | 0  | 0.000 | none | none |
| triclinic | JARVIS OptB88vdW | 5 | 2340 | 0  | 0.000 | none | none |
| triclinic | JARVIS TB-mBJ | 3 | 424 | 0  | 0.000 | none | none |
| triclinic | SNUMAT GGA | 5 | 1052 | 0  | 0.000 | none | none |
| triclinic | SNUMAT HSE06 | 5 | 1074 | 0  | 0.000 | none | none |
| monoclinic | MP GGA | 8 | 25514 | 0  | 0.010 | none | carbide->nitride->other (0.01) |
| monoclinic | MP r2SCAN | 6 | 3849 | 0  | 0.000 | none | none |
| monoclinic | MP GGA (ICSD-matched) | 8 | 9723 | 0  | 0.000 | none | none |
| monoclinic | JARVIS OptB88vdW | 7 | 5465 | 0  | 0.010 | none | halide->nitride->other (0.01) |
| monoclinic | JARVIS TB-mBJ | 5 | 1586 | 0  | 0.000 | none | none |
| monoclinic | SNUMAT GGA | 5 | 2807 | 0  | 0.000 | none | none |
| monoclinic | SNUMAT HSE06 | 6 | 2904 | 0  | 0.000 | none | none |
| orthorhombic | MP GGA | 8 | 16270 | 0  | 0.000 | none | none |
| orthorhombic | MP r2SCAN | 7 | 3730 | 0  | 0.070 | none | boride->other->pnictide (0.04) |
| orthorhombic | MP GGA (ICSD-matched) | 8 | 7308 | 1 (boride->halide->oxide) | 0.075 | none | boride->halide->oxide (0.06) |
| orthorhombic | JARVIS OptB88vdW | 8 | 4572 | 0  | 0.050 | none | boride->carbide->nitride (0.02) |
| orthorhombic | JARVIS TB-mBJ | 7 | 1614 | 0  | 0.025 | none | carbide->chalcogenide->other (0.01) |
| orthorhombic | SNUMAT GGA | 6 | 2692 | 0  | 0.000 | none | none |
| orthorhombic | SNUMAT HSE06 | 6 | 2772 | 0  | 0.000 | none | none |
| tetragonal | MP GGA | 8 | 6040 | 0  | 0.000 | none | none |
| tetragonal | MP r2SCAN | 6 | 1567 | 0  | 0.000 | none | none |
| tetragonal | MP GGA (ICSD-matched) | 8 | 2368 | 0  | 0.005 | none | carbide->halide->oxide (0.01) |
| tetragonal | JARVIS OptB88vdW | 8 | 2980 | 0  | 0.020 | none | halide->nitride->oxide (0.01) |
| tetragonal | JARVIS TB-mBJ | 7 | 1612 | 0  | 0.005 | none | carbide->chalcogenide->other (0.01) |
| tetragonal | SNUMAT GGA | 6 | 1196 | 0  | 0.000 | none | none |
| tetragonal | SNUMAT HSE06 | 6 | 1269 | 0  | 0.000 | none | none |
| trigonal | MP GGA | 7 | 6474 | 0  | 0.000 | none | none |
| trigonal | MP r2SCAN | 7 | 1848 | 0  | 0.015 | none | carbide->nitride->other (0.01) |
| trigonal | MP GGA (ICSD-matched) | 7 | 2361 | 0  | 0.005 | none | carbide->chalcogenide->other (0.01) |
| trigonal | JARVIS OptB88vdW | 7 | 3173 | 0  | 0.005 | none | carbide->halide->oxide (0.01) |
| trigonal | JARVIS TB-mBJ | 6 | 1420 | 0  | 0.000 | none | none |
| trigonal | SNUMAT GGA | 7 | 1048 | 0  | 0.005 | none | chalcogenide->other->pnictide (0.01) |
| trigonal | SNUMAT HSE06 | 7 | 1080 | 0  | 0.005 | none | chalcogenide->other->pnictide (0.01) |
| hexagonal | MP GGA | 7 | 2209 | 0  | 0.050 | none | carbide->chalcogenide->nitride (0.05) |
| hexagonal | MP r2SCAN | 7 | 720 | 0  | 0.050 | none | carbide->chalcogenide->nitride (0.04) |
| hexagonal | MP GGA (ICSD-matched) | 6 | 1075 | 0  | 0.000 | none | none |
| hexagonal | JARVIS OptB88vdW | 7 | 1295 | 0  | 0.000 | none | none |
| hexagonal | JARVIS TB-mBJ | 7 | 655 | 0  | 0.025 | none | carbide->nitride->other (0.01) |
| hexagonal | SNUMAT GGA | 5 | 349 | 0  | 0.000 | none | none |
| hexagonal | SNUMAT HSE06 | 5 | 373 | 0  | 0.000 | none | none |
| cubic | MP GGA | 8 | 5897 | 0  | 0.020 | none | boride->nitride->oxide (0.01) |
| cubic | MP r2SCAN | 6 | 1364 | 0  | 0.000 | none | none |
| cubic | MP GGA (ICSD-matched) | 6 | 1712 | 0  | 0.000 | none | none |
| cubic | JARVIS OptB88vdW | 8 | 2903 | 0  | 0.035 | none | boride->carbide->nitride (0.01) |
| cubic | JARVIS TB-mBJ | 7 | 1583 | 0  | 0.005 | none | carbide->chalcogenide->nitride (0.01) |
| cubic | SNUMAT GGA | 6 | 793 | 0  | 0.000 | none | none |
| cubic | SNUMAT HSE06 | 6 | 874 | 0  | 0.010 | none | chalcogenide->nitride->other (0.01) |
| binary | MP GGA | 8 | 4768 | 0  | 0.005 | none | boride->halide->oxide (0.01) |
| binary | MP r2SCAN | 7 | 2296 | 0  | 0.000 | none | none |
| binary | MP GGA (ICSD-matched) | 8 | 2259 | 0  | 0.060 | none | boride->halide->oxide (0.04) |
| binary | Experiment | 7 | 1491 | 0  | 0.000 | none | none |
| binary | JARVIS OptB88vdW | 8 | 3304 | 0  | 0.010 | none | boride->chalcogenide->other (0.01) |
| binary | JARVIS TB-mBJ | 8 | 1977 | 0  | 0.005 | none | boride->other->pnictide (0.01) |
| binary | SNUMAT GGA | 7 | 1766 | 0  | 0.000 | none | none |
| binary | SNUMAT HSE06 | 7 | 1899 | 0  | 0.000 | none | none |
| ternary | MP GGA | 8 | 31313 | 0  | 0.000 | none | none |
| ternary | MP r2SCAN | 8 | 10144 | 0  | 0.015 | none | carbide->chalcogenide->other (0.01) |
| ternary | MP GGA (ICSD-matched) | 8 | 10415 | 0  | 0.050 | none | boride->halide->oxide (0.03) |
| ternary | Experiment | 4 | 1321 | 0  | 0.000 | none | none |
| ternary | JARVIS OptB88vdW | 8 | 11512 | 0  | 0.015 | none | carbide->chalcogenide->other (0.01) |
| ternary | JARVIS TB-mBJ | 8 | 5296 | 0  | 0.010 | none | carbide->other->pnictide (0.01) |
| ternary | SNUMAT GGA | 8 | 4860 | 0  | 0.025 | none | boride->other->oxide (0.01) |
| ternary | SNUMAT HSE06 | 8 | 5069 | 0  | 0.065 | none | carbide->chalcogenide->other (0.04) |
| quaternary+ | MP GGA | 8 | 58413 | 0  | 0.005 | none | carbide->chalcogenide->other (0.01) |
| quaternary+ | MP r2SCAN | 4 | 1545 | 0  | 0.010 | none | halide->nitride->oxide (0.01) |
| quaternary+ | MP GGA (ICSD-matched) | 8 | 14341 | 0  | 0.010 | none | carbide->halide->nitride (0.01) |
| quaternary+ | Experiment | 4 | 1038 | 0  | 0.000 | none | none |
| quaternary+ | JARVIS OptB88vdW | 6 | 7782 | 0  | 0.000 | none | none |
| quaternary+ | JARVIS TB-mBJ | 5 | 1615 | 0  | 0.000 | none | none |
| quaternary+ | SNUMAT GGA | 5 | 3353 | 0  | 0.000 | none | none |
| quaternary+ | SNUMAT HSE06 | 6 | 3422 | 0  | 0.000 | none | none |

Summary: 73 stratum tournaments, 2400 triangles; 1 tournaments contain a point-estimate 3-cycle; 0 triangles are cyclic in at least half of the bootstrap resamples.

## Every step of a random campaign (L3)

Hodge energy fractions of the empirical flow on the observed comparison graph; cycles are directed 3-cycles among fully observed triangles.

| k | budget | observed pairs | gradient | local curl | global harmonic | triangles cyclic | P(any cycle) | cycles touching a near-tie pair |
|---|---|---|---|---|---|---|---|---|
| 2 | 8 | 3.7 | 0.983 | 0.017 | 0.000 | 0.013 | 0.013 | 0.20 |
| 2 | 16 | 6.9 | 0.882 | 0.065 | 0.053 | 0.070 | 0.087 | 0.34 |
| 2 | 24 | 9.7 | 0.758 | 0.158 | 0.085 | 0.113 | 0.240 | 0.37 |
| 2 | 32 | 12.2 | 0.663 | 0.244 | 0.093 | 0.122 | 0.388 | 0.40 |
| 2 | 48 | 15.9 | 0.581 | 0.381 | 0.038 | 0.097 | 0.560 | 0.46 |
| 2 | 64 | 18.6 | 0.555 | 0.435 | 0.010 | 0.084 | 0.647 | 0.45 |
| 2 | 96 | 22.4 | 0.510 | 0.490 | 0.000 | 0.083 | 0.770 | 0.44 |
| 2 | 128 | 24.7 | 0.505 | 0.495 | 0.000 | 0.080 | 0.792 | 0.45 |
| 2 | 176 | 26.1 | 0.526 | 0.474 | 0.000 | 0.074 | 0.877 | 0.46 |
| 2 | 240 | 26.8 | 0.565 | 0.435 | 0.000 | 0.068 | 0.870 | 0.50 |
| 2 | 480 | 27.1 | 0.681 | 0.319 | 0.000 | 0.055 | 0.805 | 0.56 |
| 2 | 960 | 27.1 | 0.779 | 0.221 | 0.000 | 0.041 | 0.725 | 0.65 |
| 2 | 2400 | 27.1 | 0.878 | 0.122 | 0.000 | 0.019 | 0.495 | 0.82 |
| 4 | 8 | 10.7 | 0.779 | 0.221 | 0.000 | 0.000 | 0.000 | 0.00 |
| 4 | 16 | 17.1 | 0.682 | 0.318 | 0.000 | 0.024 | 0.278 | 0.44 |
| 4 | 24 | 21.1 | 0.633 | 0.367 | 0.000 | 0.042 | 0.540 | 0.45 |
| 4 | 32 | 23.5 | 0.608 | 0.392 | 0.000 | 0.043 | 0.660 | 0.47 |
| 4 | 48 | 25.9 | 0.606 | 0.394 | 0.000 | 0.047 | 0.710 | 0.50 |
| 4 | 64 | 26.6 | 0.627 | 0.373 | 0.000 | 0.045 | 0.713 | 0.49 |
| 4 | 96 | 27.1 | 0.662 | 0.338 | 0.000 | 0.048 | 0.748 | 0.51 |
| 4 | 128 | 27.1 | 0.697 | 0.303 | 0.000 | 0.045 | 0.740 | 0.53 |
| 4 | 176 | 27.1 | 0.740 | 0.260 | 0.000 | 0.043 | 0.740 | 0.57 |
| 4 | 240 | 27.1 | 0.782 | 0.218 | 0.000 | 0.036 | 0.645 | 0.63 |
| 4 | 480 | 27.1 | 0.854 | 0.146 | 0.000 | 0.022 | 0.520 | 0.72 |
| 4 | 960 | 27.1 | 0.908 | 0.092 | 0.000 | 0.013 | 0.393 | 0.85 |
| 4 | 2400 | 27.1 | 0.944 | 0.056 | 0.000 | 0.006 | 0.223 | 0.98 |

## Data sources as players (L4)

| scope | sources | pairs with >= 30 shared compounds | closed triangles | 3-cycles (bootstrap persistence; cyclic on the common compounds?) |
|---|---|---|---|---|
| all | 8 | 28 | 56 | 2 (jarvis_optb88vdw->mp_gga->mp_gga_expsynth (0.73; no, n = 6475), jarvis_optb88vdw->mp_gga_expsynth->snumat_gga (1.00; no, n = 3683)) |
| boride | 7 | 10 | 5 | 0  |
| carbide | 7 | 15 | 16 | 0  |
| chalcogenide | 8 | 28 | 56 | 0  |
| halide | 8 | 28 | 56 | 1 (expt_zhuo->jarvis_mbj->snumat_hse (0.55; no, n = 89)) |
| nitride | 7 | 21 | 35 | 0  |
| other | 7 | 21 | 35 | 0  |
| oxide | 8 | 28 | 56 | 2 (jarvis_optb88vdw->mp_gga->snumat_gga (0.49; no, n = 1148), jarvis_optb88vdw->mp_gga_expsynth->snumat_gga (0.60; no, n = 1094)) |
| pnictide | 8 | 28 | 56 | 0  |

Edges of the source-level cycles (P = share of shared compounds where the first source reports the wider gap):

- all: jarvis_optb88vdw>mp_gga P = 0.497 over 10811 compounds
- all: mp_gga>mp_gga_expsynth P = 0.452 over 23595 compounds
- all: mp_gga_expsynth>jarvis_optb88vdw P = 0.462 over 6475 compounds
- all: jarvis_optb88vdw>mp_gga_expsynth P = 0.538 over 6475 compounds
- all: mp_gga_expsynth>snumat_gga P = 0.559 over 7344 compounds
- all: snumat_gga>jarvis_optb88vdw P = 0.528 over 3937 compounds
- halide: expt_zhuo>jarvis_mbj P = 0.473 over 93 compounds
- halide: jarvis_mbj>snumat_hse P = 0.387 over 636 compounds
- halide: snumat_hse>expt_zhuo P = 0.473 over 129 compounds
- oxide: jarvis_optb88vdw>mp_gga P = 0.532 over 4184 compounds
- oxide: mp_gga>snumat_gga P = 0.507 over 2661 compounds
- oxide: snumat_gga>jarvis_optb88vdw P = 0.504 over 1192 compounds
- oxide: jarvis_optb88vdw>mp_gga_expsynth P = 0.588 over 2392 compounds
- oxide: mp_gga_expsynth>snumat_gga P = 0.518 over 2510 compounds
