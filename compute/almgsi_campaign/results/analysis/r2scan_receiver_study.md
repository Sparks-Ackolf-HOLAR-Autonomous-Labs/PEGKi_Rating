# GPAW r2SCAN as a third DFT receiver

180 pool compounds with a GPAW r2SCAN formation energy (PAW, 700 eV, libxc r2SCAN).

| Donor | n | Kendall tau | MAD (eV/atom) | top-10 overlap |
|---|---|---|---|---|
| QE PBE | 168 | +0.93 | 0.029 | 9 |
| QE HSE06 | 168 | +0.89 | 0.012 | 9 |
| MACE-MP-0, 0 K | 168 | +0.86 | 0.039 | 7 |
| MACE-MP-0 MD, 300 K | 16 | +0.58 | 0.017 | 10 |
| NEP89 MD, 300 K | 168 | -0.12 | 0.131 | 4 |
| MEAM MD, 300 K | 165 | +0.47 | 0.059 | 1 |

MP r2SCAN (VASP) against GPAW r2SCAN: tau +0.90, MAD 0.012 eV/atom (n = 15).
GPAW PBE against QE PBE (code check): tau +0.89, MAD 0.013 eV/atom (n = 168).

Mg2Si (eV/atom): GPAW r2SCAN -0.197, MP r2SCAN -0.194, calorimetric -0.269.

Ten most stable r2SCAN compounds: mp_1367_Mg2Si, mp_1101851_Mg2Si, mp_1074652_Mg2Si, mp_1094116_MgAl2, mp_3833_Mg(AlSi)2, mp_4499_MgAlSi, fcc_MgAl7_19, mp_1073794_MgSi, fcc_MgAl7_112, fcc_MgAl3_02; not in the MP archive: 3.
