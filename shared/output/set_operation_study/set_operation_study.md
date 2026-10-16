# Set-operation studies on the band-gap archive

Archive sha256[:16] = `f2158ecdc78e0229`; scissor tolerance 0.3 eV; 50 record-bootstrap replicates (95% intervals).

## H1: does a donor transfer better where its scissor-corrected gaps agree?

| donor to receiver | shared compounds | scissor shift (eV) | agree / disagree records | Delta into agree | Delta into disagree | difference [95% CI] |
|---|---|---|---|---|---|---|
| MP GGA to Experiment | 1288 | +0.65 | 741 / 1582 | +0.829 | +0.874 | -0.044 [-0.124, +0.019] |
| SNUMAT GGA to SNUMAT HSE06 | 8412 | +1.01 | 4767 / 5330 | +0.575 | +1.081 | -0.507 [-0.659, -0.267] |
| JARVIS OptB88vdW to JARVIS TB-mBJ | 5865 | +0.86 | 2527 / 5804 | +1.056 | +1.074 | -0.018 [-0.134, +0.033] |

## H2: best base donor for each derived region

| region | best donor | Delta(best) minus Delta(defining donor) [95% CI] |
|---|---|---|
| Experiment where MP GGA agrees with | SNUMAT HSE06 | +0.250 [+0.172, +0.341] |
| Experiment where MP GGA disagrees with | SNUMAT HSE06 | +0.265 [+0.166, +0.351] |
| SNUMAT HSE06 where SNUMAT GGA agrees with | JARVIS OptB88vdW | +0.201 [+0.032, +0.285] |
| SNUMAT HSE06 where SNUMAT GGA disagrees with | JARVIS TB-mBJ | +0.027 [-0.056, +0.112] |
| JARVIS TB-mBJ where JARVIS OptB88vdW agrees with | MP r2SCAN | +0.168 [+0.027, +0.206] |
| JARVIS TB-mBJ where JARVIS OptB88vdW disagrees with | SNUMAT HSE06 | +0.258 [+0.136, +0.333] |

## H3: does a donor's ordering carry to compounds it never reported?

| donor to receiver | compounds covered / uncovered | Delta into covered | Delta into uncovered | difference [95% CI] |
|---|---|---|---|---|
| MP GGA to Experiment | 1288 / 1183 | +0.864 | +0.838 | +0.026 [-0.054, +0.116] |
| MP GGA to JARVIS TB-mBJ | 4864 / 1575 | +1.184 | +1.001 | +0.183 [+0.046, +0.371] |
| JARVIS OptB88vdW to SNUMAT HSE06 | 4004 / 4725 | +1.012 | +0.901 | +0.111 [-0.001, +0.216] |

## Maximin routing over all studies

Weight on derived (set-operation) studies: 0.381

- JARVIS OptB88vdW: 0.348
- JARVIS TB-mBJ where JARVIS OptB88vdW agrees with: 0.273
- MP r2SCAN: 0.271
- Experiment where MP GGA disagrees with: 0.108
