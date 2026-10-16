# Literature candidates for measured receivers (Al-Mg-Si)

Search date: 4 Oct 2026. Tools: Asta snippet search and paper lookup,
OpenAlex, Crossref,
Europe PMC, Mendeley Data API. Every DOI below was resolved at
`https://doi.org/<doi>` with `Accept: application/vnd.citationstyles.csl+json`;
titles, authors and years in this file come from those records. BibTeX for
all 27 DOIs (fetched from doi.org with `Accept: application/x-bibtex`) is in
`refs_candidates.bib`.

Conversion used: 1 eV/atom = 96.485 kJ/mol of atoms.

Value status key: **verified** = read by us in retrieved text (page or table
named); **not extracted** = the source is relevant but the number was not
read (paywall, bot block, or figure only).

## A. Measured formation enthalpies

### A1. Values extracted

| Phase | Value, kJ/mol-atoms | Value, eV/atom | Method | Temperature | Source DOI | Where in source | Access |
|---|---|---|---|---|---|---|---|
| beta, Al140Mg89 (Al3Mg2) | -2.47 | -0.0256 | experimental (Brown and Pratt 1970, ref. 138 of the review); method not extracted | not extracted | review 10.1155/2014/704283; primary 10.1007/bf03037810 | review Table 4, p. 5 (read from rendered page image) | review CC BY (publisher site bot-blocked; read from Wayback copy of 13 Feb 2020); primary paywalled |
| beta, Al140Mg89 (Al3Mg2) | -1.22 | -0.0126 | experimental (Predel and Huelse 1978, Z. Metallkd. 69, 661, ref. 139); method not extracted | not extracted | review 10.1155/2014/704283; primary has no DOI found in Crossref | review Table 4, p. 5 | review CC BY; primary print only |
| gamma, Al12Mg17 | -3.80 | -0.0394 | experimental (Brown and Pratt 1970) | not extracted | review 10.1155/2014/704283; primary 10.1007/bf03037810 | review Table 4, p. 5 | as above |
| gamma, Al12Mg17 | -1.01 | -0.0105 | experimental (Predel and Huelse 1978) | not extracted | review 10.1155/2014/704283 | review Table 4, p. 5 | as above |
| epsilon, Al30Mg23 | no experimental value in the review | | only first-principles (-3.42) and CALPHAD (-0.99) values listed | | 10.1155/2014/704283 | review Table 4, p. 5 | |
| Mg2Si | -21.2 to -29.9 (range) | -0.22 to -0.31 (range) | range of literature experiments, as summarised by Jund et al. | not stated | 10.1088/0953-8984/25/3/035403 | Section on formation energies, text citing refs. 36 (Feufel et al. 1997) and 37 (Yuan et al. 2009); arXiv:1309.7246 full text | free (arXiv) |
| Mg2Si | about -22.2 | about -0.23 | "most accurate experiments" (solution calorimetry; mass-spectrometric thermal dissociation gives -0.22 to -0.25 eV/atom) | not stated | 10.1088/0953-8984/25/3/035403, citing 10.1016/s0925-8388(96)02655-2 | same paragraph | free (arXiv) |
| Mg2Si (already used) | -25.94 | -0.2689 | NIST-JANAF assessment | 298.15 K | 10.18434/T42S31 | table Mg-034 (see `LAB_RECEIVERS.md`) | open |

Remarks on A1:

- Reading caution: the plain-text extraction of the review's Table 4
  misaligns the phase labels; the grouping above was checked against the
  rendered page image. Brown and Pratt values are about two to four times
  more negative than Predel and Huelse for the same phase, so a receiver
  built from them should carry this spread as an uncertainty band (for
  gamma, -1.0 to -3.8 kJ/mol-atoms).
- With the existing PBE result for Mg2Si (-0.160 eV/atom), the second Mg2Si
  source (about -0.23 eV/atom) shrinks the PBE error from +0.109 to about
  +0.07 eV/atom; the measured range itself (0.09 eV/atom wide) is comparable
  to the spread between DFT sources.
- With beta and gamma added, an Al-Mg-Si enthalpy receiver would have three
  compounds (Mg2Si, beta, gamma); epsilon has no measured value in the
  sources read.

### A2. Sources that tabulate or report further measured enthalpies (values not extracted)

| Source (DOI) | Content relevant to A | Access |
|---|---|---|
| Zhong, Yang, Liu, Calphad 29 (2005) 303, 10.1016/j.calphad.2005.08.004 | Al-Mg assessment combining first-principles and experiment; compares formation enthalpies of beta, gamma, epsilon with measurements | paywalled |
| Saunders, Calphad 14 (1990) 61, 10.1016/0364-5916(90)90040-7 | review and assessment of Al-Mg, including experimental thermochemistry | paywalled |
| Liang et al., Thermochim. Acta 314 (1998) 87, 10.1016/s0040-6031(97)00458-9 | Al-Mg-Zn experiments and assessment (Al-Mg binary data reviewed) | paywalled |
| Brown and Pratt, Metall. Trans. 1 (1970) 2743, 10.1007/bf03037810 | primary measurement for beta and gamma (values in A1 via review) | paywalled |
| Feufel, Goedecke, Lukas, Sommer, J. Alloys Compd. 247 (1997) 31, 10.1016/s0925-8388(96)02655-2 | Al-Mg-Si experiments and calculation; calorimetric Mg2Si enthalpy (cited by Jund et al. as the most accurate) | paywalled |
| Yuan, Sun, Du, Zhao, Yang, Calphad 33 (2009) 673, 10.1016/j.calphad.2009.08.004 | Mg-Si assessment with a compilation of Mg2Si data | paywalled |
| Povoden-Karadeniz et al., Calphad 43 (2013) 94, 10.1016/j.calphad.2013.03.004 | CALPHAD Gibbs energies of metastable beta'', beta', B', U1, U2; would give assessed (not measured) enthalpies for the metastable phases | paywalled |
| Starink, Int. Mater. Rev. 49 (2004) 191, 10.1179/095066004225010532 | review of quantitative DSC of Al alloys, including heats of precipitation reactions | paywalled (Southampton eprint listed but bot-blocked) |
| Giersberg, Milkereit, Schick, Kessler, Mater. Sci. Forum 794-796 (2014) 939, 10.4028/www.scientific.net/msf.794-796.939 | in situ isothermal calorimetry of precipitation heat in Al-Mg-Si | open (hybrid) |
| Miyauchi, Fujikawa, Hirano, J. Jpn. Inst. Light Met. 21 (1971) 565, 10.2464/jilm.21.565 | specific-heat (heat evolution) curves of four quasi-binary Al-Mg2Si alloys (0.60 to 1.31 wt% Mg2Si) | open (J-STAGE, Japanese) |

No source found reports a standard enthalpy of formation for beta'',
beta', U1, U2 or B' from measurement. DSC papers report heats of
precipitation per gram of alloy, which convert to an enthalpy per mole of
precipitate only with an assumed precipitate composition and volume
fraction; such a receiver would be model-dependent.

## B. Age-hardening datasets with many compositions

| Source (DOI) | Number of alloys | Composition range, wt% | Ageing temperature(s) | Property | Data access |
|---|---|---|---|---|---|
| Bhat, Barnard, Birbilis, R. Soc. Open Sci. 10 (2023), 10.1098/rsos.220360; dataset 10.17632/tvtg7gs59p.1 | 1154 alloys in total; 217 rows (68 distinct Mg-Si-Cu compositions) in a 6xxx-like window (Mg 0.2 to 1.6, Si 0.2 to 1.8, Cu at most 1.0, Zn at most 0.5) by our filter of the file | Mg 0.38 to 1.45, Si 0.30 to 1.40 in that window (file stores fractions; we read them as mass fractions, unit not confirmed in the file) | temper labels (T6 82 rows, T4 29, T5 8, others); condition text per row | **no property columns**; the paper states mechanical properties are available on request | Mendeley Data xlsx, CC BY, downloaded (sha256 5628a949...5f01) |
| Bhat, Barnard, Birbilis, J. Mater. Sci. 59 (2024) 1448, 10.1007/s10853-023-09317-2 | same 1154-alloy set with yield strength, tensile strength, elongation | as above | as above | YS, UTS, elongation (property values not public in the Mendeley file) | open article; Springer site bot-blocked |
| Tiwari, Jalalian, Mendis, Eskin, JOM 75 (2023) 4526, 10.1007/s11837-023-06025-9 | 50 extruded T6 6xxx alloys | Mg:Si ratio 0.34 to 1.71 (wt ratio as stated) | T6 only (temperature not stated per alloy) | YS 170 to 430 MPa, UTS 210 to 483 MPa, elongation, modulus of toughness | supplementary data file (Brunel repository copy of article read; supplement not downloaded) |
| Tiwari, Mendis, Eskin, J. Sustain. Metall. 11 (2025) 2323, 10.1007/s40831-025-01112-4 | 292 data sets on 42 wrought 6xxx grades | not extracted | T5, T6, T7 tempers | mechanical and service properties | supplementary material (not downloaded) |
| Poznak, Thole, Sanders, Metals 8 (2018) 309, 10.3390/met8050309 | several alloys varying Mg + Si and Mg/Si (number not extracted) | not extracted | not extracted | full hardness and conductivity ageing curves, with natural-ageing delay and cooling rate as variables | CC BY; MDPI and Michigan Tech copies bot-blocked |
| Alshibli, Gebril, Shuaeib, Diffus. Found. Mater. Appl. 38 (2025) 3, 10.4028/p-Om12I7 | 15 runs (Box-Behnken design, three factors at three levels) | Mg 0.4 to 1.0, Si 0.4 to 0.8, varied independently | 170, 190 (centre, inferred from design), 210 C | hardness, HRB (cast A6063) | CC BY; values per run not extracted |
| Khangholi, Javidani, Maltais, Chen, MATEC Web Conf. 326 (2020) 08002, 10.1051/matecconf/202032608002 | 3 | Mg 0.65 constant; Si 0.40, 0.66, 0.89 | 180 C | peak HV 114, 127, 140 (abstract); conductivity | open (UQAC repository PDF read) |
| Zheng, Luo, Bai, Wang, Yin, Metals 7 (2017) 387, 10.3390/met7100387 | 3 | Si 1.21 (alloy A) to 1.73 (alloy C); Mg and alloy B Si not extracted | 170 C | peak HV 123, 127, 129 (Section 3.1 text via snippet) | CC BY |
| Deng et al., Mater. Res. Express 12 (2025) 106508, 10.1088/2053-1591/ae0f7e | 5 (Al-xMg-0.5Si, x = 0.55 to 0.95) | Mg 0.55 to 0.95, Si 0.5 | ageing of drawn wire (300 C for 10 h reported as optimum) | tensile strength, conductivity | CC BY; not read beyond abstract |
| Marioara, Andersen, Zandbergen, Holmestad, Metall. Mater. Trans. A 36 (2005) 691, 10.1007/s11661-005-1001-7 | several alloys of varied Mg/Si (number not extracted) | not extracted | not extracted | hardness and precipitate statistics | paywalled |
| Gupta, Lloyd, Court, Mater. Sci. Eng. A 316 (2001) 11, 10.1016/s0921-5093(01)01247-3 | balanced and excess-Si alloys (number not extracted) | not extracted | not extracted | hardness / strength ageing curves | paywalled |
| Yamane, Kawasaki, Mima, J. Jpn. Inst. Light Met. 19 (1969) 305, 10.2464/jilm.19.7_305 | alloys with excess Si or excess Mg over Mg2Si (number not extracted) | not extracted | not extracted | specific heat and hardness | open (J-STAGE, Japanese) |
| Shen et al., Materials 18 (2025) 3663, 10.3390/ma18153663 | continuous composition gradients (kinetic diffusion multiple 6016/Al-3Si/Al-1.2Mg/Al) | gradient | diffusion annealed 530 C for 72 h; tested at room temperature | nanoindentation hardness and creep | CC BY; not an ageing study |
| Oien, Myhr, Ringen, Metals 16 (2026) 142, 10.3390/met16020142 | "existing datasets from a range of Al-Mg-Si alloys" (abstract) | not extracted | not extracted | tensile properties | CC BY; dataset access not checked |

## Candidate additions, ordered by manual digitising effort (assessed 4 Oct 2026)

1. **Add beta (Al3Mg2) and gamma (Al12Mg17) to the enthalpy receiver** from
   Table 4 of the CC BY review (10.1155/2014/704283), keeping both
   experimental values per phase as a band (beta -1.22 to -2.47, gamma -1.01
   to -3.80 kJ/mol-atoms), and add the Jund et al. Mg2Si range (-0.22 to
   -0.31 eV/atom) as a second Mg2Si source. No digitising needed; this
   gives three compounds, enough for an order check of Mg2Si versus the
   Al-Mg phases.
2. **Khangholi et al. (180 C, three Si levels at fixed Mg)** and **Zheng
   et al. (170 C, three Si levels)**: peak hardness is stated in the text,
   so six further alloys can be added with no digitising (time to peak
   would need the figures).
3. **Alshibli et al. Box-Behnken study**: 15 runs with Mg and Si varied
   independently at three ageing temperatures; if the design matrix table
   carries the measured HRB per run, this is the best single
   multi-composition source. Needs a browser download (publisher blocks
   scripted access).
4. **Tiwari et al. JOM supplement (50 T6 alloys)**: tabular, so no
   digitising, but only T6 strength of commercial grades (one ageing state,
   other elements vary); useful as a strength-versus-composition receiver.
5. **Poznak et al. Metals 8, 309**: full ageing curves for alloys varying
   Mg + Si and Mg/Si independently; requires digitising curves and a
   browser download.
6. **Paywalled primaries** (Brown and Pratt 1970; Feufel et al. 1997; Zhong
   et al. 2005; Marioara et al. 2005; Gupta et al. 2001) through the
   university library, to confirm the measurement method and temperature
   behind the values in step 1.
7. The 1154-alloy Mendeley file is not a receiver as published (no
   property columns); the property values would have to be requested from
   the corresponding author.

## Download links for the two sources without scripted access (4 Oct 2026)

Both are CC BY; scripted access returns HTTP 403, so they need a manual download.

| Source | DOI | Article page | PDF |
|---|---|---|---|
| Alshibli, Gebril, Shuaeib, "Investigating the Effect of Aging Temperature, Magnesium, and Silicon Content on the Hardness of A6063 Alloy", Diffus. Found. Mater. Appl. 38 (2025) 3-13 | https://doi.org/10.4028/p-Om12I7 | https://www.scientific.net/DFMA.38.3 | https://www.scientific.net/DFMA.38.3.pdf |
| Poznak, Thole, Sanders, "The Natural Aging Effect on Hardenability in Al-Mg-Si: A Complex Interaction between Composition and Heat Treatment Parameters", Metals 8 (2018) 309 | https://doi.org/10.3390/met8050309 | https://www.mdpi.com/2075-4701/8/5/309 | https://www.mdpi.com/2075-4701/8/5/309/pdf ; also https://digitalcommons.mtu.edu/michigantech-p/1931 and https://doaj.org/article/00591d413f8a489499d0c5c54f817009 |

## Use of these candidates (4 Oct 2026)

The beta, gamma and Mg2Si values of section A1 are used as a measured
receiver over all eight combinations of one measurement per compound:
`MEASURED_ENTHALPY_COMBINATIONS.md` (script `../../measured_enthalpy_combinations.py`).

Of the two sources without scripted access, the full text of Poznak et al.
2018 was obtained and that of Alshibli et al. 2025 was not; Alshibli is
not used. Peak hardness of the six Poznak alloys was read from Figs. 3
and 10 into `poznak2018_ageing.json` and evaluated by
`compute/almgsi_campaign/poznak_receiver_study.py` (`poznak_receiver.json`);
PEGKi Annex S4.13 and Fig. S-poznak report it.
