# Laboratory receivers for the computed Al-Mg-Si campaign

Two receivers built from laboratory measurements, in addition to the observed
precipitation sequence (`../analysis/lab_sequence_study.md`). Script:
`../../lab_receivers_study.py` (run with a Python that has ase and scipy). Metrics: `lab_receivers_metrics.json`.
References (BibTeX fetched from doi.org): `refs.bib`.

## Receiver 1: measured formation enthalpies

**Sources searched (all machine-readable, all downloaded):**

| Compilation | Access | Entries | Al-Mg-Si compounds found |
|---|---|---|---|
| Kim, Meschel, Nash, Chen, Sci. Data 4, 170162 (2017), doi:10.1038/sdata.2017.162; data doi:10.6084/m9.figshare.c.3822835.v1 | figshare file 9218491 (raw, with method and citation) and matminer `expt_formation_enthalpy` (figshare file 13464467) | 1276 | none (only Al-Fe-Si ternaries) |
| Wang, Kingsbury et al., Sci. Rep. 11 (2021), doi:10.1038/s41598-021-94550-5; matminer `expt_formation_enthalpy_kingsbury` (Ward et al., doi:10.1016/j.commatsci.2018.05.018) | figshare file 28232538, sha256 a2d2ced9... (matches matminer metadata) | 2135 | Mg2Si only (from NIST-JANAF); the rest are oxides or Fe-bearing |
| NIST-JANAF Thermochemical Tables, SRD 13, doi:10.18434/T42S31 | table Mg-034 (text file saved in `raw/`) | | Mg2Si (cr) |
| OQMD | site returned HTTP 502 on every attempt during data collection; one archived page (2015) for Mg17Al12 lists no experimental value | | not available |

The hits are saved in `raw/compilation_hits.json`. No measured formation
enthalpy for gamma (Mg17Al12), epsilon (Mg23Al30), beta (Al3Mg2), beta'',
beta', U1, U2 or any Al-Si compound was found in these open compilations, so
the observed Al-Mg phases could not be added to this receiver.

**Value.** Mg2Si(cr), standard enthalpy of formation at 298.15 K:
-77.822 kJ/mol of formula unit (JANAF Mg-034), i.e. -25.94 kJ/mol of atoms,
**-0.2689 eV/atom**. The Wang/Kingsbury compilation carries the same value
(-0.26886 eV/atom) with no uncertainty, and the machine-readable JANAF table
gives none, so `sem` is null.

**Mapping rule.** The compilation's likely MP id (mp-1367, antifluorite
Fm-3m) equals the candidate's MP id, so the item is `mp_1367_Mg2Si`. This is
also the lowest-PBE Mg2Si candidate. Under the alternative rule "lowest-energy
Mg2Si polymorph per source", MEAM and NEP89 would select other polymorphs
(recorded in the metrics file).

**Result (n = 1; Kendall tau is undefined, error only):**

| Source | computed, eV/atom | error vs measured, eV/atom |
|---|---|---|
| QE PBE, 0 K | -0.160 | +0.109 |
| QE HSE06, 0 K | -0.130 | +0.139 |
| MP archive (GGA) | -0.109 | +0.160 |
| MACE-MP-0, 0 K | -0.081 | +0.188 |
| MACE-MP-0 MD, 300 K | -0.082 | +0.186 |
| MEAM MD, 300 K | +0.044 | +0.312 |
| NEP89 MD, 300 K | +0.089 | +0.358 |

Every source underbinds Mg2Si; the DFT sources are closest. MEAM and NEP89
give the wrong sign. With one compound this receiver checks the absolute
scale of each source but cannot rank sources by ordering.

## Receiver 2: age hardening and DSC of alloys with different Mg/Si ratio

**Source.** Ma, Tang, Tu, Ye, Jiang, Zhan, Zhao, "Effect of a Trace Addition
of Sn on the Aging Behavior of Al-Mg-Si Alloy with a Different Mg/Si Ratio",
Materials 13 (2020) 913, doi:10.3390/ma13040913 (CC BY 4.0; full text fetched
from PubMed Central PMC7078789, saved in `raw/`). Six alloys at Mg + Si of
about 1.6 wt% with three Mg/Si ratios, each with and without about 0.1 wt%
Sn. Sheets solution treated at 570 C for 30 min, water quenched, aged at
**180 C**; Vickers hardness at 4.9 N (mean of five indents); DSC of
as-quenched samples at 10 C/min.

| Alloy | Mg, wt% | Si, wt% | Sn, wt% | Mg/Si (atomic) | peak HV | time to peak, h | beta'' DSC peak, C |
|---|---|---|---|---|---|---|---|
| H-M | 1.01 | 0.60 | 0.00 | 1.95 | 117.5 | 12 | 249.6 |
| B-M | 0.81 | 0.80 | 0.00 | 1.17 | 128 | 8 | 252.8 |
| L-M | 0.60 | 1.03 | 0.00 | 0.67 | 119.6 | 5 | 244.7 |
| H-M-Sn | 1.05 | 0.62 | 0.10 | 1.96 | 126 | 8 | 243.9 |
| B-M-Sn | 0.80 | 0.81 | 0.11 | 1.14 | 128 | 8 | 242.3 |
| L-M-Sn | 0.59 | 1.01 | 0.11 | 0.68 | 116 | 3 | 236.7 |

Compositions: Table 1 (at% computed with Al as balance; Fe and Ti also
tabulated, see `age_hardening.json`). Peak hardness and time to peak: Section
3.1 text (curves in Figure 1). DSC beta'' peak: Table 4. The beta' exotherm
appears only in a figure and was not digitised.

**Mapping rule (stated, simple).** Each alloy is an item. A source's
prediction for an alloy is its formation energy (eV/atom) of the most stable
Mg- and Si-containing candidate whose atomic Mg/Si ratio is closest, in log
ratio, to the alloy's. Observed phases beta'' (Mg5Si6) and beta' (Mg9Si5) are
included. Expected direction: a more negative predicted energy (stronger
driving force) should give higher peak hardness, shorter time to peak and a
lower beta'' DSC temperature. Agreement = -tau for hardness and +tau for time
and DSC temperature (+1 = measured order reproduced).

**Agreement (Kendall tau-b, sign applied):**

| Source | peak HV, n = 3 | time to peak, n = 3 | DSC beta'', n = 3 | peak HV, n = 6 | time to peak, n = 6 | DSC beta'', n = 6 |
|---|---|---|---|---|---|---|
| QE PBE, 0 K | -0.33 | -1.00 | -0.33 | +0.15 | -0.83 | -0.30 |
| QE HSE06, 0 K | -0.33 | -1.00 | -0.33 | +0.15 | -0.83 | -0.30 |
| MACE-MP-0 (0 K and 300 K MD) | -0.33 | -1.00 | -0.33 | +0.15 | -0.83 | -0.30 |
| MEAM MD, 300 K | +0.33 | -0.33 | -1.00 | +0.77 | -0.50 | -0.30 |
| NEP89 MD, 300 K | +0.33 | +1.00 | +0.33 | -0.15 | +0.83 | +0.30 |

n = 3 is the clean Sn-free series; n = 6 adds the Sn alloys, which share the
prediction of their Sn-free partner (Sn is not represented), so they add
ties rather than information. For PBE, HSE06 and MACE the mapping picks
Mg2Si (H-M), MgAlSi (B-M) and beta'' Mg5Si6 (L-M); their energy decreases
with Mg/Si ratio, while the measured time to peak rises with Mg/Si ratio, so
the DFT-like sources come out reversed. NEP89 places beta'' lowest and so
matches the kinetic order; this coincides with its incorrect energy for
equilibrium Mg2Si (Receiver 1), so the agreement is not evidence of better
physics.

## Limits

- Receiver 1 has a single compound: it measures absolute error, not order.
  Measured enthalpies of the Al-Mg phases exist in the calorimetry
  literature, but none were present in the open machine-readable
  compilations checked here, and OQMD was unreachable.
- Receiver 2 has three independent compositions (six alloys with Sn). With
  n = 3 the only possible tau values are -1, -1/3, +1/3 and +1, and none is
  statistically distinguishable from zero.
- Hardness and DSC peak temperatures reflect kinetics (vacancies, clustering,
  precipitate number density), whereas the mapping uses a 0 K or 300 K
  formation energy of one phase. The mapping is a stated cold-start default;
  a mapping through solute-cluster binding energies
  (`../analysis/cluster_binding.json`) or a composition-weighted precipitate
  fraction is a possible alternative and would change the taus.
- The natural-ageing delay before artificial ageing is not reported in the
  source paper.

## Further measured sources (literature search, 4 Oct 2026)

`LITERATURE_CANDIDATES.md` (references in `refs_candidates.bib`) lists
measured formation enthalpies for beta (Al3Mg2) and gamma (Al12Mg17) from two
calorimetric studies quoted in a CC BY review, a second Mg2Si source, and
multi-composition age-hardening sources. They are not included in Receiver 1:
the two studies disagree by a factor of 2 to 4 and order beta and gamma
oppositely (Brown and Pratt 1970: gamma more negative; Predel and Huelse 1978:
beta more negative), so no single measurement is preferred. All eight
combinations of one measurement per compound are evaluated separately in
`MEASURED_ENTHALPY_COMBINATIONS.md`.
