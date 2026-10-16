# PEGKi

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23413455.svg)](https://doi.org/10.5281/zenodo.23413455)

**Paper.** Evaluating transfer across modalities in autonomous materials systems (submitted to Digital Discovery (RSC)).

This repository holds the code and data behind the paper: the scripts that draw its figures, the result files they read, the code that produced those results, and the archived or computed inputs that code reads. The manuscript itself is not part of the repository. Every script runs unchanged from the repository root.

**Authors.** Kinston Ackölf, Taylor D. Sparks, Materials Science and Engineering, The University of Utah.

**Citation.** `CITATION.cff` (software citation) and `.zenodo.json` (archive metadata).

## Layout

- `shared/`, `compute/`, `campaign_protocols/`: result files and the code that produced them.
- `writings/papers_2026/regen_figs/` and `paper_fig_style.py`: figure scripts (each reads only committed result files and draws at print size into `writings/papers_2026/<paper>/figures/`).

## Running

Figures: run a script from the folder root, for example `python writings/papers_2026/regen_figs/<script>.py` (environment: `writings/papers_2026/regen_figs/requirements.txt`). `REPRODUCE.md` gives the full procedure: software versions, the rebuild of the band-gap archive, the commands that regenerate every result file, and the DFT, MD and potential calculations behind the computed records.

## Figure scripts

| Figure | Script | Files read |
|---|---|---|
| `source_combinability.png` | `writings/papers_2026/regen_figs/pegki_source_combinability.py` | `shared/output/source_combinability/source_combinability.json` |
| `ladder.png` | `writings/papers_2026/regen_figs/pegki_ladder.py` | `shared/output/hypergraph_trial_bandgap/hypergraph_trial_bandgap.json`, `shared/output/nonadditive_pooled/nonadditive_pooled.json` |
| `family_offsets.png` | `writings/papers_2026/regen_figs/pegki_family_offsets.py` | `shared/output/family_offsets/family_offsets.json` |
| `omega_regret.png` | `writings/papers_2026/regen_figs/pegki_omega_regret.py` | `shared/output/omega_regret/omega_regret.json` |
| `set_operations.png` | `writings/papers_2026/regen_figs/pegki_set_operations.py` | `shared/output/set_operation_study/set_operation_study.json` |
| `routing_k_sweep.png` | `writings/papers_2026/regen_figs/pegki_routing_k_sweep.py` | `shared/output/routing_k_sweep_50seed/routing_k_sweep.json`, `shared/output/routing_k_sweep_50seed/routing_k_sweep_bootstrap.json` |
| `decision_body_cycles.png` | `writings/papers_2026/regen_figs/pegki_decision_body_cycles.py` | `shared/output/decision_body_cycles/decision_body_cycles.json` |
| `live_transfer.png` | `writings/papers_2026/regen_figs/live_almgsi_figs.py` | `compute/almgsi_campaign/results/analysis/source_agreement.json`, `compute/almgsi_campaign/results/campaigns/qe_hse06_all/summary.json`, `compute/almgsi_campaign/results/campaigns/qe_pbe_all/summary.json`, `compute/almgsi_campaign/results_T/analysis/temperature_study.json` |
| `bridge_budget.png` | `writings/papers_2026/regen_figs/pegki_bridge_budget.py` | `shared/output/bridge_budget/bridge_budget.json` |
| `pair_combinatorics.png` | `writings/papers_2026/regen_figs/pegki_pair_combinatorics.py` | `shared/output/pair_combinatorics/pair_combinatorics.json` |
| `local_cycles.png` | `writings/papers_2026/regen_figs/pegki_local_cycles.py` | `shared/output/campaign_loop_bands/campaign_loop_bands.json`, `shared/output/decision_body_cycles/decision_body_cycles.json` |
| `lab_sequence.png` | `writings/papers_2026/regen_figs/pegki_lab_sequence.py` | `compute/almgsi_campaign/results/analysis/lab_sequence_study.json` |
| `cluster_transfer.png` | `writings/papers_2026/regen_figs/pegki_cluster_transfer.py` | `compute/almgsi_campaign/results/analysis/cluster_transfer.json` |
| `measured_enthalpy.png` | `writings/papers_2026/regen_figs/pegki_measured_enthalpy.py` | `compute/almgsi_campaign/results/lab_receivers/measured_enthalpy_combinations.json` |
| `receiver_matrix.png` | `writings/papers_2026/regen_figs/pegki_receiver_matrix.py` | `compute/almgsi_campaign/candidates_index.json`, `compute/almgsi_campaign/results/analysis/hull_receiver.json`, `compute/almgsi_campaign/results/analysis/vacancy_binding_n3.json`, `compute/almgsi_campaign/results/calphad/calphad_T300.json` and 9 more |
| `poznak_ageing.png` | `writings/papers_2026/regen_figs/pegki_poznak_ageing.py` | `compute/almgsi_campaign/results/lab_receivers/poznak2018_ageing.json`, `compute/almgsi_campaign/results/lab_receivers/poznak_receiver.json` |
| `almgsi_campaign.png` | `writings/papers_2026/regen_figs/pegki_almgsi_campaign.py` | `compute/almgsi_campaign/results/analysis/prisms_comparison/prisms_comparison_qe_hse06.json`, `compute/almgsi_campaign/results/analysis/prisms_comparison/prisms_comparison_qe_pbe.json`, `compute/almgsi_campaign/results/campaigns/gpaw_r2scan_all/summary.json`, `compute/almgsi_campaign/results/campaigns/qe_hse06_all/summary.json` and 1 more |

## Producing code

The result files above were written by these scripts (local modules they import are included as well):

- `shared/scripts/source_combinability.py`
- `shared/scripts/hypergraph_trial_bandgap.py`
- `shared/scripts/nonadditive_pooled.py`
- `shared/scripts/family_offsets.py`
- `shared/scripts/omega_regret.py`
- `shared/scripts/set_operation_study.py`
- `shared/scripts/routing_k_sweep.py`
- `shared/scripts/decision_body_cycles.py`
- `shared/scripts/bridge_budget.py`
- `shared/scripts/pair_combinatorics.py`
- `shared/scripts/campaign_loop_bands.py`
- `shared/scripts/local_cycles.py`
- `shared/scripts/learned_transfer_weights.py`
- `shared/scripts/transfer_score_checks.py`
- `shared/scripts/property_family_receivers.py`
- `shared/scripts/materials_summaries.py`
- `campaign_protocols/code/experiments/24_hyperedge_budget.py`
- `campaign_protocols/code/experiments/26_paired_campaign.py`
- `campaign_protocols/code/experiments/27_learned_acquisition.py`
- `compute/bandgap_archive/build_inventory.py`
- `compute/bandgap_archive/make_splits.py`
- `compute/bandgap_archive/refetch_channels.py`
- `compute/almgsi_campaign/generate_candidates.py`
- `compute/almgsi_campaign/add_observed_phases.py`
- `compute/almgsi_campaign/dgx_sources.py`
- `compute/almgsi_campaign/aggregate_sources.py`
- `compute/almgsi_campaign/analyze_live.py`
- `compute/almgsi_campaign/campaign.py`
- `compute/almgsi_campaign/live_campaign.py`
- `compute/almgsi_campaign/compare_live_replay.py`
- `compute/almgsi_campaign/temperature_study.py`
- `compute/almgsi_campaign/dft_receiver_study.py`
- `compute/almgsi_campaign/lab_sequence_study.py`
- `compute/almgsi_campaign/lab_receivers_study.py`
- `compute/almgsi_campaign/measured_enthalpy_combinations.py`
- `compute/almgsi_campaign/cluster_transfer_study.py`
- `compute/almgsi_campaign/hull_receiver_study.py`
- `compute/almgsi_campaign/vacancy_study.py`
- `compute/almgsi_campaign/vacancy_size_study.py`
- `compute/almgsi_campaign/qha_study.py`
- `compute/almgsi_campaign/calphad/calphad_receiver.py`
- `compute/almgsi_campaign/receiver_extension_study.py`
- `compute/almgsi_campaign/r2scan_receiver_study.py`
- `compute/almgsi_campaign/mp_r2scan_study.py`
- `compute/almgsi_campaign/fetch_mp_r2scan.py`
- `compute/almgsi_campaign/poznak_receiver_study.py`
- `compute/semiconductor_campaign/generate_pool.py`
- `compute/semiconductor_campaign/archive_donors.py`
- `compute/semiconductor_campaign/campaign_gap.py`
- `compute/semiconductor_campaign/summarise.py`
- `compute/alsicu_campaign/generate_candidates.py`
- `compute/alsicu_campaign/dgx_sources.py`
- `compute/alsicu_campaign/aggregate_sources.py`
- `compute/alsicu_campaign/campaign.py`
- `compute/almgsi_campaign/prisms_comparison.py`
- `compute/almgsi_campaign/prisms_sensitivity.py`
- `compute/almgsi_archive/audit_almgsi_coverage.py`
- `compute/almgsi_campaign/rerun_md_noise.py`
- `compute/almgsi_campaign/temperature_coalitions.py`
- `compute/almgsi_campaign/shard_cids.py`
- `compute/almgsi_campaign/interface_study.py`
- `compute/almgsi_campaign/migration_study.py`
- `compute/almgsi_campaign/r2scan_cluster_study.py`

## Carried-over code and data

Some results in the paper come from campaign protocols and analyses that also serve other projects of the group. The scripts and the result files they need are copied here, so this repository depends on no other repository. Their files are marked `(carried over)` in `CONTENTS.tsv`. Folders:

- `campaign_protocols/code/`
- `campaign_protocols/results/`

## Notes on inputs

- The version-locked band-gap archive (`shared/output/bandgap_archive_build/bandgap_universe.jsonl`, 261,053 records) is not stored in git; rebuild it from the included source caches with `python compute/bandgap_archive/build_inventory.py` and check it against the manifest hash in `build_stats.json`.
- Computed Al-Mg-Si records (DFT, MACE-MP-0, MD) are under `compute/almgsi_campaign/results*`; the HSE06 values with the half q grid, set aside after the convergence test, are left out.

## Contents by role

- code: 152 files
- input data: 3456 files
- figure data: 5522 files
- result named in paper: 19 files
- environment: 4 files
- code (carried over): 5 files
- result named in paper (carried over): 3 files
- result read by code (carried over): 21 files

`FILES.sha256` lists every file with its SHA-256 checksum; `CONTENTS.tsv` gives the role of each file.

## Provenance

The folder was assembled from the `singular_obj/` project of the Sparks group repository by `build_paper_repos.py`, which keeps the project's relative layout and gives working folders their public names.

## Licence

Code (Python, shell and SLURM scripts) is released under the MIT License (LICENSE). Data, result files and figures are released under the Creative Commons Attribution 4.0 International licence (LICENSE-CC-BY-4.0.txt). The archived third-party datasets (Materials Project, JARVIS-DFT, SNUMAT, the Zhuo et al. experimental compilation) keep their own terms; their sources are cited in the paper.

## Acknowledgements

Taylor D. Sparks acknowledges support from the Army Research Office Materials Design, under Contract No. W911NF-23-1-0333.
