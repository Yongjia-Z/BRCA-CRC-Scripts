# BRCA-CRC-Scripts

Analysis code for:

> Zeng Y, Shi S, Wang F, Zhao L, Ye L, Wu X, Yan D. Subtype-specific Core Transcriptional Regulatory Circuitry Links Transcriptional States, Chromatin Context, and Functional Dependencies in Breast Cancer. *Functional & Integrative Genomics* (under review).

## Overview

The repository contains the code that analyses the public data and produces the results reported in the paper, including the analyses added during peer review. It does not download data, select datasets, draw figures, or format supplementary tables:

- **Public data.** Obtain the datasets from the sources listed in Tables S2 and S5e of the paper and place them under `datasets/<resource>/`. Each script names the files it reads.
- **Dataset selection.** The lists of datasets used are given in `data/dataset_lists/`: the ReMap breast tier-1 datasets, the ENCODE MCF-7 H3K27ac files, and the curated ReMap, ChIP-Atlas, and GEO profiles (Tables S5c–S5e).
- **Results.** Result tables are in Online Resource 1 of the paper.

## Layout

| Path | Content |
|---|---|
| `pipeline/00_data/` | Preprocessing of public data and the evidence layers of TF scoring (TCGA, ATAC, H3K27ac, ReMap, scATAC, external cohorts) |
| `pipeline/01_tf_scoring/` | Evidence scoring of 1,639 TFs and selection of the 13 CRC TFs (Fig. 1A, Table S1); TF–context association (Fig. 1E) |
| `pipeline/02_chromatin/` | Local ATAC–RNA coupling (Fig. 3D, Table S6) |
| `pipeline/03_bulk_modules/` | TCGA-BRCA module scores, tumor-state signatures, and MCP-counter correlations (Fig. 6) |
| `pipeline/04_single_cell/` | GSE176078 malignant-cell states (Fig. 7) |
| `pipeline/05_outcome_dependency/` | State scores, Cox models, external cohorts, and DepMap (Fig. 8) |
| `revision_analyses/` | Analyses added during peer review (A0–A22, first round R1) |
| `_repro_common/` | Helper functions |
| `data/` | Inputs that public data cannot regenerate: the context membership that defined the modules of Figs. 6–8 (`crc_tf_table_modules.tsv`), the network of the submitted manuscript (`submitted_network/`), and the dataset lists (`dataset_lists/`) |
| `source_tables/` | Panel source tables of the original analysis read by the review analyses |

Scripts read public data from `CRC_DATA_ROOT` (default `<repo>/datasets`) and write processed tables to `CRC_CORE_ROOT` (default `<repo>/core_data`) and `source_tables/`. The review analyses take the project root (the repository) as their first argument.

## Run order

1. `pipeline/00_data`, scripts numbered within each folder:
   1. `tf_resources`, `gencode`;
   2. `tcga` 03–08;
   3. `cistrome_cancer` 03, then `tcga` 09;
   4. `cistrome_cancer` 05 and `encode` 01, then `tcga` 10–11;
   5. `encode` 02 and `sedb` 01, then `tcga` 12 and `trrust` 00;
   6. `tcga` 14 and 17;
   7. `gse85158_h3k27ac`, `scatac_multiome`, `remap`;
   8. `gse96058`, `gse176078`, `immune_deconv`.
2. TF scoring and Fig. 1E:
   ```bash
   python pipeline/01_tf_scoring/tf_evidence_scoring.py --raw-root datasets --remap-result-dir core_data/remap_reanalysis/results \
       --output-dir source_tables/figure1 --intermediate-dir source_tables/intermediate/figure1a
   python pipeline/01_tf_scoring/prepare_subtype_association_inputs.py
   python pipeline/01_tf_scoring/subtype_association.py --input-dir source_tables/figure1/subtype_association_inputs --output-dir source_tables/figure1
   ```
3. `02_chromatin`, then `03_bulk_modules`.
4. `04_single_cell` 01 → 02 → 03.
5. `05_outcome_dependency`: `tcga_state_scores.py`, `external_cohort_state_genes.py`, `depmap_state_scores.py`, `tcga_state_cox.py`.
6. Review analyses:
   1. `R1_first_round`;
   2. A7;
   3. A8, A15, A18, A21;
   4. A14 (02–06, 09);
   5. A16 (`scripts/`, in the order given below);
   6. A19, then A17, A22, A13;
   7. A10, then A1 with A10's `A10_fig6_scores_corrected_membership.tsv` as third argument, then A0, A12.

Files placed by hand:
- GSE251868 HER2 H3K27ac BED files: `datasets/GSE251868_HER2_H3K27ac/data`.
- GSE245830 TNBC H3K27ac (`GSE245830_RAW.tar`): `datasets/GSE245830_TNBC_H3K27ac/data`. A7 takes the folder into which it is extracted as second argument.
- GSE213840 FOXC1 peaks: `datasets/GSE213840_FOXC1_TNBC_ChIPseq/`. TNBC binding is read from folders named `GSE*_<TF>_TNBC*`.
- Aran et al. (2015) Supplementary Data 1: `revision_analyses/A1_purity_adjustment/Aran2015_SupData1.xlsx`.
- ChIP-Atlas `experimentList.tab`: `revision_analyses/A16_supplementary_ChIP/raw/metadata/` (read by A15 and A18).
- Peak files of the curated ChIP-Atlas PGR profiles: `revision_analyses/A14_occupancy_and_network_checks/raw_chipatlas_PGR/<SRX>.05.bed`.
- A16 inputs: peak files, GEO SOFT records, and chain files under `revision_analyses/A16_supplementary_ChIP/raw/` and `derived/`, at the paths listed in `data/dataset_lists/A16_peak_file_manifest.tsv`; the FOXC1 ChIP and input bigWig files of GSE213840 (GSM6595097–GSM6595101) in `datasets/Figure4_browser_tracks/data`. The A16 scripts read the project root from `CRC_PROJECT_ROOT` (default: the repository root).

## Scripts

Scripts without arguments read their inputs from `CRC_DATA_ROOT` and `CRC_CORE_ROOT`.

### Helper functions (`_repro_common/`)

| Script | Purpose | Arguments |
|---|---|---|
| `crc_repro.py` | Shared helpers (table input/output, PAM50-to-context mapping, BH correction, gene-by-sample matrices) used by subtype_association.py and A12. | imported by other scripts |

### Public data preprocessing (`pipeline/00_data/`)

| Script | Purpose | Arguments |
|---|---|---|
| `cistrome_cancer/03_candidate_tf_tf_edges.py` | Builds candidate TF-to-TF edges from the Cistrome Cancer CancerTarget target tables. | – |
| `cistrome_cancer/05_enhancer_signal_near_tf_genes.py` | Summarizes the TCGA enhancer signal in the gene windows of the candidate TFs. | – |
| `encode/01_mcf7_h3k27ac_overlap.py` | Counts overlaps of the ENCODE MCF-7 H3K27ac peaks (file list in data/dataset_lists) with the candidate TF gene windows. | – |
| `encode/02_mcf7_rose_like_super_enhancers.py` | Stitches MCF-7 H3K27ac peaks in the manner of ROSE to obtain a super-enhancer proxy and links it to the candidate TF loci. | – |
| `gencode/00_parse_gencode_v36.py` | Parses the GENCODE v36 annotation into a gene table. | – |
| `gse176078/00_convert_mtx_to_h5ad.py` | Converts the GSE176078 count matrix and metadata into an h5ad file. | – |
| `gse85158_h3k27ac/02_quantify_hg19_to_hg38.py` | Quantifies GSE85158 H3K27ac signal in TF gene windows (GENCODE v19), lifts the windows to hg38 and summarizes by subtype. | `[--limit-files N] [--max-workers N]` |
| `gse96058/00_annotation_and_module_scores.py` | Reads the GSE96058 sample annotation and expression and computes CRC module scores (used by A18 and R1). | – |
| `immune_deconv/00_mcp_counter.R` | MCP-counter immune and stromal estimates for TCGA-BRCA (input of Fig. 6C, D). Requires the MCPcounter package. | – |
| `remap/02_breast_promoter_overlap.py` | Counts peak overlaps of each ReMap breast tier-1 dataset (list in data/dataset_lists) with the candidate TF promoters. | – |
| `remap/03_tnbc_direct_binding_overlap.py` | Overlaps the TNBC FOXC1 and TP63 ChIP-seq/CUT&Tag peaks with the candidate TF promoters and ±100 kb loci (TNBC part of Fig. 4B, C). | – |
| `scatac_multiome/00_tcga_brca_scatac_locus_support.py` | Counts accessible peaks and motifs at the CRC TF loci in breast cancer and healthy breast scATAC peak sets (Fig. 3A–C; scATAC layer). | – |
| `sedb/01_dbsuper_breast_super_enhancers.py` | Lifts the dbSUPER breast super-enhancers to hg38 and counts overlaps within ±100 kb of the TF loci. | – |
| `tcga/03_atac_peak_summary.py` | Summarizes sample- and peak-level accessibility of the ATAC peak matrix. | – |
| `tcga/04_atac_footprinting.py` | Reads the TCGA ATAC footprinting and motif tables and correlates footprint depth with TF expression (footprint layer). | – |
| `tcga/06_tf_expression_activity_screen.py` | Builds the TF expression matrix and sample annotation, tests expression and subtype association, and computes regulon-based TF activity from DoRothEA (RNA and activity layers; matrices of Fig. 1E, F). | – |
| `tcga/07_atac_peak_annotation_subtype.py` | Annotates the ATAC peaks and tests for subtype-specific peaks. | – |
| `tcga/08_atac_candidate_tf_integration.py` | Adds subtype-specific peaks, peak-to-gene links and motif enrichment to the candidate TF evidence table. | – |
| `tcga/09_cistrome_tf_tf_edges.py` | Adds the Cistrome Cancer TF-TF edges to the candidate TF evidence table. | – |
| `tcga/10_enhancer_h3k27ac_evidence.py` | Adds the TCGA enhancer and ENCODE MCF-7 H3K27ac evidence. | – |
| `tcga/11_crc_member_evidence_scores.py` | Computes CRC member evidence scores of the candidate TFs and lists the high-confidence candidates. | – |
| `tcga/12_tf_tf_edge_evidence_scores.py` | Scores the evidence of every candidate TF-TF edge (edge evidence table of the Fig. 2 network). | – |
| `tcga/14_crc_activity_scores.py` | Computes per-sample expression and activity module scores of the CRC TFs. | – |
| `tcga/17_clinical_outcome_dataset.py` | Merges the CRC module scores with TCGA-BRCA clinical follow-up and runs clinical association and survival (log-rank, Cox) analyses. | – |
| `tf_resources/00_lambert_tf_list.py` | Extracts the 1,639 human TFs from the Lambert et al. (2018) database table. | – |
| `tf_resources/01_convert_dorothea_rda.R` | Converts the DoRothEA regulon RDA file to a TSV table. | `[input.rda] [output.tsv]` |
| `trrust/00_trrust_edges.py` | Adds the TRRUST v2 human TF-target relations to the TF-TF edge evidence table. | – |

### TF scoring and Fig. 1E (`pipeline/01_tf_scoring/`)

| Script | Purpose | Arguments |
|---|---|---|
| `tf_evidence_scoring.py` | Fig. 1A and Table S1: scores 1,639 TFs on seven core and three supportive layers, selects the 13 CRC TFs by the tier rules and assigns context membership. | `--raw-root --remap-result-dir --output-dir --intermediate-dir [--max-workers ...] (see README)` |
| `prepare_subtype_association_inputs.py` | Writes the three inputs of subtype_association.py (TF expression matrix, sample annotation, CRC TF table). | `[out_dir]` |
| `subtype_association.py` | Fig. 1E: Kruskal–Wallis test of each of 1,631 TFs across the Luminal, HER2, Basal/TNBC and Normal-like groups with BH correction; writes all TFs and the top 20. | `--input-dir DIR --output-dir DIR [--top-n 20]` |

### Chromatin (`pipeline/02_chromatin/`)

| Script | Purpose | Arguments |
|---|---|---|
| `atac_rna_coupling.py` | Fig. 3D and Table S6: correlation of ATAC accessibility within ±50 kb of each TF locus with the TF's expression. | – |

### Bulk modules (`pipeline/03_bulk_modules/`)

| Script | Purpose | Arguments |
|---|---|---|
| `module_scores_tumor_state_mcp.py` | Fig. 6: CRC module scores (expression and activity) in TCGA-BRCA and their correlations with tumor-state signatures and MCP-counter estimates, before and after context adjustment. | – |

### Single cell (`pipeline/04_single_cell/`)

| Script | Purpose | Arguments |
|---|---|---|
| `01_process_and_score_cells.py` | Fig. 7A: processes GSE176078, computes the UMAP embedding and per-cell CRC TF module scores. | – |
| `02_state_discovery.py` | Fig. 7B, C: clusters TF co-expression within ER+, HER2+ and TNBC malignant cells and defines the six malignant states and their signature genes. | – |
| `03_state_scores_per_cell.py` | Fig. 7C: per-cell state scores and dominant state. | – |

### Outcome and dependency (`pipeline/05_outcome_dependency/`)

| Script | Purpose | Arguments |
|---|---|---|
| `tcga_state_scores.py` | Fig. 8A: the six state scores in TCGA-BRCA and Kruskal–Wallis tests across contexts. | – |
| `external_cohort_state_genes.py` | Projects the state signatures onto TCGA, METABRIC and GSE96058 and fits continuous Cox models in each cohort. | – |
| `depmap_state_scores.py` | Fig. 8C: state expression scores in DepMap breast cancer lines and their correlation with CRISPR gene effect. | – |
| `tcga_state_cox.py` | Fig. 8B and Table S7: continuous Cox models of each state score in TCGA-BRCA (OS, PFS; per 1 SD; 120-month truncation), via tcga_state_cox.R. | `[out_dir]` |
| `tcga_state_cox.R` | Called by tcga_state_cox.py: Cox models and proportional-hazards test with the survival package. | `<model input tsv> <results tsv>` |

### Review analyses (`revision_analyses/`)

| Script | Purpose | Arguments |
|---|---|---|
| `R1_first_round/R1_analyze_revision.py` | First review round: selection evidence of XBP1 and ZNF703, layer omission, weight perturbation and random sets, Fig. 3D peak audit, Fig. 6A/8A medians, self-locus occupancy, Cox models with BH correction and clinical adjustment, DepMap matched null and permutations (tables R1_1–R1_6). | `[project_root] [out_dir]` |
| `R1_first_round/R1_external_projection.py` | First review round: projection of the six state signatures onto METABRIC and GSE96058 (tables R1_7; Fig. S4, Table S10). | `[project_root] [out_dir]` |
| `A0_reanalysis_checks/A0_verify.py` | Independent re-computation of the review analyses: selection rules, random sets, Fig. 3D statistics, Cox models, DepMap matched null, external projection, Fig. 6A/8A medians. | `<project_root> <out_dir>` |
| `A1_purity_adjustment/A1_purity_adjusted_immune.py` | Fig. S3 and Tables S9a, S9b: purity-adjusted (CPE) correlations of module scores with MCP-counter estimates. | `<project_root> <out_dir> [scores_tsv]` |
| `A7_Basal_TNBC_H3K27ac_layer/A7_basal_h3k27ac_audit.py` | Adds the GSE245830 TNBC H3K27ac peaks and recomputes the Basal/TNBC layer and context membership (Tables S1, S4d, S4e). | `<project_root> <extracted_GSE245830_dir> <out_dir>` |
| `A8_selection_robustness/A8_selection_robustness_corrected.py` | Selection robustness on the revised score table: weight perturbation, layer omission, threshold shifts, random sets (Fig. S1; Tables S4a–c). | `<A7_corrected_csv> <out_dir>` |
| `A10_Fig6_modules/A10_fig6_modules_corrected.py` | Fig. 6 module scores with the revised Luminal membership (without RUNX1) (Tables S9c, S9d). | `<project_root> <A7_corrected_csv> <out_dir>` |
| `A12_Fig3D_mean_accessibility_script/figure3d_corrected.py` | Corrected Fig. 3D statistics: mean local accessibility versus expression, and BH correction over all local peak tests. | `--input-dir DIR --output-dir DIR [--flank-bp 50000]` |
| `A13_motif_context_membership/A13_motif_membership.py` | Classifies the Fig. 2C motifs by the context membership of their nodes (Table S14). | `<figure2c.csv> <A7_CRC_TF_calls_stored_vs_corrected.tsv> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_02_pgr_target_and_network.py` | Adds PGR as an occupancy target, rebuilds the edge set with ESR1→PGR and enumerates motifs. | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_03_chipatlas_crosscheck.py` | Cross-checks promoter occupancy of the 13 TFs with local ChIP-Atlas breast peak files. | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_04_knockdown_label_sensitivity.py` | Effect of removing five tier-1 datasets with knockdown labels missed by the metadata pattern. | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_05_tables_and_control.py` | Updated occupancy matrix (PGR target column from ReMap, PGR row from ChIP-Atlas), self-occupancy and network-rule tables. | `<project_root> <S5_xlsx> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_06_edges_table.py` | Tidies the table of edges scoring ≥ 6 (Table S14b). | `<out_dir>` |
| `A14_occupancy_and_network_checks/A14_09_pgr_row_binding_knockdown.py` | Final occupancy tables from the curated ReMap and ChIP-Atlas PGR datasets (curation in data/dataset_lists): PGR column, row and self-locus value, and the effect of the PGR binding score on selection. | `<project_root> <Table S5 workbook> <out_dir>` |
| `A15_selection_coverage_audit/A15_selection_bias_audit.py` | Read-only audit of the selection: layer coverage, sparse-layer membership, ranking without sparse layers, uniform footprint evidence, rule replay under alternative scoring, TFs with breast ChIP but a zero binding layer. | `<project_root> <revision_analyses dir> <out_dir> <experimentList.tab>` |
| `A16_supplementary_ChIP/scripts/a16_common.py` | Shared A16 paths, file input/output and checksums, promoter table and peak-overlap helpers. | imported by other scripts |
| `A16_supplementary_ChIP/scripts/a16_liftover.py` | (order 1) Chain-based interval lift-over (≥ 95% of bases in one chain), cross-checked with pyliftover. | – |
| `A16_supplementary_ChIP/scripts/a16_bootstrap.py` | (order 2) Reproduces the PGR baseline of the submitted manuscript and writes the promoter windows. | – |
| `A16_supplementary_ChIP/scripts/a16_analyze_peaks.py` | (order 3) Reproduces the earlier occupancy matrix and analyses the supplementary ChIP peaks by source (datasets in data/dataset_lists). | – |
| `A16_supplementary_ChIP/scripts/a16_signal_quantify.mjs` | (order 4) Quantifies the FOXC1 genome-wide bigWig signal at the TF promoters (no peak calling). Requires Node.js with @gmod/bbi. | `node a16_signal_quantify.mjs <project_root>` |
| `A16_supplementary_ChIP/scripts/a16_quality.py` | (order 5) Interval audit of the peak files and correspondence between processing pipelines. | – |
| `A16_supplementary_ChIP/scripts/a16_sensitivity.py` | (order 6) Conditional replay on the original 13-TF inputs, changing only the binding layer (Table S4g). | – |
| `A17_score5_motifs/A17_01_score5_motifs.py` | Whether the HER2 or Basal/TNBC sets form closed loops once edges scoring 5 are admitted (Tables S14d, S14e). | `<A16 results dir> <out_dir>` |
| `A18_selection_and_subtype_checks/A18_selection_subtype_checks.py` | Checks for comments R1.5 and R5.2–R5.5: Fig. 1E ranks and reasons, five-layer reference, balanced subtype composition, PAM50 versus receptor groups, reused resources (Tables S2c, S2d, S4h, S4k, S4l). | `<project_root> <revision_analyses dir> <A16 dir> <out_dir>` |
| `A19_occupancy_network_self_locus/A19_occupancy_network.py` | Occupancy matrix, TF network and motifs with the supplementary ChIP, and replication in studies not used for selection (Figs. 2, 4C; Tables S5, S14). | `<project_root> <revision_analyses dir> <A16 dir> <out_dir>` |
| `A19_occupancy_network_self_locus/A19_02_self_locus_100kb.py` | Self-locus occupancy within TSS ±100 kb, ranked against the same window of all protein-coding genes (Fig. S2B; Table S5h). | `<project_root> <A16 dir> <out_dir>` |
| `A21_RNA_grouping_subsampling/A21_rna_grouping_balanced.py` | Grouping of the RNA-layer subtype test and its recalculation in balanced and proportional subsamples (Tables S4i, S4i2). | `<project_root> <revision_analyses dir> <out_dir>` |
| `A22_binding_by_cell_line/A22_cell_line_binding.py` | Cell lines behind the binding of each network edge, and closed structures within one cell line (Tables S14f, S14g). | `<A16 dir> <A19 out dir> <out_dir>` |

## Requirements

- Python ≥ 3.10 with the packages in `requirements.txt` (`pyranges` is optional). The runs in the paper used Python 3.12–3.13; see `software_versions.tsv`.
- R with `survival` (Cox models) and `MCPcounter` (R 4.5.2 in the paper).
- Node.js with `@gmod/bbi` for the bigWig step of A16.

## Notes on interpretation

- H3K27ac signals are enhancer and super-enhancer proxies; no formal super-enhancer calling (for example ROSE) was performed.
- Hi-C and other contact data show physical proximity and are not used as TF–TF regulatory edges.
- DoRothEA-derived values are regulon-based activity proxies.
- Survival analyses are association tests, not risk models.

## License

MIT (see `LICENSE`). Third-party datasets and software are subject to their own terms and should be cited as in the paper.
