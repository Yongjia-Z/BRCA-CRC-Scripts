# BRCA-CRC-Scripts

[English](#english) | [简体中文](#chinese)

---

<a id="english"></a>

## Overview

This repository provides the computational analysis code accompanying:

> Zeng Y, Shi S, Wang F, Zhao L, Ye L, Wu X, Yan D. **Subtype-specific Core Transcriptional Regulatory Circuitry Links Transcriptional States, Chromatin Context, and Functional Dependencies in Breast Cancer.** *Functional & Integrative Genomics* (under review).

The workflow integrates multi-omics public breast cancer resources—bulk RNA-seq, DoRothEA regulon activity, ATAC-seq accessibility/motifs/footprints, subtype-matched H3K27ac enhancer proxies, breast-context direct TF binding (ReMap, ChIP-Atlas, GEO), scATAC-seq, scRNA-seq, external expression cohorts (METABRIC, GSE96058), and DepMap CRISPR dependency screens—to score 1,639 human transcription factors (TFs), prioritize 13 breast cancer core regulatory circuitry (CRC) TFs across Luminal, HER2, and Basal/TNBC contexts, resolve six single-cell malignant transcriptional states, and evaluate their clinical and functional relevance. All additional analyses performed during peer review (`R1` and `A0–A22`) are included.

### Scope

- **Included:** End-to-end data preprocessing, multi-layer TF evidence scoring, chromatin coupling, bulk module/immune correlations, single-cell state discovery, survival and DepMap dependency modeling, peer-review revision analyses, curated dataset manifests (`data/dataset_lists/`), non-regenerable reference inputs (`data/`), and original panel source tables (`source_tables/`).
- **Not included:** Automated dataset downloads, raw/controlled-access public data matrices, figure plotting/vector styling scripts, and formatted supplementary Excel tables (published as **Online Resource 1** of the paper).

---

## Repository Structure

```text
BRCA-CRC-Scripts/
├── pipeline/
│   ├── 00_data/                 # Public dataset preprocessing & TF scoring evidence layers (27 scripts)
│   ├── 01_tf_scoring/           # 1,639-TF evidence scoring, 13 CRC TFs & subtype tests (3 scripts)
│   ├── 02_chromatin/            # Local ±50 kb ATAC–RNA coupling (1 script)
│   ├── 03_bulk_modules/         # TCGA-BRCA module scores, tumor states & MCP-counter (1 script)
│   ├── 04_single_cell/          # GSE176078 processing & 6 malignant-cell states (3 scripts)
│   └── 05_outcome_dependency/   # TCGA/external cohort state scores, Cox models & DepMap (5 scripts)
├── revision_analyses/           # Peer-review sensitivity, audit & supplementary ChIP analyses (29 scripts)
├── _repro_common/               # Shared I/O, PAM50-to-context mapping & BH correction helpers (1 script)
├── data/
│   ├── crc_tf_table_modules.tsv # Context membership defining the modules in Figs. 6–8
│   ├── dataset_lists/           # Curated ReMap, ENCODE, ChIP-Atlas & GEO manifests (Tables S5c–S5e)
│   └── submitted_network/       # Submitted manuscript network edges & first-round occupancy table
├── source_tables/               # Original panel source tables read by peer-review analyses
├── requirements.txt             # Python package requirements
├── software_versions.tsv        # Exact Python, R & Node.js versions audited for the paper
├── CITATION.cff                 # Machine-readable citation metadata
└── LICENSE                      # MIT License
```

| Directory | Paper Mapping | Summary |
|---|---|---|
| `pipeline/00_data/` | Upstream layers (Figs. 1–4, 6) | Preprocessing of public resources (TCGA, ATAC, H3K27ac, ReMap, scATAC, external cohorts, immune deconvolution) |
| `pipeline/01_tf_scoring/` | Fig. 1A, 1E; Table S1 | Evidence scoring of 1,639 TFs, selection of 13 CRC TFs, and TF–context Kruskal–Wallis tests |
| `pipeline/02_chromatin/` | Fig. 3D; Table S6 | Local ATAC accessibility vs. RNA expression coupling (±50 kb) |
| `pipeline/03_bulk_modules/` | Fig. 6 | TCGA-BRCA CRC module scores, tumor-state signatures, and MCP-counter correlations |
| `pipeline/04_single_cell/` | Fig. 7 | GSE176078 malignant epithelial cell clustering and six CRC transcriptional states |
| `pipeline/05_outcome_dependency/` | Fig. 8; Table S7 | State scores, continuous Cox models (OS/PFS), METABRIC/GSE96058 projection, and DepMap CRISPR dependencies |
| `revision_analyses/` | Figs. 2, 4C, S1–S4; Tables S2, S4, S5, S9, S10, S14 | First-round (`R1`) and second-round (`A0–A22`) peer-review audits, controls, and supplementary ChIP analyses |
| `data/` & `source_tables/` | Reference inputs | Curated dataset manifests, submitted network reference, context membership table, and original panel source tables |

---

## Requirements & Setup

The analyses reported in the paper were audited under **Python 3.12–3.13**, **R 4.5.2**, and **Node.js** (see [`software_versions.tsv`](software_versions.tsv) for exact versions).

### 1. Python (≥ 3.10)

```bash
git clone https://github.com/Yongjia-Z/BRCA-CRC-Scripts.git
cd BRCA-CRC-Scripts

# Create an isolated environment (recommended)
conda create -n brca-crc python=3.12 -y
conda activate brca-crc

# Install core Python dependencies
pip install -r requirements.txt

# Optional: faster genomic interval overlaps (vectorized fallback is used if absent)
# pip install pyranges
```

### 2. R (for Cox models and immune deconvolution)

Install `survival` (used by `pipeline/05_outcome_dependency/tcga_state_cox.R`) and `MCPcounter` (used by `pipeline/00_data/immune_deconv/00_mcp_counter.R`):

```r
install.packages("survival")
# Install MCPcounter from GitHub if not already installed:
# remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")
```

### 3. Node.js (only for A16 bigWig quantification)

```bash
npm install @gmod/bbi
```

---

## Data Preparation & Environment Variables

Obtain public datasets from the accessions listed in **Tables S2 and S5e** of the paper (with curated dataset lists provided in [`data/dataset_lists/`](data/dataset_lists/)) and organize them under `datasets/<resource>/`.

### Environment Variables

| Variable | Default | Purpose |
|---|---|---|
| `CRC_DATA_ROOT` | `<repo>/datasets` | Root directory for downloaded public datasets read by pipeline scripts |
| `CRC_CORE_ROOT` | `<repo>/core_data` | Output directory for processed intermediate tables |
| `CRC_PROJECT_ROOT` | `<repo>` | Repository root used by `revision_analyses/A16_supplementary_ChIP/scripts/` |

> [!NOTE]
> Most scripts in `revision_analyses/` take the repository root (`<project_root>`) as their first positional CLI argument.

### Manually Placed Files

In addition to standard resources under `datasets/<resource>/`, place the following files at their expected locations before running the corresponding analyses:

| Resource | Expected Path | Consumed By |
|---|---|---|
| GSE251868 HER2 H3K27ac BED files | `datasets/GSE251868_HER2_H3K27ac/data` | `pipeline/01_tf_scoring/tf_evidence_scoring.py` |
| GSE245830 TNBC H3K27ac (`GSE245830_RAW.tar`) | `datasets/GSE245830_TNBC_H3K27ac/data` (pass extracted folder as 2nd arg) | `revision_analyses/A7_Basal_TNBC_H3K27ac_layer/` |
| GSE213840 FOXC1 TNBC peaks | `datasets/GSE213840_FOXC1_TNBC_ChIPseq/` (folders matching `GSE*_<TF>_TNBC*`) | `pipeline/00_data/remap/03_tnbc_direct_binding_overlap.py` |
| GSE213840 FOXC1 ChIP & input bigWig (`GSM6595097–GSM6595101`) | `datasets/Figure4_browser_tracks/data` | `revision_analyses/A16_supplementary_ChIP/scripts/a16_signal_quantify.mjs` |
| Aran et al. (2015) Supplementary Data 1 | `revision_analyses/A1_purity_adjustment/Aran2015_SupData1.xlsx` | `revision_analyses/A1_purity_adjustment/` |
| ChIP-Atlas `experimentList.tab` | `revision_analyses/A16_supplementary_ChIP/raw/metadata/` | `A15` and `A18` |
| Curated ChIP-Atlas PGR peak files | `revision_analyses/A14_occupancy_and_network_checks/raw_chipatlas_PGR/<SRX>.05.bed` | `A14` |
| A16 peak files, GEO SOFT records & liftOver chains | `revision_analyses/A16_supplementary_ChIP/raw/` and `derived/` (exact paths in `data/dataset_lists/A16_peak_file_manifest.tsv`) | `A16` |

---

## Execution Order

Run the stages in the sequence below (within each folder, execute scripts in numerical prefix order):

1. **Stage 0 — Public data preprocessing (`pipeline/00_data/`):**
   1. `tf_resources`, `gencode`
   2. `tcga` (`03` → `08`)
   3. `cistrome_cancer/03`, then `tcga/09`
   4. `cistrome_cancer/05` and `encode/01`, then `tcga` (`10` → `11`)
   5. `encode/02` and `sedb/01`, then `tcga/12` and `trrust/00`
   6. `tcga/14` and `tcga/17`
   7. `gse85158_h3k27ac`, `scatac_multiome`, `remap`
   8. `gse96058`, `gse176078`, `immune_deconv`

2. **Stage 1 — TF evidence scoring & subtype association (`pipeline/01_tf_scoring/`):**
   ```bash
   python pipeline/01_tf_scoring/tf_evidence_scoring.py \
       --raw-root datasets \
       --remap-result-dir core_data/remap_reanalysis/results \
       --output-dir source_tables/figure1 \
       --intermediate-dir source_tables/intermediate/figure1a

   python pipeline/01_tf_scoring/prepare_subtype_association_inputs.py

   python pipeline/01_tf_scoring/subtype_association.py \
       --input-dir source_tables/figure1/subtype_association_inputs \
       --output-dir source_tables/figure1
   ```

3. **Stages 2 & 3 — Chromatin coupling & bulk modules:**
   `pipeline/02_chromatin/atac_rna_coupling.py` → `pipeline/03_bulk_modules/module_scores_tumor_state_mcp.py`

4. **Stage 4 — Single-cell states (`pipeline/04_single_cell/`):**
   `01_process_and_score_cells.py` → `02_state_discovery.py` → `03_state_scores_per_cell.py`

5. **Stage 5 — Clinical outcome & functional dependency (`pipeline/05_outcome_dependency/`):**
   `tcga_state_scores.py` → `external_cohort_state_genes.py` → `depmap_state_scores.py` → `tcga_state_cox.py` (invokes `tcga_state_cox.R`)

6. **Stage 6 — Peer-review analyses (`revision_analyses/`):**
   1. `R1_first_round` (`R1_analyze_revision.py`, `R1_external_projection.py`)
   2. `A7` (Basal/TNBC H3K27ac layer update)
   3. `A8`, `A15`, `A18`, `A21`
   4. `A14` (`02` → `06`, `09`)
   5. `A16/scripts/` in order: `a16_liftover.py` → `a16_bootstrap.py` → `a16_analyze_peaks.py` → `a16_signal_quantify.mjs` → `a16_quality.py` → `a16_sensitivity.py`
   6. `A19` (`A19_occupancy_network.py`, `A19_02_self_locus_100kb.py`), then `A17`, `A22`, `A13`
   7. `A10`, then `A1` (passing `A10`'s `A10_fig6_scores_corrected_membership.tsv` as the 3rd argument), then `A0`, `A12`

---

## Script Reference

Scripts without CLI arguments read paths from `CRC_DATA_ROOT` and `CRC_CORE_ROOT`. Expand each module below for the full script catalog (70 scripts total).

<details>
<summary><strong>Helper Module (<code>_repro_common/</code> — 1 script)</strong></summary>

| Script | Purpose | Arguments |
|---|---|---|
| `crc_repro.py` | Shared helpers (table I/O, PAM50-to-context mapping, BH FDR correction, gene-by-sample matrices) used by `subtype_association.py` and `A12`. | Imported by other scripts |

</details>

<details>
<summary><strong>Stage 0: Public Data Preprocessing (<code>pipeline/00_data/</code> — 27 scripts)</strong></summary>

| Script | Purpose | Arguments |
|---|---|---|
| `cistrome_cancer/03_candidate_tf_tf_edges.py` | Builds candidate TF-to-TF edges from Cistrome Cancer CancerTarget target tables. | – |
| `cistrome_cancer/05_enhancer_signal_near_tf_genes.py` | Summarizes TCGA enhancer signal within candidate TF gene windows. | – |
| `encode/01_mcf7_h3k27ac_overlap.py` | Counts overlaps of ENCODE MCF-7 H3K27ac peaks (`data/dataset_lists/`) with candidate TF gene windows. | – |
| `encode/02_mcf7_rose_like_super_enhancers.py` | Stitches MCF-7 H3K27ac peaks in a ROSE-like manner as a super-enhancer proxy and links them to candidate TF loci. | – |
| `gencode/00_parse_gencode_v36.py` | Parses GENCODE v36 annotation into a gene reference table. | – |
| `gse176078/00_convert_mtx_to_h5ad.py` | Converts the GSE176078 count matrix and metadata into an `.h5ad` file. | – |
| `gse85158_h3k27ac/02_quantify_hg19_to_hg38.py` | Quantifies GSE85158 H3K27ac signal in TF gene windows (GENCODE v19), lifts windows to hg38, and summarizes by subtype. | `[--limit-files N] [--max-workers N]` |
| `gse96058/00_annotation_and_module_scores.py` | Reads GSE96058 sample annotation and expression and computes CRC module scores (used by `A18` and `R1`). | – |
| `immune_deconv/00_mcp_counter.R` | Computes MCP-counter immune and stromal estimates for TCGA-BRCA (input for Fig. 6C, D; requires `MCPcounter`). | – |
| `remap/02_breast_promoter_overlap.py` | Counts peak overlaps of each ReMap breast tier-1 dataset (`data/dataset_lists/`) with candidate TF promoters. | – |
| `remap/03_tnbc_direct_binding_overlap.py` | Overlaps TNBC FOXC1 and TP63 ChIP-seq/CUT&Tag peaks with candidate TF promoters and ±100 kb loci (TNBC part of Fig. 4B, C). | – |
| `scatac_multiome/00_tcga_brca_scatac_locus_support.py` | Counts accessible peaks and motifs at CRC TF loci in breast cancer and healthy breast scATAC peak sets (Fig. 3A–C; scATAC layer). | – |
| `sedb/01_dbsuper_breast_super_enhancers.py` | Lifts dbSUPER breast super-enhancers to hg38 and counts overlaps within ±100 kb of TF loci. | – |
| `tcga/03_atac_peak_summary.py` | Summarizes sample- and peak-level accessibility of the TCGA ATAC peak matrix. | – |
| `tcga/04_atac_footprinting.py` | Reads TCGA ATAC footprinting and motif tables and correlates footprint depth with TF expression (footprint layer). | – |
| `tcga/06_tf_expression_activity_screen.py` | Builds the TF expression matrix and sample annotation, tests expression/subtype association, and computes DoRothEA regulon TF activity (RNA & activity layers; matrices of Fig. 1E, F). | – |
| `tcga/07_atac_peak_annotation_subtype.py` | Annotates ATAC peaks and tests for subtype-specific accessibility peaks. | – |
| `tcga/08_atac_candidate_tf_integration.py` | Integrates subtype-specific peaks, peak-to-gene links, and motif enrichment into the candidate TF evidence table. | – |
| `tcga/09_cistrome_tf_tf_edges.py` | Adds Cistrome Cancer TF–TF edges to the candidate TF evidence table. | – |
| `tcga/10_enhancer_h3k27ac_evidence.py` | Adds TCGA enhancer and ENCODE MCF-7 H3K27ac evidence. | – |
| `tcga/11_crc_member_evidence_scores.py` | Computes CRC member evidence scores of candidate TFs and lists high-confidence candidates. | – |
| `tcga/12_tf_tf_edge_evidence_scores.py` | Scores the evidence of every candidate TF–TF edge (edge evidence table of the Fig. 2 network). | – |
| `tcga/14_crc_activity_scores.py` | Computes per-sample expression and activity module scores of the CRC TFs. | – |
| `tcga/17_clinical_outcome_dataset.py` | Merges CRC module scores with TCGA-BRCA clinical follow-up and runs clinical association and survival (log-rank, Cox) analyses. | – |
| `tf_resources/00_lambert_tf_list.py` | Extracts the 1,639 human TFs from the Lambert et al. (2018) table. | – |
| `tf_resources/01_convert_dorothea_rda.R` | Converts the DoRothEA regulon `.rda` file to a `.tsv` table. | `[input.rda] [output.tsv]` |
| `trrust/00_trrust_edges.py` | Adds TRRUST v2 human TF–target interactions to the TF–TF edge evidence table. | – |

</details>

<details>
<summary><strong>Stages 1–5: Core Manuscript Pipeline (<code>pipeline/01_tf_scoring/</code> – <code>05_outcome_dependency/</code> — 13 scripts)</strong></summary>

| Stage / Script | Purpose | Arguments |
|---|---|---|
| `01_tf_scoring/tf_evidence_scoring.py` | **Fig. 1A, Table S1:** Scores 1,639 TFs across 7 core and 3 supportive layers, selects 13 CRC TFs by tier rules, and assigns context membership. | `--raw-root --remap-result-dir --output-dir --intermediate-dir [--max-workers ...]` |
| `01_tf_scoring/prepare_subtype_association_inputs.py` | Exports the three inputs required by `subtype_association.py` (TF expression matrix, sample annotation, CRC TF table). | `[out_dir]` |
| `01_tf_scoring/subtype_association.py` | **Fig. 1E:** Kruskal–Wallis test of 1,631 TFs across Luminal, HER2, Basal/TNBC, and Normal-like groups with BH correction; outputs all TFs and top 20. | `--input-dir DIR --output-dir DIR [--top-n 20]` |
| `02_chromatin/atac_rna_coupling.py` | **Fig. 3D, Table S6:** Correlation between ATAC accessibility within ±50 kb of each TF locus and TF expression. | – |
| `03_bulk_modules/module_scores_tumor_state_mcp.py` | **Fig. 6:** TCGA-BRCA CRC module scores (expression & activity) and correlations with tumor-state signatures and MCP-counter estimates before/after context adjustment. | – |
| `04_single_cell/01_process_and_score_cells.py` | **Fig. 7A:** Processes GSE176078, computes UMAP embedding, and scores per-cell CRC TF modules. | – |
| `04_single_cell/02_state_discovery.py` | **Fig. 7B, C:** Clusters TF co-expression within ER+, HER2+, and TNBC malignant cells to define the six malignant states and signature genes. | – |
| `04_single_cell/03_state_scores_per_cell.py` | **Fig. 7C:** Computes per-cell state scores and dominant state assignments. | – |
| `05_outcome_dependency/tcga_state_scores.py` | **Fig. 8A:** Computes the six state scores in TCGA-BRCA and Kruskal–Wallis tests across contexts. | – |
| `05_outcome_dependency/external_cohort_state_genes.py` | Projects the six state signatures onto TCGA, METABRIC, and GSE96058 and fits continuous Cox models. | – |
| `05_outcome_dependency/depmap_state_scores.py` | **Fig. 8C:** Computes state expression scores in DepMap breast cancer cell lines and correlates them with CRISPR gene effect. | – |
| `05_outcome_dependency/tcga_state_cox.py` | **Fig. 8B, Table S7:** Continuous Cox models of each state score in TCGA-BRCA (OS, PFS; per 1 SD; 120-month truncation) via `tcga_state_cox.R`. | `[out_dir]` |
| `05_outcome_dependency/tcga_state_cox.R` | Invoked by `tcga_state_cox.py`: fits Cox models and proportional-hazards tests using `survival`. | `<model_input.tsv> <results.tsv>` |

</details>

<details>
<summary><strong>Stage 6: Peer-Review Revision Analyses (<code>revision_analyses/</code> — 29 scripts)</strong></summary>

| Script | Purpose | Arguments |
|---|---|---|
| `R1_first_round/R1_analyze_revision.py` | First-round review: XBP1/ZNF703 selection evidence, layer omission, weight perturbation & random sets, Fig. 3D peak audit, Fig. 6A/8A medians, self-locus occupancy, BH-corrected/clinically adjusted Cox models, and DepMap matched null/permutations (tables `R1_1–R1_6`). | `[project_root] [out_dir]` |
| `R1_first_round/R1_external_projection.py` | First-round review: projection of the six state signatures onto METABRIC and GSE96058 (`R1_7`; Fig. S4, Table S10). | `[project_root] [out_dir]` |
| `A0_reanalysis_checks/A0_verify.py` | Independent re-computation verifying selection rules, random sets, Fig. 3D statistics, Cox models, DepMap matched null, external projection, and Fig. 6A/8A medians. | `<project_root> <out_dir>` |
| `A1_purity_adjustment/A1_purity_adjusted_immune.py` | **Fig. S3, Tables S9a–b:** Purity-adjusted (CPE) correlations of module scores with MCP-counter estimates. | `<project_root> <out_dir> [scores_tsv]` |
| `A7_Basal_TNBC_H3K27ac_layer/A7_basal_h3k27ac_audit.py` | Incorporates GSE245830 TNBC H3K27ac peaks and recomputes the Basal/TNBC layer and context membership (Tables S1, S4d, S4e). | `<project_root> <extracted_GSE245830_dir> <out_dir>` |
| `A8_selection_robustness/A8_selection_robustness_corrected.py` | **Fig. S1, Tables S4a–c:** Selection robustness on the revised score table (weight perturbation, layer omission, threshold shifts, random sets). | `<A7_corrected_csv> <out_dir>` |
| `A10_Fig6_modules/A10_fig6_modules_corrected.py` | **Tables S9c–d:** Fig. 6 module scores using the revised Luminal membership (excluding RUNX1). | `<project_root> <A7_corrected_csv> <out_dir>` |
| `A12_Fig3D_mean_accessibility_script/figure3d_corrected.py` | Corrected Fig. 3D statistics: mean local accessibility vs. expression and BH correction across all local peak tests. | `--input-dir DIR --output-dir DIR [--flank-bp 50000]` |
| `A13_motif_context_membership/A13_motif_membership.py` | **Table S14:** Classifies Fig. 2C network motifs by node context membership. | `<figure2c.csv> <A7_CRC_TF_calls_stored_vs_corrected.tsv> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_02_pgr_target_and_network.py` | Adds PGR as an occupancy target, rebuilds edges with ESR1→PGR, and enumerates motifs. | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_03_chipatlas_crosscheck.py` | Cross-checks promoter occupancy of the 13 TFs against local ChIP-Atlas breast peak files. | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_04_knockdown_label_sensitivity.py` | Evaluates the impact of excluding five tier-1 datasets with knockdown labels missed by metadata filtering. | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_05_tables_and_control.py` | Generates updated occupancy matrix (ReMap PGR column, ChIP-Atlas PGR row), self-occupancy, and network-rule tables. | `<project_root> <S5_xlsx> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_06_edges_table.py` | Formats the table of network edges scoring ≥ 6 (Table S14b). | `<out_dir>` |
| `A14_occupancy_and_network_checks/A14_09_pgr_row_binding_knockdown.py` | Final occupancy tables from curated ReMap and ChIP-Atlas PGR datasets (`data/dataset_lists/`): PGR column/row/self-locus values and binding-score effect on selection. | `<project_root> <Table_S5_workbook> <out_dir>` |
| `A15_selection_coverage_audit/A15_selection_bias_audit.py` | Read-only audit of layer coverage, sparse-layer membership, rankings without sparse layers, uniform footprint evidence, alternative-scoring rule replay, and TFs with breast ChIP but zero binding score. | `<project_root> <revision_analyses_dir> <out_dir> <experimentList.tab>` |
| `A16_supplementary_ChIP/scripts/a16_common.py` | Shared A16 paths, file I/O, checksums, promoter table, and peak-overlap helpers. | Imported by A16 scripts |
| `A16_supplementary_ChIP/scripts/a16_liftover.py` | *(Order 1)* Chain-based interval liftOver (≥ 95% bases in one chain), cross-checked with `pyliftover`. | – |
| `A16_supplementary_ChIP/scripts/a16_bootstrap.py` | *(Order 2)* Reproduces the submitted manuscript PGR baseline and writes promoter windows. | – |
| `A16_supplementary_ChIP/scripts/a16_analyze_peaks.py` | *(Order 3)* Reproduces the earlier occupancy matrix and analyzes supplementary ChIP peaks by source (`data/dataset_lists/`). | – |
| `A16_supplementary_ChIP/scripts/a16_signal_quantify.mjs` | *(Order 4)* Quantifies FOXC1 genome-wide bigWig signal at TF promoters without peak calling (Node.js + `@gmod/bbi`). | `node a16_signal_quantify.mjs <project_root>` |
| `A16_supplementary_ChIP/scripts/a16_quality.py` | *(Order 5)* Interval audit of peak files and cross-pipeline correspondence. | – |
| `A16_supplementary_ChIP/scripts/a16_sensitivity.py` | *(Order 6)* **Table S4g:** Conditional replay on original 13-TF inputs modifying only the binding layer. | – |
| `A17_score5_motifs/A17_01_score5_motifs.py` | **Tables S14d–e:** Tests whether HER2 or Basal/TNBC sets form closed loops when edges scoring 5 are admitted. | `<A16_results_dir> <out_dir>` |
| `A18_selection_and_subtype_checks/A18_selection_subtype_checks.py` | **Tables S2c–d, S4h, S4k–l:** Reviewer checks (R1.5, R5.2–R5.5): Fig. 1E ranks, 5-layer reference, balanced subtype composition, PAM50 vs. receptor groups, and resource reuse. | `<project_root> <revision_analyses_dir> <A16_dir> <out_dir>` |
| `A19_occupancy_network_self_locus/A19_occupancy_network.py` | **Figs. 2, 4C; Tables S5, S14:** Occupancy matrix, TF network/motifs with supplementary ChIP, and replication in studies excluded from TF selection. | `<project_root> <revision_analyses_dir> <A16_dir> <out_dir>` |
| `A19_occupancy_network_self_locus/A19_02_self_locus_100kb.py` | **Fig. S2B, Table S5h:** Self-locus occupancy within TSS ±100 kb ranked against all protein-coding genes. | `<project_root> <A16_dir> <out_dir>` |
| `A21_RNA_grouping_subsampling/A21_rna_grouping_balanced.py` | **Tables S4i, S4i2:** Grouping definition of the RNA-layer subtype test and recalculation in balanced/proportional subsamples. | `<project_root> <revision_analyses_dir> <out_dir>` |
| `A22_binding_by_cell_line/A22_cell_line_binding.py` | **Tables S14f–g:** Cell-line provenance of each network edge and single-cell-line closed structures. | `<A16_dir> <A19_out_dir> <out_dir>` |

</details>

---

## Citation

If you use this code, please cite the accompanying manuscript (see also [`CITATION.cff`](CITATION.cff)):

```bibtex
@article{Zeng2026BRCACRC,
  author  = {Zeng, Yongjia and Shi, Songfu and Wang, Fan and Zhao, Lixin and Ye, Lingxi and Wu, Xiaodi and Yan, Dong},
  title   = {Subtype-specific Core Transcriptional Regulatory Circuitry Links Transcriptional States, Chromatin Context, and Functional Dependencies in Breast Cancer},
  journal = {Functional \& Integrative Genomics},
  year    = {2026},
  note    = {Under review}
}
```

## License

Released under the [MIT License](LICENSE). Third-party public datasets and software tools remain subject to their respective licenses and data-use terms and should be cited as in the paper.

---
---

<a id="chinese"></a>

## 项目简介（简体中文）

本仓库为以下研究论文的配套生物信息学分析代码：

> Zeng Y, Shi S, Wang F, Zhao L, Ye L, Wu X, Yan D. **Subtype-specific Core Transcriptional Regulatory Circuitry Links Transcriptional States, Chromatin Context, and Functional Dependencies in Breast Cancer.** *Functional & Integrative Genomics*（审稿中）。

本项目整合多组学乳腺癌公共资源——包括 bulk RNA-seq、DoRothEA 调控子（regulon）活性、ATAC-seq 染色质开放性/基序/足迹（footprint）、亚型匹配的 H3K27ac 增强子代理信号、乳腺背景下的直接 TF 结合证据（ReMap、ChIP-Atlas、GEO）、scATAC-seq、scRNA-seq、外部队列（METABRIC、GSE96058）以及 DepMap CRISPR 基因依赖性数据——对 1,639 个人类转录因子（TF）进行多维度证据评分，筛选出覆盖 Luminal、HER2 与 Basal/TNBC 背景的 13 个乳腺癌核心转录调控环路（CRC）TF，并在单细胞层面解析 6 种恶性上皮细胞转录状态及其临床预后与功能依赖特征。本仓库同时完整收录了同行评审阶段新增的审稿补充分析（`R1` 及 `A0–A22`）。

### 仓库范围说明

- **包含内容：** 公共数据预处理、多维 TF 证据评分、局部染色质-表达偶联、bulk 模块评分与免疫浸润相关性、单细胞恶性状态识别、生存与 DepMap 依赖性分析、审稿敏感性与审计脚本、经人工核查的数据集清单（`data/dataset_lists/`）、公共数据无法重新生成的参考输入（`data/`）以及原始面板数据表（`source_tables/`）。
- **不包含内容：** 自动下载脚本、原始或受控访问的大型公共数据集、最终发表级配图绘制与矢量排版代码，以及格式化补充表格（最终结果表已作为论文 **Online Resource 1** 提供）。

---

## 目录结构与论文对应关系

```text
BRCA-CRC-Scripts/
├── pipeline/
│   ├── 00_data/                 # 公共数据预处理与 TF 评分各证据层（27 个脚本）
│   ├── 01_tf_scoring/           # 1,639 TF 证据评分、13 CRC TF 筛选及亚型关联检验（3 个脚本）
│   ├── 02_chromatin/            # TF 位点 ±50 kb 局部 ATAC–RNA 偶联分析（1 个脚本）
│   ├── 03_bulk_modules/         # TCGA-BRCA 模块评分、肿瘤状态及 MCP-counter 相关性（1 个脚本）
│   ├── 04_single_cell/          # GSE176078 预处理与 6 种恶性上皮细胞状态识别（3 个脚本）
│   └── 05_outcome_dependency/   # TCGA/外部队列状态评分、Cox 模型与 DepMap 依赖分析（5 个脚本）
├── revision_analyses/           # 同行评审阶段补充分析、鲁棒性验证与补充 ChIP 分析（29 个脚本）
├── _repro_common/               # 公共辅助函数：表读写、PAM50 背景映射、BH 校正（1 个脚本）
├── data/
│   ├── crc_tf_table_modules.tsv # 定义 Figs. 6–8 模块归属的上下文成员表
│   ├── dataset_lists/           # 筛选后的 ReMap、ENCODE、ChIP-Atlas 与 GEO 清单（Tables S5c–S5e）
│   └── submitted_network/       # 初稿网络边集合与首轮启动子结合表
├── source_tables/               # 审稿分析所读取的初版分析面板结果表
├── requirements.txt             # Python 依赖包列表
├── software_versions.tsv        # 论文实际核查使用的 Python、R 与 Node.js 软件版本
├── CITATION.cff                 # 机器可读引用元数据
└── LICENSE                      # MIT 开源许可证
```

| 目录路径 | 对应论文图表 | 主要内容 |
|---|---|---|
| `pipeline/00_data/` | 上游证据层（Figs. 1–4, 6） | 公共数据预处理（TCGA、ATAC、H3K27ac、ReMap、scATAC、外部队列与免疫浸润反卷积） |
| `pipeline/01_tf_scoring/` | Fig. 1A, 1E; Table S1 | 1,639 个 TF 多维证据评分、13 个 CRC TF 筛选及亚型 Kruskal–Wallis 检验 |
| `pipeline/02_chromatin/` | Fig. 3D; Table S6 | TF 位点 ±50 kb 局部 ATAC 开放性与表达量的偶联分析 |
| `pipeline/03_bulk_modules/` | Fig. 6 | TCGA-BRCA CRC 模块评分、肿瘤状态特征及 MCP-counter 相关性（含背景校正前后） |
| `pipeline/04_single_cell/` | Fig. 7 | GSE176078 恶性上皮细胞聚类与 6 种 CRC 转录状态发现 |
| `pipeline/05_outcome_dependency/` | Fig. 8; Table S7 | 状态评分、连续 Cox 模型（OS/PFS）、METABRIC/GSE96058 外部映射及 DepMap CRISPR 依赖性 |
| `revision_analyses/` | Figs. 2, 4C, S1–S4; Tables S2, S4, S5, S9, S10, S14 | 第一轮（`R1`）及后续审稿补充分析（`A0–A22`）：纯度校正、TNBC H3K27ac 更新、稳健性检验、补充 ChIP 与细胞系溯源等 |
| `data/` 与 `source_tables/` | 参考输入与基础表 | 数据集清单、初稿网络参考表、模块成员表及初版面板数据表 |

---

## 运行环境与依赖安装

论文分析环境经核查使用 **Python 3.12–3.13**、**R 4.5.2** 及 **Node.js**，详细版本清单见 [`software_versions.tsv`](software_versions.tsv)。

### 1. Python（≥ 3.10）

```bash
git clone https://github.com/Yongjia-Z/BRCA-CRC-Scripts.git
cd BRCA-CRC-Scripts

# 建议创建独立虚拟环境
conda create -n brca-crc python=3.12 -y
conda activate brca-crc

# 安装核心 Python 依赖
pip install -r requirements.txt

# 可选：安装 pyranges 以加速基因组区间重叠计算（缺失时自动回退至向量化实现）
# pip install pyranges
```

### 2. R（用于 Cox 生存分析与 MCP-counter 免疫反卷积）

需安装 `survival`（供 `pipeline/05_outcome_dependency/tcga_state_cox.R` 调用）与 `MCPcounter`（供 `pipeline/00_data/immune_deconv/00_mcp_counter.R` 调用）：

```r
install.packages("survival")
# 如尚未安装 MCPcounter，可通过 GitHub 安装：
# remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")
```

### 3. Node.js（仅 A16 的 bigWig 信号定量步骤需要）

```bash
npm install @gmod/bbi
```

---

## 数据准备与环境变量

请依据论文 **Table S2 与 Table S5e** 所列数据来源自行下载公共数据集（精选数据集清单见 [`data/dataset_lists/`](data/dataset_lists/)），并按 `datasets/<resource>/` 目录结构存放。

### 环境变量说明

| 环境变量 | 默认值 | 用途 |
|---|---|---|
| `CRC_DATA_ROOT` | `<repo>/datasets` | 公共数据集根目录，供主流程脚本读取原始数据 |
| `CRC_CORE_ROOT` | `<repo>/core_data` | 中间处理结果表输出目录 |
| `CRC_PROJECT_ROOT` | `<repo>` | 仓库根目录，供 `revision_analyses/A16_supplementary_ChIP/scripts/` 读取 |

> [!NOTE]
> `revision_analyses/` 下的多数脚本默认将仓库根目录（`<project_root>`）作为首个位置命令行参数传入。

### 需手动放置的数据文件

除常规 `datasets/<resource>/` 外，运行对应模块前需将以下文件放置至指定路径：

| 数据资源 | 预期存放路径 | 读取模块 |
|---|---|---|
| GSE251868 HER2 H3K27ac BED 文件 | `datasets/GSE251868_HER2_H3K27ac/data` | `pipeline/01_tf_scoring/tf_evidence_scoring.py` |
| GSE245830 TNBC H3K27ac（`GSE245830_RAW.tar`） | `datasets/GSE245830_TNBC_H3K27ac/data`（将解压目录作为 A7 第二个参数） | `revision_analyses/A7_Basal_TNBC_H3K27ac_layer/` |
| GSE213840 FOXC1 TNBC 结合峰 | `datasets/GSE213840_FOXC1_TNBC_ChIPseq/`（匹配 `GSE*_<TF>_TNBC*` 目录名） | `pipeline/00_data/remap/03_tnbc_direct_binding_overlap.py` |
| GSE213840 FOXC1 ChIP 与 input bigWig（`GSM6595097–GSM6595101`） | `datasets/Figure4_browser_tracks/data` | `revision_analyses/A16_supplementary_ChIP/scripts/a16_signal_quantify.mjs` |
| Aran et al. (2015) Supplementary Data 1 | `revision_analyses/A1_purity_adjustment/Aran2015_SupData1.xlsx` | `revision_analyses/A1_purity_adjustment/` |
| ChIP-Atlas `experimentList.tab` | `revision_analyses/A16_supplementary_ChIP/raw/metadata/` | `A15` 与 `A18` |
| 筛选后的 ChIP-Atlas PGR 结合峰文件 | `revision_analyses/A14_occupancy_and_network_checks/raw_chipatlas_PGR/<SRX>.05.bed` | `A14` |
| A16 结合峰文件、GEO SOFT 记录及 liftOver chain 文件 | `revision_analyses/A16_supplementary_ChIP/raw/` 与 `derived/`（具体路径见 `data/dataset_lists/A16_peak_file_manifest.tsv`） | `A16` |

---

## 运行顺序

请按以下顺序依次执行（同一文件夹内按脚本数字编号顺序运行）：

1. **阶段 0 — 公共数据预处理（`pipeline/00_data/`）：**
   1. `tf_resources`、`gencode`
   2. `tcga`（`03` → `08`）
   3. `cistrome_cancer/03`，随后 `tcga/09`
   4. `cistrome_cancer/05` 与 `encode/01`，随后 `tcga`（`10` → `11`）
   5. `encode/02` 与 `sedb/01`，随后 `tcga/12` 与 `trrust/00`
   6. `tcga/14` 与 `tcga/17`
   7. `gse85158_h3k27ac`、`scatac_multiome`、`remap`
   8. `gse96058`、`gse176078`、`immune_deconv`

2. **阶段 1 — TF 多维证据评分与 Fig. 1E 亚型关联（`pipeline/01_tf_scoring/`）：**
   ```bash
   python pipeline/01_tf_scoring/tf_evidence_scoring.py \
       --raw-root datasets \
       --remap-result-dir core_data/remap_reanalysis/results \
       --output-dir source_tables/figure1 \
       --intermediate-dir source_tables/intermediate/figure1a

   python pipeline/01_tf_scoring/prepare_subtype_association_inputs.py

   python pipeline/01_tf_scoring/subtype_association.py \
       --input-dir source_tables/figure1/subtype_association_inputs \
       --output-dir source_tables/figure1
   ```

3. **阶段 2 与 3 — 染色质偶联与 bulk 模块分析：**
   `pipeline/02_chromatin/atac_rna_coupling.py` → `pipeline/03_bulk_modules/module_scores_tumor_state_mcp.py`

4. **阶段 4 — 单细胞状态分析（`pipeline/04_single_cell/`）：**
   `01_process_and_score_cells.py` → `02_state_discovery.py` → `03_state_scores_per_cell.py`

5. **阶段 5 — 临床预后与功能依赖分析（`pipeline/05_outcome_dependency/`）：**
   `tcga_state_scores.py` → `external_cohort_state_genes.py` → `depmap_state_scores.py` → `tcga_state_cox.py`（自动调用 `tcga_state_cox.R`）

6. **阶段 6 — 同行评审补充分析（`revision_analyses/`）：**
   1. `R1_first_round`（`R1_analyze_revision.py`、`R1_external_projection.py`）
   2. `A7`（更新 Basal/TNBC H3K27ac 证据层）
   3. `A8`、`A15`、`A18`、`A21`
   4. `A14`（`02` → `06`、`09`）
   5. `A16/scripts/` 按顺序：`a16_liftover.py` → `a16_bootstrap.py` → `a16_analyze_peaks.py` → `a16_signal_quantify.mjs` → `a16_quality.py` → `a16_sensitivity.py`
   6. `A19`（`A19_occupancy_network.py`、`A19_02_self_locus_100kb.py`），随后 `A17`、`A22`、`A13`
   7. `A10`，随后 `A1`（将 `A10` 生成的 `A10_fig6_scores_corrected_membership.tsv` 作为第三个参数传入），最后运行 `A0`、`A12`

---

## 完整脚本清单

未列出参数的脚本默认从 `CRC_DATA_ROOT` 与 `CRC_CORE_ROOT` 读取路径。点击下方折叠面板可查看全部 70 个脚本的功能与命令行参数说明。

<details>
<summary><strong>公共辅助模块（<code>_repro_common/</code> — 1 个脚本）</strong></summary>

| 脚本 | 功能说明 | 参数 |
|---|---|---|
| `crc_repro.py` | 公共工具函数（表格读写、PAM50 至分析背景映射、BH 校正、基因×样本矩阵构建），供 `subtype_association.py` 与 `A12` 导入。 | 被其他脚本导入 |

</details>

<details>
<summary><strong>阶段 0：公共数据预处理（<code>pipeline/00_data/</code> — 27 个脚本）</strong></summary>

| 脚本 | 功能说明 | 参数 |
|---|---|---|
| `cistrome_cancer/03_candidate_tf_tf_edges.py` | 基于 Cistrome Cancer CancerTarget 靶基因表构建候选 TF–TF 调控边。 | – |
| `cistrome_cancer/05_enhancer_signal_near_tf_genes.py` | 汇总候选 TF 基因窗口内的 TCGA 增强子信号。 | – |
| `encode/01_mcf7_h3k27ac_overlap.py` | 统计 ENCODE MCF-7 H3K27ac 峰（文件清单见 `data/dataset_lists/`）与候选 TF 基因窗口的重叠数。 | – |
| `encode/02_mcf7_rose_like_super_enhancers.py` | 采用类 ROSE 峰拼接策略构建 MCF-7 H3K27ac 超级增强子代理信号，并关联至候选 TF 位点。 | – |
| `gencode/00_parse_gencode_v36.py` | 解析 GENCODE v36 注释文件生成基因参考表。 | – |
| `gse176078/00_convert_mtx_to_h5ad.py` | 将 GSE176078 表达矩阵与元数据转换为 `.h5ad` 文件。 | – |
| `gse85158_h3k27ac/02_quantify_hg19_to_hg38.py` | 在 GENCODE v19 基因窗口中定量 GSE85158 H3K27ac 信号，转换坐标至 hg38 并按亚型汇总。 | `[--limit-files N] [--max-workers N]` |
| `gse96058/00_annotation_and_module_scores.py` | 读取 GSE96058 样本注释与表达矩阵并计算 CRC 模块评分（供 `A18` 与 `R1` 使用）。 | – |
| `immune_deconv/00_mcp_counter.R` | 计算 TCGA-BRCA 的 MCP-counter 免疫与基质细胞评分（Fig. 6C, D 输入；依赖 `MCPcounter`）。 | – |
| `remap/02_breast_promoter_overlap.py` | 统计各 ReMap 乳腺 tier-1 数据集（清单见 `data/dataset_lists/`）与候选 TF 启动子的峰重叠。 | – |
| `remap/03_tnbc_direct_binding_overlap.py` | 将 TNBC FOXC1 和 TP63 的 ChIP-seq/CUT&Tag 峰与候选 TF 启动子及 ±100 kb 位点取交集（Fig. 4B, C 中 TNBC 部分）。 | – |
| `scatac_multiome/00_tcga_brca_scatac_locus_support.py` | 统计乳腺癌与健康乳腺 scATAC 峰集合中 CRC TF 位点的开放峰与基序支持（Fig. 3A–C；scATAC 层）。 | – |
| `sedb/01_dbsuper_breast_super_enhancers.py` | 将 dbSUPER 乳腺超级增强子坐标转换至 hg38，并统计 TF 位点 ±100 kb 内的重叠数。 | – |
| `tcga/03_atac_peak_summary.py` | 汇总 TCGA ATAC 峰矩阵在样本与峰维度的可及性统计。 | – |
| `tcga/04_atac_footprinting.py` | 读取 TCGA ATAC 足迹与基序表，计算足迹深度与 TF 表达相关性（footprint 层）。 | – |
| `tcga/06_tf_expression_activity_screen.py` | 构建 TF 表达矩阵与样本注释，执行表达/亚型关联检验，并基于 DoRothEA 计算调控子活性（RNA 与活性层；Fig. 1E, F 矩阵）。 | – |
| `tcga/07_atac_peak_annotation_subtype.py` | 对 ATAC 峰进行基因组注释并检验亚型特异性开放峰。 | – |
| `tcga/08_atac_candidate_tf_integration.py` | 将亚型特异性峰、峰-基因关联及基序富集整合至候选 TF 证据表。 | – |
| `tcga/09_cistrome_tf_tf_edges.py` | 将 Cistrome Cancer TF–TF 边整合至候选 TF 证据表。 | – |
| `tcga/10_enhancer_h3k27ac_evidence.py` | 整合 TCGA 增强子与 ENCODE MCF-7 H3K27ac 证据。 | – |
| `tcga/11_crc_member_evidence_scores.py` | 计算候选 TF 的 CRC 成员证据评分并筛选高置信候选因子。 | – |
| `tcga/12_tf_tf_edge_evidence_scores.py` | 为每条候选 TF–TF 边计算多维证据评分（Fig. 2 网络边证据表）。 | – |
| `tcga/14_crc_activity_scores.py` | 计算每个样本的 CRC TF 表达与活性模块评分。 | – |
| `tcga/17_clinical_outcome_dataset.py` | 合并 CRC 模块评分与 TCGA-BRCA 临床随访数据，运行临床关联及生存分析（log-rank、Cox）。 | – |
| `tf_resources/00_lambert_tf_list.py` | 从 Lambert et al. (2018) 数据库表中提取 1,639 个人类 TF 列表。 | – |
| `tf_resources/01_convert_dorothea_rda.R` | 将 DoRothEA 调控子 `.rda` 文件转换为 `.tsv` 格式。 | `[input.rda] [output.tsv]` |
| `trrust/00_trrust_edges.py` | 将 TRRUST v2 人类 TF–靶基因关系补充至 TF–TF 边证据表。 | – |

</details>

<details>
<summary><strong>阶段 1–5：正文核心分析流程（<code>pipeline/01_tf_scoring/</code> – <code>05_outcome_dependency/</code> — 13 个脚本）</strong></summary>

| 模块 / 脚本 | 功能说明 | 参数 |
|---|---|---|
| `01_tf_scoring/tf_evidence_scoring.py` | **Fig. 1A, Table S1：** 基于 7 个核心层与 3 个支持层对 1,639 个 TF 评分，按分层规则筛选 13 个 CRC TF 并分配亚型背景归属。 | `--raw-root --remap-result-dir --output-dir --intermediate-dir [--max-workers ...]` |
| `01_tf_scoring/prepare_subtype_association_inputs.py` | 生成 `subtype_association.py` 所需的 3 个输入文件（TF 表达矩阵、样本注释、CRC TF 表）。 | `[out_dir]` |
| `01_tf_scoring/subtype_association.py` | **Fig. 1E：** 对 1,631 个 TF 在 Luminal、HER2、Basal/TNBC 与 Normal-like 组间执行 Kruskal–Wallis 检验及 BH 校正，输出全部 TF 及 Top 20。 | `--input-dir DIR --output-dir DIR [--top-n 20]` |
| `02_chromatin/atac_rna_coupling.py` | **Fig. 3D, Table S6：** 计算每个 TF 位点 ±50 kb 范围内的 ATAC 可及性与自身表达量之间的相关性。 | – |
| `03_bulk_modules/module_scores_tumor_state_mcp.py` | **Fig. 6：** 计算 TCGA-BRCA 的 CRC 模块评分（表达与活性）及其与肿瘤状态特征、MCP-counter 评分在背景校正前后的相关性。 | – |
| `04_single_cell/01_process_and_score_cells.py` | **Fig. 7A：** 处理 GSE176078 单细胞数据，构建 UMAP 降维并计算单细胞 CRC TF 模块评分。 | – |
| `04_single_cell/02_state_discovery.py` | **Fig. 7B, C：** 在 ER+、HER2+ 与 TNBC 恶性细胞内对 TF 共表达聚类，定义 6 种恶性细胞状态及特征基因。 | – |
| `04_single_cell/03_state_scores_per_cell.py` | **Fig. 7C：** 计算每个恶性细胞的状态评分及主导状态归属。 | – |
| `05_outcome_dependency/tcga_state_scores.py` | **Fig. 8A：** 在 TCGA-BRCA 中计算 6 种状态评分并进行跨背景 Kruskal–Wallis 检验。 | – |
| `05_outcome_dependency/external_cohort_state_genes.py` | 将 6 种状态特征映射至 TCGA、METABRIC 与 GSE96058 队列并拟合连续变量 Cox 模型。 | – |
| `05_outcome_dependency/depmap_state_scores.py` | **Fig. 8C：** 计算 DepMap 乳腺癌细胞系的状态表达评分，并检验其与 CRISPR 基因敲除效应的相关性。 | – |
| `05_outcome_dependency/tcga_state_cox.py` | **Fig. 8B, Table S7：** 通过 `tcga_state_cox.R` 在 TCGA-BRCA 中拟合各状态评分的连续 Cox 模型（OS、PFS；按每 1 SD；120 个月截尾）。 | `[out_dir]` |
| `05_outcome_dependency/tcga_state_cox.R` | 由 `tcga_state_cox.py` 调用：使用 `survival` 包拟合 Cox 模型并执行比例风险假设检验。 | `<model_input.tsv> <results.tsv>` |

</details>

<details>
<summary><strong>阶段 6：同行评审补充分析（<code>revision_analyses/</code> — 29 个脚本）</strong></summary>

| 脚本 | 功能说明 | 参数 |
|---|---|---|
| `R1_first_round/R1_analyze_revision.py` | 首轮审稿综合分析：XBP1/ZNF703 筛选证据、单层剔除、权重扰动与随机集合对照、Fig. 3D 峰审计、Fig. 6A/8A 中位数、自身位点结合、BH 校正及临床协变量调整 Cox 模型、DepMap 匹配零模型与置换检验（表 `R1_1–R1_6`）。 | `[project_root] [out_dir]` |
| `R1_first_round/R1_external_projection.py` | 首轮审稿：将 6 种恶性状态特征映射至 METABRIC 与 GSE96058（表 `R1_7`；Fig. S4, Table S10）。 | `[project_root] [out_dir]` |
| `A0_reanalysis_checks/A0_verify.py` | 独立复算与核验：筛选规则、随机集合、Fig. 3D 统计、Cox 模型、DepMap 匹配零模型、外部队列映射及 Fig. 6A/8A 中位数。 | `<project_root> <out_dir>` |
| `A1_purity_adjustment/A1_purity_adjusted_immune.py` | **Fig. S3, Tables S9a–b：** 经肿瘤纯度（CPE）校正后的模块评分与 MCP-counter 相关性分析。 | `<project_root> <out_dir> [scores_tsv]` |
| `A7_Basal_TNBC_H3K27ac_layer/A7_basal_h3k27ac_audit.py` | 纳入 GSE245830 TNBC H3K27ac 峰，重算 Basal/TNBC 证据层与背景成员归属（Tables S1, S4d, S4e）。 | `<project_root> <extracted_GSE245830_dir> <out_dir>` |
| `A8_selection_robustness/A8_selection_robustness_corrected.py` | **Fig. S1, Tables S4a–c：** 基于修订后评分表的筛选稳健性检验（权重扰动、逐层剔除、阈值偏移、随机集合对照）。 | `<A7_corrected_csv> <out_dir>` |
| `A10_Fig6_modules/A10_fig6_modules_corrected.py` | **Tables S9c–d：** 使用修订后 Luminal 成员集合（剔除 RUNX1）重算 Fig. 6 模块评分。 | `<project_root> <A7_corrected_csv> <out_dir>` |
| `A12_Fig3D_mean_accessibility_script/figure3d_corrected.py` | 修订后的 Fig. 3D 统计：局部平均可及性与表达量相关性，以及覆盖全部局部峰检验的 BH 校正。 | `--input-dir DIR --output-dir DIR [--flank-bp 50000]` |
| `A13_motif_context_membership/A13_motif_membership.py` | **Table S14：** 按节点的亚型背景归属对 Fig. 2C 网络基序进行分类汇总。 | `<figure2c.csv> <A7_CRC_TF_calls_stored_vs_corrected.tsv> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_02_pgr_target_and_network.py` | 将 PGR 纳入结合靶基因，加入 ESR1→PGR 重建网络边集并枚举调控基序。 | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_03_chipatlas_crosscheck.py` | 利用本地 ChIP-Atlas 乳腺峰文件交叉验证 13 个 TF 的启动子结合情况。 | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_04_knockdown_label_sensitivity.py` | 评估剔除 5 个元数据正则遗漏的敲低（knockdown）标签 tier-1 数据集对结果的影响。 | `<project_root> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_05_tables_and_control.py` | 生成更新版结合矩阵（ReMap 的 PGR 列、ChIP-Atlas 的 PGR 行）、自身结合及网络规则对照表。 | `<project_root> <S5_xlsx> <out_dir>` |
| `A14_occupancy_and_network_checks/A14_06_edges_table.py` | 整理证据分 ≥ 6 的网络边汇总表（Table S14b）。 | `<out_dir>` |
| `A14_occupancy_and_network_checks/A14_09_pgr_row_binding_knockdown.py` | 基于筛选后的 ReMap 与 ChIP-Atlas PGR 数据集（清单见 `data/dataset_lists/`）生成最终结合表，并评估 PGR 结合分对筛选结果的影响。 | `<project_root> <Table_S5_workbook> <out_dir>` |
| `A15_selection_coverage_audit/A15_selection_bias_audit.py` | 只读式筛选偏差审计：各层覆盖率、稀疏层归属、去除稀疏层后的排名、一致足迹证据、替代评分规则重放及有乳腺 ChIP 但结合层为 0 的 TF。 | `<project_root> <revision_analyses_dir> <out_dir> <experimentList.tab>` |
| `A16_supplementary_ChIP/scripts/a16_common.py` | A16 共享路径、文件读写、校验和、启动子表及峰重叠辅助函数。 | 被 A16 脚本导入 |
| `A16_supplementary_ChIP/scripts/a16_liftover.py` | *（顺序 1）* 基于 chain 文件的区间坐标转换（单 chain 覆盖 ≥ 95% 碱基），并与 `pyliftover` 交叉验证。 | – |
| `A16_supplementary_ChIP/scripts/a16_bootstrap.py` | *（顺序 2）* 复现初稿中的 PGR 基线并输出启动子窗口表。 | – |
| `A16_supplementary_ChIP/scripts/a16_analyze_peaks.py` | *（顺序 3）* 复现初版结合矩阵并按数据来源分析补充 ChIP 峰（清单见 `data/dataset_lists/`）。 | – |
| `A16_supplementary_ChIP/scripts/a16_signal_quantify.mjs` | *（顺序 4）* 直接在 TF 启动子处定量 FOXC1 全基因组 bigWig 连续信号（不经峰识别；需 Node.js 与 `@gmod/bbi`）。 | `node a16_signal_quantify.mjs <project_root>` |
| `A16_supplementary_ChIP/scripts/a16_quality.py` | *（顺序 5）* 峰文件区间质量审计及不同处理流程间的一致性核对。 | – |
| `A16_supplementary_ChIP/scripts/a16_sensitivity.py` | *（顺序 6）* **Table S4g：** 仅替换结合证据层后对原 13-TF 输入的条件重放测试。 | – |
| `A17_score5_motifs/A17_01_score5_motifs.py` | **Tables S14d–e：** 检验放宽纳入证据分为 5 的边后，HER2 或 Basal/TNBC 集合是否形成闭合环路。 | `<A16_results_dir> <out_dir>` |
| `A18_selection_and_subtype_checks/A18_selection_subtype_checks.py` | **Tables S2c–d, S4h, S4k–l：** 针对审稿意见 R1.5 及 R5.2–R5.5 的专项核查（Fig. 1E 排名与原因、五层参考、均衡亚型构成、PAM50 与受体分组对照、资源复用核查）。 | `<project_root> <revision_analyses_dir> <A16_dir> <out_dir>` |
| `A19_occupancy_network_self_locus/A19_occupancy_network.py` | **Figs. 2, 4C; Tables S5, S14：** 纳入补充 ChIP 后的结合矩阵、TF 调控网络与基序分析，以及未用于 TF 筛选的独立研究复现检验。 | `<project_root> <revision_analyses_dir> <A16_dir> <out_dir>` |
| `A19_occupancy_network_self_locus/A19_02_self_locus_100kb.py` | **Fig. S2B, Table S5h：** TSS ±100 kb 内的自身位点结合强度，并在全部蛋白编码基因同等窗口中排名。 | `<project_root> <A16_dir> <out_dir>` |
| `A21_RNA_grouping_subsampling/A21_rna_grouping_balanced.py` | **Tables S4i, S4i2：** RNA 层亚型检验分组定义及其在均衡抽样与按比例抽样下的重算。 | `<project_root> <revision_analyses_dir> <out_dir>` |
| `A22_binding_by_cell_line/A22_cell_line_binding.py` | **Tables S14f–g：** 追溯每条网络边结合证据对应的细胞系来源，并检验单一细胞系内的闭合环路结构。 | `<A16_dir> <A19_out_dir> <out_dir>` |

</details>

---

## 引用格式

如您在研究中使用了本仓库代码，请引用对应论文（亦可直接使用 GitHub 的 `Cite this repository` 按钮读取 [`CITATION.cff`](CITATION.cff)）：

```bibtex
@article{Zeng2026BRCACRC,
  author  = {Zeng, Yongjia and Shi, Songfu and Wang, Fan and Zhao, Lixin and Ye, Lingxi and Wu, Xiaodi and Yan, Dong},
  title   = {Subtype-specific Core Transcriptional Regulatory Circuitry Links Transcriptional States, Chromatin Context, and Functional Dependencies in Breast Cancer},
  journal = {Functional \& Integrative Genomics},
  year    = {2026},
  note    = {Under review}
}
```

## 开源许可

本仓库代码基于 [MIT License](LICENSE) 开源发布。所涉及的第三方公共数据集与软件工具仍遵循其各自的使用条款与许可协议，使用时请按论文规范引用原始文献。
