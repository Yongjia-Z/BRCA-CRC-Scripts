# BRCA-CRC-Scripts
Analysis code (Python) for the core transcriptional regulatory circuit in breast cancer

# Breast Cancer CRC Atlas Figure Scripts

## Overview

This repository contains the figure-level analysis scripts used to reproduce the processed data tables underlying Figures 1-8 of this study on subtype-aware core regulatory circuitry in breast cancer.

Only the analysis scripts are provided here. Raw public datasets, controlled-access clinical data, large sequencing files, processed intermediate matrices, manuscript files, and final edited figure files are not included. Users should download the required datasets from the original sources described in the manuscript and supplementary materials, then update the input paths or configuration blocks in each script before running the analysis.

The scripts are designed to reproduce processed analysis tables used for figure generation, rather than to produce final publication-ready artwork. Final figure styling, panel assembly, and manual vector editing were performed separately.

## Repository Contents

```text
.
├── figure1.py
├── figure2.py
├── figure3.py
├── figure4.py
├── figure5.py
├── figure6.py
├── figure7.py
├── figure8.py
└── README.md
```

## What Is Included

- Data-processing code for Figures 1-8.
- Figure-level analysis workflows from input data to processed output tables.
- Evidence integration, subtype-aware scoring, network analysis, chromatin-accessibility analysis, bulk cohort validation, single-cell analysis summaries, clinical association analysis, and functional dependency analysis.
- Reproducibility notes for dataset usage and interpretation.

## What Is Not Included

- Raw sequencing data.
- Large public expression, ATAC-seq, ChIP-seq, CUT&Tag, Hi-C, single-cell, or CRISPR dependency matrices.
- Controlled-access clinical data.
- Final edited figure panels.
- Manuscript files.
- Local absolute paths.
- Software environments tied to a specific user machine.

## Required Input Data

The scripts assume that the required datasets have already been downloaded locally. Depending on the figure, these may include:

- Bulk breast cancer transcriptomic datasets.
- Breast cancer molecular subtype annotations.
- Clinical outcome and pathological annotation tables.
- TF activity or regulon activity matrices.
- Bulk ATAC-seq or chromatin accessibility peak files.
- H3K27ac, enhancer, or super-enhancer proxy evidence.
- Breast-context TF binding evidence from ChIP-seq, CUT&Tag, ReMap, ChIP-Atlas, Cistrome, GEO, or related resources.
- Hi-C or enhancer-promoter contact evidence when available.
- Single-cell RNA-seq expression matrices and cell annotations.
- DepMap or Project Achilles CRISPR dependency matrices.
- External validation cohorts such as METABRIC and large breast cancer expression cohorts.

Exact dataset accessions, versions, and usage context should be checked against the manuscript, figure legends, and supplementary tables.

## Important Scientific Notes

### CRC TF Input

The CRC TF candidate set used in the figure scripts is supplied as an external input table or configuration file. The gene list is not hard-coded in this README to avoid version mismatch between documentation, scripts, and manuscript tables.

### Breast-Context Evidence

TF binding evidence should be interpreted in a breast cancer or breast epithelial context whenever possible. Pan-context ReMap or other all-tissue database evidence should not be treated as direct breast cancer evidence unless explicitly labeled as auxiliary.

### Enhancer and Super-Enhancer Evidence

H3K27ac and enhancer-based signals are treated as enhancer or super-enhancer proxy evidence unless formal ROSE, LILY, SEdb, dbSUPER, SEA, or equivalent super-enhancer calling is performed.

### Chromatin Contact Evidence

Hi-C, promoter-capture Hi-C, HiChIP, ChIA-PET, or ABC-like links represent physical or predicted enhancer-promoter contact evidence. These data should not be interpreted as TF-TF regulatory edges.

### Single-Cell Evidence

Single-cell RNA-seq analyses in these scripts are used to localize CRC TF expression and module scores across malignant epithelial and non-malignant cell states. Unless explicitly stated, they do not represent SCENIC, chromVAR, or single-cell ATAC-derived regulon activity.

### Survival Analysis

Clinical association and survival analyses are intended to test whether CRC-related expression states have clinical relevance. They are not intended to build or claim a clinical risk prediction model.

## Software Requirements

Recommended environment:

```text
Python >= 3.10
pandas
numpy
scipy
matplotlib
seaborn
scikit-learn
statsmodels
lifelines
openpyxl
networkx
```

Additional packages may be required for specific figures:

```text
scanpy
anndata
pyranges
pyBigWig
pybedtools
h5py
tqdm
```

For genomic interval overlap analysis, optimized genomic interval libraries such as `pyranges` are recommended. Avoid row-wise interval matching with slow table iteration methods.

## Installation Example

Create a clean Python environment:

```bash
conda create -n breast-crc-atlas python=3.10
conda activate breast-crc-atlas
pip install pandas numpy scipy matplotlib seaborn scikit-learn statsmodels lifelines openpyxl networkx tqdm
```

Install optional packages if needed:

```bash
pip install scanpy anndata pyranges pyBigWig h5py
```

Some genomic tools may require installation through conda:

```bash
conda install -c bioconda bedtools
```

## How to Run

Each script corresponds to one main figure. Before running a script, update the user-editable input and output paths near the top of the script, or provide paths through the command-line interface if the script supports command-line arguments.

Typical usage:

```bash
python figure1.py
python figure2.py
python figure3.py
python figure4.py
python figure5.py
python figure6.py
python figure7.py
python figure8.py
```

If a script supports command-line arguments, a typical pattern is:

```bash
python figure1.py --input-dir /path/to/input_data --output-dir /path/to/results
```

Outputs are usually written as `.csv`, `.tsv`, `.json`, or intermediate table files that can be used for figure plotting and downstream manual figure assembly.

## Figure-Level Description

### Figure 1

Multi-evidence discovery and cohort contextualization of breast cancer CRC TF candidates.

Main analysis components include:

- Integration of RNA expression evidence.
- TF activity or regulon-level evidence.
- ATAC-seq, motif, footprint, or accessibility evidence.
- Subtype-aware enhancer or H3K27ac proxy evidence.
- Breast-context direct TF binding evidence.
- External expression validation in independent breast cancer cohorts.
- Evidence-intersection summaries and subtype-associated TF ranking.

### Figure 2

Subtype-aware CRC network topology and TF-TF regulatory evidence.

Main analysis components include:

- Construction of CRC TF network nodes and edges.
- Calculation of network centrality metrics.
- Evidence-based TF-TF edge scoring.
- Subtype-aware circuit or loop motif summaries.
- Separation of regulatory network evidence from physical chromatin contact evidence.

### Figure 3

Chromatin accessibility and regulatory potential around CRC TF loci.

Main analysis components include:

- CRC TF locus-level accessibility summaries.
- Peak-count or accessibility burden analysis.
- Genomic annotation of open chromatin regions.
- Motif-family burden analysis.
- ATAC-RNA coupling analysis where matched or comparable data are available.

### Figure 4

Breast-context TF binding and chromatin evidence.

Main analysis components include:

- Filtering TF binding evidence to breast cancer or breast epithelial contexts when possible.
- Incorporation of TNBC-specific ChIP-seq, CUT&Tag, or related binding datasets when available.
- Quantitative overlap analysis between TF binding, chromatin accessibility, enhancer evidence, and CRC TF loci.
- Documentation of data gaps where direct subtype-matched binding evidence is unavailable.

### Figure 5

Locus-level chromatin and 3D genome support.

Main analysis components include:

- Browser-track level summaries for selected CRC TF loci.
- Continuous signal track processing when bigWig files are available.
- Peak-based evidence summaries when only BED-like peak files are available.
- Hi-C or contact-matrix evidence for selected loci.
- Context-aware interpretation of cell-line matched and non-matched 3D genome data.

### Figure 6

Bulk tumor-state projection of CRC TF modules.

Main analysis components include:

- CRC TF module score calculation in bulk breast cancer cohorts.
- Association with molecular subtype.
- Correlation with tumor-state, immune, stromal, epithelial, or microenvironment-related signatures.
- Adjustment or stratification analyses when appropriate.

### Figure 7

Clinical, external cohort, and functional validation.

Main analysis components include:

- External validation of CRC TF module or state scores.
- Cohort-level comparison across independent breast cancer datasets.
- Survival and Cox association analyses where clinical endpoints are available.
- DepMap or Project Achilles CRISPR dependency analysis.
- Assessment of whether subtype-specific CRC states show reproducible clinical or functional relevance.

### Figure 8

Single-cell localization and CRC state decomposition.

Main analysis components include:

- Single-cell RNA-seq preprocessing or reuse of processed single-cell matrices.
- Cell-type and malignant epithelial cell annotation.
- CRC TF expression localization across cell states.
- Malignant-cell subcluster analysis.
- Projection of single-cell-derived CRC states back to bulk cohorts when applicable.

## Output Files

The scripts generate processed analysis outputs such as:

```text
figure1_processed_table.csv
figure2_network_edges.csv
figure3_accessibility_summary.csv
figure4_binding_overlap_summary.csv
figure5_locus_contact_summary.csv
figure6_bulk_module_scores.csv
figure7_survival_cox_summary.csv
figure8_single_cell_crc_state_summary.csv
```

Actual file names may differ by script. Users should check each script's output section.

## Reproducibility Recommendations

To improve reproducibility:

1. Record dataset accession numbers and download dates.
2. Keep genome builds explicit, especially when combining hg19 and hg38 resources.
3. Use liftover when integrating datasets from different genome builds.
4. Store raw downloaded files separately from processed outputs.
5. Keep figure-ready tables separate from final edited artwork.
6. Do not mix all-context TF binding evidence with subtype-matched breast cancer evidence without explicit labeling.
7. Record missing evidence rather than silently removing genes or datasets.
8. Use optimized interval-overlap methods for genomic coordinate analyses.

## Data Availability

All datasets used by the scripts should be obtained from their original public repositories or access-controlled portals. Users are responsible for complying with the data-use terms of each resource.

This repository does not redistribute raw or controlled-access data.

## AI-Assisted Code Development

AI-based coding tools were used to assist script writing, code refactoring, and documentation preparation. All analytical decisions, dataset selection, biological interpretation, code review, and final validation were performed by the authors.

## Citation

If you use these scripts, please cite the associated manuscript and the original data resources used in each analysis.

---

# 乳腺癌 CRC 图表分析脚本

## 项目简介

本仓库包含本研究 Figure 1-8 所用的图级分析脚本，用于复现乳腺癌亚型相关核心转录调控环路分析中的原始结果表格。

本仓库仅提供分析脚本。原始公共数据、受控访问临床数据、大型测序文件、处理后的中间矩阵、手稿文件以及最终精修图片均不包含在本仓库中。使用者需要根据正文和补充材料中列出的数据来源，自行从原始数据库下载所需数据，并在运行脚本前修改各脚本中的输入路径或配置参数。

这些脚本的目标是复现用于作图的分析结果表格，而不是直接生成最终投稿用的精修图片。最终图片排版、面板组合和矢量图细节调整是在脚本输出结果之后单独完成的。

## 仓库内容

```text
.
├── figure1.py
├── figure2.py
├── figure3.py
├── figure4.py
├── figure5.py
├── figure6.py
├── figure7.py
├── figure8.py
└── README.md
```

## 包含内容

- Figure 1-8 的数据处理代码。
- 从输入数据到图表结果表格的分析流程。
- 多组学证据整合、亚型相关评分、网络分析、染色质开放性分析、bulk 队列验证、单细胞分析摘要、临床相关分析和功能依赖分析。
- 数据使用和结果解释的可重复性说明。

## 不包含内容

- 原始测序数据。
- 大型公共表达矩阵、ATAC-seq、ChIP-seq、CUT&Tag、Hi-C、单细胞或 CRISPR dependency 矩阵。
- 受控访问临床数据。
- 最终精修图片。
- 手稿文件。
- 本地绝对路径。
- 绑定特定个人电脑的软件环境。

## 所需输入数据

脚本默认使用者已经在本地下载好所需数据。根据不同图，可能需要以下数据：

- 乳腺癌 bulk 转录组数据。
- 乳腺癌分子亚型注释。
- 临床结局和病理注释表。
- TF activity 或 regulon activity 矩阵。
- Bulk ATAC-seq 或染色质开放峰文件。
- H3K27ac、enhancer 或 super-enhancer proxy 证据。
- 来自 ChIP-seq、CUT&Tag、ReMap、ChIP-Atlas、Cistrome、GEO 或相关资源的乳腺癌背景 TF binding 证据。
- 可用的 Hi-C 或 enhancer-promoter contact 证据。
- 单细胞 RNA-seq 表达矩阵和细胞注释。
- DepMap 或 Project Achilles CRISPR dependency 矩阵。
- METABRIC 和其他大型乳腺癌外部验证队列。

具体数据 accession、版本和使用场景应以正文、图注和补充表为准。

## 重要科学说明

### CRC TF 输入

图表脚本使用的 CRC TF 候选集合应作为外部输入表格或配置文件提供。为避免 README、脚本和正文表格之间出现版本不一致，本 README 不直接写入具体基因名单。

### 乳腺癌背景证据

TF binding 证据应尽可能限定在乳腺癌或乳腺上皮背景中解释。Pan-context ReMap 或其他跨组织数据库证据不能直接等同于乳腺癌直接证据，除非明确标注为辅助证据。

### Enhancer 和 Super-Enhancer 证据

在未进行正式 ROSE、LILY、SEdb、dbSUPER、SEA 或等价 super-enhancer calling 的情况下，H3K27ac 和 enhancer 相关结果应解释为 enhancer 或 super-enhancer proxy evidence。

### 染色质空间接触证据

Hi-C、promoter-capture Hi-C、HiChIP、ChIA-PET 或 ABC-like links 表示物理或预测的 enhancer-promoter contact 证据。这类数据不能直接解释为 TF-TF 调控边。

### 单细胞证据

这些脚本中的单细胞 RNA-seq 分析主要用于定位 CRC TF 表达和 module score 在恶性上皮细胞及非恶性细胞状态中的分布。除非特别说明，这些结果不代表 SCENIC、chromVAR 或单细胞 ATAC 推断的 regulon activity。

### 生存分析

临床相关和生存分析用于检验 CRC 相关表达状态是否具有临床相关性，不用于构建或声称临床风险预测模型。

## 软件需求

推荐环境：

```text
Python >= 3.10
pandas
numpy
scipy
matplotlib
seaborn
scikit-learn
statsmodels
lifelines
openpyxl
networkx
```

部分图可能需要额外软件包：

```text
scanpy
anndata
pyranges
pyBigWig
pybedtools
h5py
tqdm
```

对于基因组区间重叠分析，推荐使用 `pyranges` 等专门的基因组区间处理工具。应避免使用低效的逐行表格循环方法进行大规模区间匹配。

## 安装示例

创建干净的 Python 环境：

```bash
conda create -n breast-crc-atlas python=3.10
conda activate breast-crc-atlas
pip install pandas numpy scipy matplotlib seaborn scikit-learn statsmodels lifelines openpyxl networkx tqdm
```

如需运行单细胞或基因组区间相关脚本，可安装额外依赖：

```bash
pip install scanpy anndata pyranges pyBigWig h5py
```

部分基因组工具可通过 conda 安装：

```bash
conda install -c bioconda bedtools
```

## 运行方式

每个脚本对应一张主要图。在运行脚本前，请修改脚本开头的用户配置区域中的输入和输出路径；如果脚本支持命令行参数，也可以通过命令行传入路径。

基本运行方式：

```bash
python figure1.py
python figure2.py
python figure3.py
python figure4.py
python figure5.py
python figure6.py
python figure7.py
python figure8.py
```

如果脚本支持命令行参数，常见形式为：

```bash
python figure1.py --input-dir /path/to/input_data --output-dir /path/to/results
```

输出通常为 `.csv`、`.tsv`、`.json` 或其他中间结果表格，可用于后续作图和人工排版。

## 各图脚本说明

### Figure 1

乳腺癌 CRC TF 候选因子的多组学证据发现与队列背景解释。

主要分析内容包括：

- RNA 表达证据整合。
- TF activity 或 regulon 层面证据。
- ATAC-seq、motif、footprint 或染色质开放性证据。
- 亚型感知的 enhancer 或 H3K27ac proxy 证据。
- 乳腺癌背景下的直接 TF binding 证据。
- 外部乳腺癌队列中的表达验证。
- 证据交集总结和亚型相关 TF 排名。

### Figure 2

亚型相关 CRC 网络拓扑和 TF-TF 调控证据。

主要分析内容包括：

- CRC TF 网络节点和边的构建。
- 网络中心性指标计算。
- 基于多层证据的 TF-TF edge scoring。
- 亚型相关 circuit 或 loop motif 总结。
- 将调控网络证据与物理染色质接触证据分开解释。

### Figure 3

CRC TF 基因座附近的染色质开放性和调控潜能。

主要分析内容包括：

- CRC TF locus-level accessibility 总结。
- Peak count 或 accessibility burden 分析。
- 开放染色质区域的基因组注释。
- Motif-family burden 分析。
- 在可用数据中进行 ATAC-RNA coupling 分析。

### Figure 4

乳腺癌背景 TF binding 和染色质证据。

主要分析内容包括：

- 尽可能将 TF binding 证据限制在乳腺癌或乳腺上皮背景。
- 在可用时整合 TNBC 特异的 ChIP-seq、CUT&Tag 或相关 binding 数据。
- 定量分析 TF binding、染色质开放性、enhancer 证据和 CRC TF loci 的重叠。
- 对缺少直接亚型匹配 binding 证据的情况进行明确记录。

### Figure 5

CRC TF 位点级染色质和 3D genome 支持证据。

主要分析内容包括：

- 选定 CRC TF loci 的浏览器轨迹总结。
- 在可获得 bigWig 文件时处理连续信号轨道。
- 仅有 BED-like peak 文件时进行 peak-based 证据总结。
- 对选定位点进行 Hi-C 或 contact matrix 证据分析。
- 区分细胞系匹配和非匹配 3D genome 数据的解释范围。

### Figure 6

Bulk 肿瘤状态中的 CRC TF module 投射。

主要分析内容包括：

- 在 bulk 乳腺癌队列中计算 CRC TF module score。
- 分析其与分子亚型的关系。
- 分析其与肿瘤状态、免疫、基质、上皮或微环境相关 signature 的关系。
- 在适用情况下进行调整或分层分析。

### Figure 7

临床、外部队列和功能验证。

主要分析内容包括：

- CRC TF module 或 state score 的外部验证。
- 独立乳腺癌数据集之间的队列层面比较。
- 在具有临床结局的数据集中进行 survival 和 Cox association 分析。
- DepMap 或 Project Achilles CRISPR dependency 分析。
- 评估亚型相关 CRC 状态是否具有可重复的临床或功能意义。

### Figure 8

单细胞定位和 CRC state 分解。

主要分析内容包括：

- 单细胞 RNA-seq 数据预处理或复用已处理矩阵。
- 细胞类型和恶性上皮细胞注释。
- CRC TF 表达在不同细胞状态中的定位。
- 恶性细胞亚群分析。
- 在适用情况下，将单细胞来源的 CRC state 回投到 bulk 队列。

## 输出文件

脚本会生成用于作图和后续分析的结果表格，例如：

```text
figure1_processed_table.csv
figure2_network_edges.csv
figure3_accessibility_summary.csv
figure4_binding_overlap_summary.csv
figure5_locus_contact_summary.csv
figure6_bulk_module_scores.csv
figure7_survival_cox_summary.csv
figure8_single_cell_crc_state_summary.csv
```

实际文件名可能因脚本而异，请以各脚本的输出部分为准。

## 可重复性建议

为提高可重复性，建议：

1. 记录数据集 accession 和下载日期。
2. 明确记录 genome build，尤其是在整合 hg19 和 hg38 数据时。
3. 跨 genome build 整合时进行 liftover 并记录方法。
4. 将原始下载文件和处理后输出分开保存。
5. 将 figure-ready tables 与最终精修图片分开保存。
6. 不要在未标注的情况下混用 all-context TF binding 证据和乳腺癌匹配证据。
7. 对缺失证据进行明确记录，而不是静默删除基因或数据集。
8. 基因组坐标分析应使用优化过的区间重叠算法。

## 数据可用性

本研究使用的数据应从原始公共数据库或受控访问平台获取。使用者需遵守各数据资源的使用条款。

本仓库不重新分发原始数据或受控访问数据。

## AI 辅助代码开发声明

本研究在脚本编写、代码重构和文档整理过程中使用了 AI 工具进行辅助。所有分析设计、数据集选择、生物学解释、代码审查和最终结果验证均由作者完成。

## 引用

如果使用本仓库中的脚本，请引用对应研究论文以及各分析所使用的原始数据资源。
