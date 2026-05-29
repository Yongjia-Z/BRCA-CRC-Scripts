# BRCA CRC Source-Table Scripts

This repository contains the cleaned Figure 1-8 source-table generation scripts for the breast cancer core transcriptional regulatory circuitry manuscript.

## Scope

The repository is intended to reproduce numerical source tables used for manuscript figures. It does not include final figure polishing code, submitted PDF/SVG/TIFF assets, raw data matrices, or large genome-track files.

Included:

- `figure1/result_scripts` to `figure8/result_scripts`: panel-level Python scripts.
- `_repro_common/crc_repro.py`: shared helper functions used by several scripts.
- `result_script_pack_manifest.csv`: panel-to-script mapping and expected output table names.
- `data_manifest.tsv`: expected input files, public sources/accessions, formats, and redistribution notes.
- `software_versions.tsv`: software versions checked during local audit.

Not included:

- raw TCGA, METABRIC, GEO, DepMap, ENCODE, ReMap, bigWig, Hi-C, or single-cell matrices;
- final figure rendering and graphical polishing scripts;
- local debug folders, temporary files, `opus`, bytecode caches, and generated outputs.

## Layout

```text
.
|-- README.md
|-- data_manifest.tsv
|-- result_script_pack_manifest.csv
|-- software_versions.tsv
|-- _repro_common/
|   `-- crc_repro.py
|-- figure1/
|   `-- result_scripts/
|-- figure2/
|   `-- result_scripts/
...
`-- figure8/
    `-- result_scripts/
```

## Running Scripts

Prepare panel-specific inputs under the `data/<panel>` paths described in `data_manifest.tsv`, then run individual scripts from the repository root. For scripts importing `_repro_common`, add it to `PYTHONPATH`.

Example:

```bash
PYTHONPATH=_repro_common python figure6/result_scripts/figure6a.py --input-dir data/figure6a --output-dir outputs/figure6a
```

Figure 1a is the discovery-stage script and requires a broader raw-data root:

```bash
python figure1/result_scripts/figure1a.py --raw-root /path/to/datasets --remap-result-dir /path/to/remap_results --output-dir outputs/figure1a --intermediate-dir outputs/figure1a/intermediate
```

## Notes

- Figure 1a reads a processed Lambert et al. human transcription-factor list, `Lambert_human_TF_list.tsv`, and does not contain a predefined final CRC TF list.
- DoRothEA-related outputs are regulon-derived TF activity proxies, not a claim of official VIPER/decoupleR processing.
- Figure 5a bigWig scripts summarize continuous signal tracks. If `pyBigWig` or referenced bigWig files are unavailable, NaN values indicate missing or unreadable track signal rather than measured zero signal.
- Figure 8b uses `Rscript` and the R `survival` package for continuous Cox models.
- Survival analyses are descriptive association analyses and should not be interpreted as optimized risk models or nomograms.

## Citation

Please cite the associated manuscript once available. Third-party datasets and software should be cited according to their original source terms.
