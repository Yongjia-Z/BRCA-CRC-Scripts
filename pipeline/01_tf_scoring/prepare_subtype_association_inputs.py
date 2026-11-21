import os
import shutil
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
RES = DATA_ROOT / "TCGA" / "results"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / "source_tables" / "figure1" / "subtype_association_inputs"
OUT.mkdir(parents=True, exist_ok=True)
pd.read_csv(RES / "TCGA_BRCA_TF_expression_matrix.tsv.gz", sep="\t").to_csv(OUT / "tcga_tf_expression_matrix.csv.gz", index=False, compression="gzip")
pd.read_csv(RES / "TCGA_BRCA_expression_sample_annotation.tsv", sep="\t").to_csv(OUT / "tcga_sample_annotation.csv", index=False)
shutil.copy2(REPO_ROOT / "source_tables" / "figure1" / "figure1a.csv", OUT / "figure1a.csv")
print("inputs written to", OUT)
