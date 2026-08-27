# =============================================================================
# 05b_run_pyscenic.sh  —  NUDT16 KO female spleen scRNA-seq
# =============================================================================
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/Users/}"
SCENIC_DIR="${SCENIC_DIR:-$PROJECT_DIR/results/}"

# ---- RESOURCES  -------------------------------------------------
# Mouse mm10 mc9nr cisTarget databases + TF list + motif annotations.
# Download once from https://resources.aertslab.org/cistarget/
RESOURCES_DIR="${RESOURCES_DIR:-$HOME/scenic_resources/mm10}"
TF_LIST="${TF_LIST:-$RESOURCES_DIR/allTFs_mm.txt}"
MOTIF_ANNOT="${MOTIF_ANNOT:-$RESOURCES_DIR/motifs-v9-nr.mgi-m0.001-o0.0.tbl}"

# Ranking databases 
RANKING_DBS="${RANKING_DBS:-\
$RESOURCES_DIR/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
$RESOURCES_DIR/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather}"

# ---- compute ----------------------------------------------------------------
NWORKERS="${NWORKERS:-8}"
SEED="${SEED:-217}"

LOOM="$SCENIC_DIR/female_counts.loom"
ADJ="$SCENIC_DIR/adjacencies.tsv"
REG="$SCENIC_DIR/regulons.csv"
AUC_LOOM="$SCENIC_DIR/aucell.loom"

cd "$SCENIC_DIR"

echo "==============================================================="
echo " pySCENIC run"
echo "   SCENIC_DIR   = $SCENIC_DIR"
echo "   TF_LIST      = $TF_LIST"
echo "   MOTIF_ANNOT  = $MOTIF_ANNOT"
echo "   RANKING_DBS  = $RANKING_DBS"
echo "   workers      = $NWORKERS   seed = $SEED"
echo "==============================================================="

# ---- ----------------------------------------------------------
command -v pyscenic >/dev/null 2>&1 || { echo "ERROR: pyscenic not on PATH."; exit 1; }
[ -f "$TF_LIST" ]     || { echo "ERROR: TF list not found: $TF_LIST"; exit 1; }
[ -f "$MOTIF_ANNOT" ] || { echo "ERROR: motif annotations not found: $MOTIF_ANNOT"; exit 1; }
for db in $RANKING_DBS; do
  [ -f "$db" ] || { echo "ERROR: ranking DB not found: $db"; exit 1; }
done

# ------------------
if [ ! -f "$LOOM" ]; then
  echo "[0] Building loom from expr_counts.mtx ..."
  python3 - "$SCENIC_DIR" <<'PY'
import sys, os
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import loompy as lp

d = sys.argv[1]

m = sio.mmread(os.path.join(d, "expr_counts.mtx")).tocsc()   # genes x cells
genes = [l.strip() for l in open(os.path.join(d, "genes.txt"))]
cells = [l.strip() for l in open(os.path.join(d, "barcodes.txt"))]
assert m.shape[0] == len(genes) and m.shape[1] == len(cells), "shape mismatch"
row_attrs = {"Gene": np.array(genes)}
col_attrs = {"CellID": np.array(cells)}
lp.create(os.path.join(d, "female_counts.loom"), m, row_attrs, col_attrs)
print("  wrote female_counts.loom:", m.shape[0], "genes x", m.shape[1], "cells")
PY
else
  echo "[0] Using existing loom: $LOOM"
fi

# ---- 1. GRNBoost2: co-expression modules (TF -> target adjacencies) ---------
echo "[1] pyscenic grn (GRNBoost2) ..."
pyscenic grn \
  "$LOOM" \
  "$TF_LIST" \
  --output "$ADJ" \
  --num_workers "$NWORKERS" \
  --seed "$SEED" \
  --method grnboost2

# ---- 2. cisTarget: prune modules into motif-enriched regulons ---------------
echo "[2] pyscenic ctx (cisTarget pruning) ..."
pyscenic ctx \
  "$ADJ" \
  $RANKING_DBS \
  --annotations_fname "$MOTIF_ANNOT" \
  --expression_mtx_fname "$LOOM" \
  --output "$REG" \
  --num_workers "$NWORKERS" \
  --mask_dropouts

# ---- 3. AUCell: per-cell regulon activity -----------------------------------
echo "[3] pyscenic aucell ..."
pyscenic aucell \
  "$LOOM" \
  "$REG" \
  --output "$AUC_LOOM" \
  --num_workers "$NWORKERS" \
  --seed "$SEED"

# ---- 4. Export tidy CSVs for the R stage (05c) ------------------------------
echo "[4] Exporting auc_mtx.csv + regulon_targets.csv ..."
python3 - "$SCENIC_DIR" <<'PY'
import sys, os, re
import numpy as np
import pandas as pd
import loompy as lp

d = sys.argv[1]

# 4a. AUCell matrix: cells x regulons
with lp.connect(os.path.join(d, "aucell.loom"), mode="r", validate=False) as ds:
    auc = pd.DataFrame(
        ds.ca["RegulonsAUC"] if "RegulonsAUC" in ds.ca.keys() else None
    )
    cells = ds.ca["CellID"]
if auc is None or auc.shape[1] == 0:
    with lp.connect(os.path.join(d, "aucell.loom"), mode="r", validate=False) as ds:
        rec = ds.ca["RegulonsAUC"]
        auc = pd.DataFrame(rec)
auc.index = cells
auc.index.name = "CellID"
auc.to_csv(os.path.join(d, "auc_mtx.csv"))
print("  wrote auc_mtx.csv:", auc.shape[0], "cells x", auc.shape[1], "regulons")

reg = pd.read_csv(os.path.join(d, "regulons.csv"), header=[0, 1], index_col=[0, 1])

tf_level = reg.index.get_level_values(0)

tgt_col = [c for c in reg.columns if c[1] == "TargetGenes"]
rows = []
if tgt_col:
    tgt = reg[tgt_col[0]]
    for tf, cell in zip(tf_level, tgt):
        if not isinstance(cell, str):
            continue
        for g in re.findall(r"\('([^']+)',\s*([0-9.eE+-]+)\)", cell):
            rows.append((f"{tf}(+)", tf, g[0], float(g[1])))
tt = pd.DataFrame(rows, columns=["regulon", "TF", "target", "weight"])

if len(tt):
    tt = (tt.sort_values("weight", ascending=False)
            .drop_duplicates(["regulon", "target"]))
tt.to_csv(os.path.join(d, "regulon_targets.csv"), index=False)
print("  wrote regulon_targets.csv:", tt.shape[0], "regulon-target edges,",
      tt["regulon"].nunique() if len(tt) else 0, "regulons")
PY

echo "==============================================================="
echo " pySCENIC complete. Now run 05c_regulons_tcell_bcell.R"
echo "==============================================================="
