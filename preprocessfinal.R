# preprocess final

# =============================================================================
# MASTER SCRIPT — NUDT16 KO Female Spleen scRNA-seq (F1 + F2, 2 WT + 2 KO)
#
# End-to-end pipeline:
#   1.  Load PIPseeker filtered matrices
#   2.  QC metrics
#   3.  Per-sample doublet removal (scDblFinder)
#   4.  MAD-based per-sample QC filtering
#   5.  SCTransform normalization per sample
#   6.  Harmony integration, clustering, UMAP
#   7.  Cell type annotation (SingleR + ImmGen)
#   8.  Cell type composition tables / bar plot
#   9.  Pseudobulk DE (DESeq2) — global + per cell type
#   10. Volcano plots
#   11. Pathway enrichment — GO BP/MF/CC + KEGG + MSigDB Hallmark
#   12. Per-gene plots — violin (by genotype, by cell type) + FeaturePlot
#   13. Heatmaps — pseudobulk Z-scored, samples × genes
#   14. Excel export of DE + enrichment tables
#   15. Save objects + sessionInfo
# =============================================================================

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  # Core
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)

  # Normalization & integration
  library(glmGamPoi)
  library(harmony)

  # QC
  library(scDblFinder)
  library(scater)
  library(SingleCellExperiment)
  library(BiocParallel)

  # Annotation
  library(SingleR)
  library(celldex)

  # Differential expression
  library(DESeq2)

  # Enrichment
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(msigdbr)

  # Visualization
  library(pheatmap)
  library(RColorBrewer)

  # Export
  library(openxlsx)
})
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflicts_prefer(
  base::rownames, base::colnames,
  base::intersect, base::setdiff, base::union,
  base::which, base::is.element, base::table,
  base::which.max, base::which.min,
  .quiet = TRUE
)

set.seed(217)
options(Seurat.object.assay.version = "v5")

# ---- Paths ------------------------------------------------------------------
project_dir <- "/Users/budankm/Desktop/Sequencing/GONG"
out_dir     <- file.path(project_dir, "SCRNASEQ", "Preprocessing")
for (sub in c("qc", "umap", "composition", "volcano", "enrichment",
              "gene_plots", "heatmap")) {
  dir.create(file.path(out_dir, sub), recursive = TRUE, showWarnings = FALSE)
}
setwd(project_dir)

# ---- Parameters (single source of truth) ------------------------------------
# QC / filtering
MAD_NMADS               <- 3
DOUBLET_SEED            <- 217
MIN_CELLS_PER_LABEL     <- 15      # ImmGen label retention threshold

# Embedding
N_VAR_FEATURES          <- 3000
N_PCS                   <- 50
N_DIMS_USE              <- 20
CLUSTER_RES             <- 0.5

# Pseudobulk DE
MIN_CELLS_PER_SAMPLE    <- 50
MIN_REPS_PER_GENOTYPE   <- 2
MIN_COUNT               <- 10
MIN_SAMPLES             <- 2
LFC_CUT                 <- 0.5
PADJ_CUT                <- 0.05

# Enrichment thresholds
ENR_PVAL                <- 0.05
ENR_QVAL                <- 0.25

# Color anchors
COL_UP   <- "#D7191C"     # red
COL_DOWN <- "#2C7FB8"     # blue
COL_NS   <- "grey80"
COL_HM   <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)

# Sample manifest (female only — 2 WT, 2 KO)
sample_dirs <- list(
  Nudt16Spleen_WTF1 = "Z_pipseeker_output/Nudt16Spleen_WTF1/filtered_matrix/sensitivity_4",
  Nudt16Spleen_KOF1 = "Z_pipseeker_output/Nudt16Spleen_KOF1/filtered_matrix/sensitivity_3",
  Nudt16Spleen_WTF2 = "Z_pipseeker_output/Nudt16Spleen_WTF2/filtered_matrix/sensitivity_2",
  Nudt16Spleen_KOF2 = "Z_pipseeker_output/Nudt16Spleen_KOF2/filtered_matrix/sensitivity_3"
)

# =============================================================================
# 1. LOAD AND MERGE
# =============================================================================
message("\n[1] Loading PIPseeker matrices...")

seurat_list <- lapply(names(sample_dirs), function(s) {
  mtx <- Read10X(data.dir = sample_dirs[[s]], gene.column = 2)
  obj <- CreateSeuratObject(counts = mtx, project = s,
                            min.cells = 3, min.features = 200)
  obj$sample     <- s
  obj$orig.ident <- s
  obj
})
names(seurat_list) <- names(sample_dirs)

seu <- merge(seurat_list[[1]], y = seurat_list[-1],
             add.cell.ids = names(seurat_list))

seu$genotype  <- factor(ifelse(grepl("KO", seu$sample), "KO", "WT"),
                        levels = c("WT", "KO"))
seu$replicate <- ifelse(grepl("F1", seu$sample), "F1", "F2")

cat("Cells loaded:", ncol(seu), "\n")

# =============================================================================
# 2. QC METRICS
# =============================================================================
message("\n[2] Computing QC metrics...")

seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rp[ls]")

VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        group.by = "sample", ncol = 4, pt.size = 0.1)
ggsave(file.path(out_dir, "qc", "QC_violins_prefilter.pdf"),
       width = 12, height = 4)

# =============================================================================
# 3. DOUBLET REMOVAL (scDblFinder, per sample)
# =============================================================================
message("\n[3] Running scDblFinder per sample...")

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

sce <- as.SingleCellExperiment(seu, assay = "RNA")
sce <- scDblFinder(sce, samples = sce$sample,
                   BPPARAM = SerialParam(RNGseed = DOUBLET_SEED))

seu$scDblFinder_class <- sce$scDblFinder.class
seu$scDblFinder_score <- sce$scDblFinder.score

cat("Doublet rates per sample:\n")
print(round(prop.table(table(seu$sample, seu$scDblFinder_class), 1), 3))

seu <- subset(seu, subset = scDblFinder_class == "singlet")
cat("Cells after doublet removal:", ncol(seu), "\n")

# =============================================================================
# 4. MAD-BASED QC (nmads = 3, per sample)
#    nFeature_RNA : both tails on log scale
#    nCount_RNA   : upper tail on log scale
#    percent.mt   : upper tail
#    percent.ribo : upper tail
# =============================================================================
message("\n[4] MAD-based QC filtering...")

discard <-
  isOutlier(seu$nFeature_RNA, nmads = MAD_NMADS, type = "both",   log = TRUE,  batch = seu$sample) |
  isOutlier(seu$nCount_RNA,   nmads = MAD_NMADS, type = "higher", log = TRUE,  batch = seu$sample) |
  isOutlier(seu$percent.mt,   nmads = MAD_NMADS, type = "higher", log = FALSE, batch = seu$sample) |
  isOutlier(seu$percent.ribo, nmads = MAD_NMADS, type = "higher", log = FALSE, batch = seu$sample)

cat("Discarded:", sum(discard), " | Retained:", sum(!discard), "\n")
seu <- seu[, !discard]

cat("\nCells per sample (post-QC):\n");   print(table(seu$sample))
cat("\nCells per genotype (post-QC):\n"); print(table(seu$genotype))

VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        group.by = "sample", ncol = 4, pt.size = 0.1)
ggsave(file.path(out_dir, "qc", "QC_violins_postfilter.pdf"),
       width = 12, height = 4)

# =============================================================================
# 5. SCTRANSFORM (per sample)
# =============================================================================
message("\n[5] SCTransform per sample...")

obj_list <- SplitObject(seu, split.by = "sample")
obj_list <- lapply(obj_list, SCTransform,
                   vars.to.regress = c("percent.mt", "percent.ribo"),
                   method     = "glmGamPoi",
                   vst.flavor = "v2",
                   verbose    = FALSE)

hvgs       <- SelectIntegrationFeatures(obj_list, nfeatures = N_VAR_FEATURES)
female_obj <- merge(obj_list[[1]], y = obj_list[-1], merge.data = TRUE)
VariableFeatures(female_obj) <- hvgs
DefaultAssay(female_obj)     <- "SCT"

# =============================================================================
# 6. PCA + HARMONY + UMAP
# =============================================================================
message("\n[6] PCA, Harmony, clustering, UMAP...")

female_obj <- RunPCA(female_obj, npcs = N_PCS, verbose = FALSE)
ggsave(file.path(out_dir, "umap", "elbow.pdf"),
       ElbowPlot(female_obj, ndims = N_PCS) + ggtitle("Females — PCA elbow"),
       width = 7, height = 5)

female_obj <- RunHarmony(female_obj,
                         group.by.vars  = "sample",
                         assay.use      = "SCT",
                         reduction      = "pca",
                         reduction.save = "harmony",
                         verbose        = FALSE)

female_obj <- FindNeighbors(female_obj, reduction = "harmony",
                            dims = 1:N_DIMS_USE, verbose = FALSE)
female_obj <- FindClusters(female_obj, resolution = CLUSTER_RES, verbose = FALSE)
female_obj <- RunUMAP(female_obj, reduction = "harmony",
                      dims = 1:N_DIMS_USE, verbose = FALSE)

female_obj <- PrepSCTFindMarkers(female_obj)

# QC UMAPs (sample / genotype / replicate / clusters)
p1 <- DimPlot(female_obj, group.by = "sample")          + ggtitle("By sample")
p2 <- DimPlot(female_obj, group.by = "genotype")        + ggtitle("By genotype")
p3 <- DimPlot(female_obj, group.by = "replicate")       + ggtitle("By replicate")
p4 <- DimPlot(female_obj, group.by = "seurat_clusters",
              label = TRUE)                              + ggtitle("Clusters")

ggsave(file.path(out_dir, "umap", "UMAP_sample.pdf"),    p1, width = 8, height = 6)
ggsave(file.path(out_dir, "umap", "UMAP_genotype.pdf"),  p2, width = 8, height = 6)
ggsave(file.path(out_dir, "umap", "UMAP_replicate.pdf"), p3, width = 8, height = 6)
ggsave(file.path(out_dir, "umap", "UMAP_clusters.pdf"),  p4, width = 8, height = 6)
ggsave(file.path(out_dir, "umap", "UMAP_all.pdf"),
       (p1 | p2) / (p3 | p4), width = 14, height = 12)

# =============================================================================
# 7. CELL TYPE ANNOTATION (SingleR + ImmGen)
# =============================================================================
message("\n[7] SingleR annotation against ImmGen...")

DefaultAssay(female_obj) <- "RNA"
female_obj <- JoinLayers(female_obj, assay = "RNA")

cts <- GetAssayData(female_obj, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(assays  = list(counts = cts),
                            colData = female_obj@meta.data)
logcounts(sce) <- log1p(counts(sce))

immgen <- celldex::ImmGenData()
sr <- SingleR(test = sce, ref = immgen,
              labels = immgen$label.main,
              BPPARAM = MulticoreParam(workers = 4))

female_obj$ImmGen.labels <- sr$labels

# Drop labels with too few cells
lab_tab     <- table(female_obj$ImmGen.labels)
keep_labels <- names(lab_tab)[lab_tab >= MIN_CELLS_PER_LABEL]
female_obj  <- subset(female_obj, subset = ImmGen.labels %in% keep_labels)

cat("\nCells per ImmGen label (post-filter):\n")
print(sort(table(female_obj$ImmGen.labels), decreasing = TRUE))

p_ct       <- DimPlot(female_obj, group.by = "ImmGen.labels",
                      label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("Females — cell types")
p_ct_split <- DimPlot(female_obj, group.by = "ImmGen.labels",
                      split.by = "genotype",
                      label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("Females — cell types by genotype")

ggsave(file.path(out_dir, "umap", "UMAP_celltypes.pdf"),
       p_ct,       width = 10, height = 8)
ggsave(file.path(out_dir, "umap", "UMAP_celltypes_by_genotype.pdf"),
       p_ct_split, width = 16, height = 6)

# =============================================================================
# 8. CELL TYPE COMPOSITION
# =============================================================================
message("\n[8] Composition analysis...")

celltype_genotype <- table(female_obj$ImmGen.labels, female_obj$genotype)
ct_counts <- as.data.frame.matrix(celltype_genotype) %>%
  tibble::rownames_to_column("CellType") %>%
  mutate(Total = KO + WT) %>% arrange(desc(Total))
write.csv(ct_counts,
          file.path(out_dir, "composition", "CellType_by_Genotype_counts.csv"),
          row.names = FALSE)

ct_prop <- as.data.frame.matrix(prop.table(celltype_genotype, 2) * 100) %>%
  tibble::rownames_to_column("CellType")
write.csv(ct_prop,
          file.path(out_dir, "composition", "CellType_by_Genotype_proportions.csv"),
          row.names = FALSE)

ct_long <- ct_prop %>%
  pivot_longer(cols = c(KO, WT), names_to = "Genotype", values_to = "Percentage")
ct_long$Genotype <- factor(ct_long$Genotype, levels = c("WT", "KO"))

p_comp <- ggplot(ct_long, aes(x = reorder(CellType, -Percentage),
                              y = Percentage, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(WT = COL_DOWN, KO = COL_UP)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Female cell type distribution by genotype",
       x = "Cell type", y = "Percentage")
ggsave(file.path(out_dir, "composition", "CellType_proportions_barplot.pdf"),
       p_comp, width = 12, height = 6)

