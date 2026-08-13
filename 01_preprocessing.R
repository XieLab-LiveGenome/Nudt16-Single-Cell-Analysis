# =============================================================================
# 01_preprocessing.R  —  NUDT16 KO female spleen scRNA-seq
#   load -> QC -> doublets -> MAD filter -> SCT -> Harmony -> cluster ->
#   UMAP -> SingleR (secondary) -> composition -> save processed RDS
#
# INPUT  : PIPseeker filtered matrices (relative to project_dir; config setwd's).
# OUTPUT : female_obj_processed.rds   (RDS_PROCESSED)
#          qc/, umap/, composition/ figures + tables
#
# =============================================================================

if (!isTRUE(getOption("nudt16.config.loaded"))) {
  cand <- c(
    "00_config.R",
    file.path("pipeline", "00_config.R"),
    file.path("results", "females", "pipeline", "00_config.R"),
    "/Users/budankm/Desktop/Sequencing/GONG/results/females/pipeline/00_config.R"
  )
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Cannot locate 00_config.R — run from the pipeline folder.")
  source(hit[1])
}

message("\n=========  STAGE 01 — PREPROCESSING  =========")

# Sample manifest (female only — 2 WT, 2 KO). PIPseeker filtered matrices.
sample_dirs <- list(
  Nudt16Spleen_WTF1 = "Z_pipseeker_output/Nudt16Spleen_WTF1/filtered_matrix/sensitivity_4",
  Nudt16Spleen_KOF1 = "Z_pipseeker_output/Nudt16Spleen_KOF1/filtered_matrix/sensitivity_3",
  Nudt16Spleen_WTF2 = "Z_pipseeker_output/Nudt16Spleen_WTF2/filtered_matrix/sensitivity_2",
  Nudt16Spleen_KOF2 = "Z_pipseeker_output/Nudt16Spleen_KOF2/filtered_matrix/sensitivity_3"
)

# ---- 1.1 Load and merge -----------------------------------------------------
message("\n[1.1] Loading PIPseeker matrices...")

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

# ---- 1.2 QC metrics ---------------------------------------------------------
message("\n[1.2] Computing QC metrics...")

seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rp[ls]")

VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        group.by = "sample", ncol = 4, pt.size = 0.1)
ggsave(rf("qc", "QC_violins_prefilter.pdf"), width = 12, height = 4)

# ---- 1.3 Doublet removal (scDblFinder, per sample) --------------------------
message("\n[1.3] Running scDblFinder per sample...")

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

sce_dbl <- as.SingleCellExperiment(seu, assay = "RNA")
sce_dbl <- scDblFinder(sce_dbl, samples = sce_dbl$sample,
                       BPPARAM = SerialParam(RNGseed = DOUBLET_SEED))

seu$scDblFinder_class <- sce_dbl$scDblFinder.class
seu$scDblFinder_score <- sce_dbl$scDblFinder.score

cat("Doublet rates per sample:\n")
print(round(prop.table(table(seu$sample, seu$scDblFinder_class), 1), 3))

seu <- subset(seu, subset = scDblFinder_class == "singlet")
cat("Cells after doublet removal:", ncol(seu), "\n")

# ---- 1.4 MAD-based QC (nmads = 3, per sample) -------------------------------
message("\n[1.4] MAD-based QC filtering...")

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
ggsave(rf("qc", "QC_violins_postfilter.pdf"), width = 12, height = 4)

# ---- 1.5 SCTransform (per sample) -------------------------------------------
message("\n[1.5] SCTransform per sample...")

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

# ---- 1.6 PCA + Harmony + UMAP -----------------------------------------------
message("\n[1.6] PCA, Harmony, clustering, UMAP...")

female_obj <- RunPCA(female_obj, npcs = N_PCS, verbose = FALSE)
ggsave(rf("umap", "elbow.pdf"),
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

p1 <- DimPlot(female_obj, group.by = "sample")          + ggtitle("By sample")
p2 <- DimPlot(female_obj, group.by = "genotype")        + ggtitle("By genotype")
p3 <- DimPlot(female_obj, group.by = "replicate")       + ggtitle("By replicate")
p4 <- DimPlot(female_obj, group.by = "seurat_clusters",
              label = TRUE)                              + ggtitle("Clusters")

ggsave(rf("umap", "UMAP_sample.pdf"),    p1, width = 8, height = 6)
ggsave(rf("umap", "UMAP_genotype.pdf"),  p2, width = 8, height = 6)
ggsave(rf("umap", "UMAP_replicate.pdf"), p3, width = 8, height = 6)
ggsave(rf("umap", "UMAP_clusters.pdf"),  p4, width = 8, height = 6)
ggsave(rf("umap", "UMAP_all.pdf"), (p1 | p2) / (p3 | p4),
       width = 14, height = 12)

# ---- 1.7 SingleR annotation (SECONDARY / validation) ------------------------
message("\n[1.7] SingleR annotation against ImmGen (secondary label)...")

DefaultAssay(female_obj) <- "RNA"
female_obj <- JoinLayers(female_obj, assay = "RNA")
# Ensure an RNA 'data' layer exists for all downstream figure panels.
female_obj <- NormalizeData(female_obj, assay = "RNA", verbose = FALSE)

cts <- GetAssayData(female_obj, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(assays  = list(counts = cts),
                            colData = female_obj@meta.data)
logcounts(sce) <- log1p(counts(sce))

immgen <- celldex::ImmGenData()
sr <- SingleR(test = sce, ref = immgen,
              labels = immgen$label.main,
              BPPARAM = MulticoreParam(workers = 4))

female_obj$ImmGen.labels <- sr$labels

lab_tab     <- table(female_obj$ImmGen.labels)
keep_labels <- names(lab_tab)[lab_tab >= MIN_CELLS_PER_LABEL]
female_obj  <- subset(female_obj, subset = ImmGen.labels %in% keep_labels)

cat("\nCells per ImmGen label (post-filter):\n")
print(sort(table(female_obj$ImmGen.labels), decreasing = TRUE))

p_ct <- DimPlot(female_obj, group.by = "ImmGen.labels",
                label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("Females — SingleR/ImmGen (secondary)")
ggsave(rf("umap", "UMAP_ImmGen_secondary.pdf"), p_ct, width = 10, height = 8)

# ---- 1.8 Composition (SingleR secondary) ------------------------------------
message("\n[1.8] Composition analysis (secondary label)...")

celltype_genotype <- table(female_obj$ImmGen.labels, female_obj$genotype)
ct_counts <- as.data.frame.matrix(celltype_genotype) %>%
  tibble::rownames_to_column("CellType") %>%
  mutate(Total = KO + WT) %>% arrange(desc(Total))
write.csv(ct_counts, rf("composition", "ImmGen_by_Genotype_counts.csv"),
          row.names = FALSE)

ct_prop <- as.data.frame.matrix(prop.table(celltype_genotype, 2) * 100) %>%
  tibble::rownames_to_column("CellType")
write.csv(ct_prop, rf("composition", "ImmGen_by_Genotype_proportions.csv"),
          row.names = FALSE)

# ---- 1.9 Save processed object ----------------------------------------------
message("\n[1.9] Saving female_obj_processed.rds ...")
saveRDS(female_obj, RDS_PROCESSED)
cat("Saved:", RDS_PROCESSED, "\n")

message("\n=========  STAGE 01 complete  =========")
