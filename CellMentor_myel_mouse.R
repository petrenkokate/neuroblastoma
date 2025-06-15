library(data.table)
library(Seurat)
library(scDblFinder)
library(DoubletFinder)
library(dplyr)
library(remotes)
# remotes::install_github("carmonalab/STACAS")
library(STACAS)
library(SingleR)
# install.packages("scGate")
library(scGate)
library(ggplot2)
library(BiocParallel)
library(harmony)
library(RColorBrewer)
library(stringr)
library(cowplot)
library(scuttle)
# install.packages("devtools")
# BiocManager::install("Homo.sapiens")
# BiocManager::install("jpmam1/scalop") #forked copy of a really good guy 
# devtools::install_github("jlaffy/infercna")
library(infercna)
library(infercnv)
library(parallel)
library(future)
library(clustree)
library(dittoSeq)
# install.packages("ggh4x")
library(ggh4x)
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db")
library(ReactomePA)
library(msigdbr)
# install.packages('qs')
library(qs)
library(DESeq2)
library(magrittr)
library(tibble)
library(tidyr)
library(pheatmap)
# BiocManager::install('apeglm')
# install.packages('ashr')

library(reticulate)
# install.packages('anndata')
library(anndata)
# devtools::install_github("cellgeni/sceasy")
# BiocManager::install(c("LoomExperiment"))
library(sceasy)
library(SeuratData)
library(CellMentor)
library(cluster)
library(aricode)
library(RcppML)
library(spam)



set.seed(100)
DATADIR = '~/neuroblastoma/data/'
WORKDIR = '~/neuroblastoma/'
PREPRDATADIR = '~/neuroblastoma/preprocessing_data/'

print("READ DATA")
seu <- qread(paste0(PREPRDATADIR, 'seu_mouse_neutrophils.qs'))
myeloid <- subset(seu, annotation == 'Myeloid_cells')

mouse_ref <- qread('mouse_ref/mouse_ref_myelsub.qs')
# mouse_ref_imm <- celldex::ImmGenData()
# logcounts_mat <- assay(mouse_ref_imm, "logcounts")
# 
# # Approximate raw counts
# approx_raw_counts <- 2^(logcounts_mat) - 1
# 
# # Create Seurat object
# mouse_ref <- CreateSeuratObject(
#   counts = logcounts_mat,
#   meta.data = as.data.frame(colData(mouse_ref_imm))
# )
# mouse_ref <- subset(mouse_ref, label.main %in% c("Macrophages", "Monocytes", "DC", "Neutrophils"))

print("CREATE CELLMENTOR OBJECT")
object = CreateCSFNMFobject(GetAssayData(subset(mouse_ref, Batch == 'Muscle_1'), slot = "counts"), subset(mouse_ref, Batch == 'Muscle_1')$MyeloidAnnotation,
                            GetAssayData(myeloid, slot = "counts"))
qsave(object, paste0(PREPRDATADIR, 'mouse_myel_CellMentor_object_v4.qs'))

#v1 - approx count
#v2 - tabula maris ref
#v3 - object = CreateCSFNMFobject(GetAssayData(subset(mouse_ref, Batch == 'PeripheralBlood_2')
#v4 - Muscle_1

# object <- qread(paste0(PREPRDATADIR, 'atlas_myel_CellMentor_object_v6.qs'))

optimal_params_meth <- CellMentor(object, num_cores = 10, gamma_range = c(0), delta_range = c(0))
print('SAVE BEST METHOD')
qsave(optimal_params_meth, paste0(PREPRDATADIR, 'mouse_myel_CellMentor_bestmethod_v4.qs'))

print('PROJECTION')
final_model <- optimal_params_meth$best_model
h_test <- project_data(
  W = final_model@W,                    # Use learned W matrix
  X = final_model@matrices@data,        # Query/test data matrix
  seed = 1,
  num_cores = 10,
  chunk_size = 10000,
  verbose = TRUE
)

print('SAVE H TEST')
qsave(h_test, paste0(PREPRDATADIR, 'mouse_myel_CellMentor_bestmethod_v4_htest.qs'))

print('DONE')

optimal_params_meth <- qread(paste0(PREPRDATADIR, 'mouse_myel_CellMentor_bestmethod_v4.qs'))
final_model <- optimal_params_meth$best_model

mouse_ref$CellMentor <- CreateDimReducObject(
  embeddings = t(as.matrix(final_model@H)),
  key = paste0('CellMentor', "_"),
  assay = DefaultAssay(mouse_ref)
)

mouse_ref <- mouse_ref %>%
  RunUMAP(reduction = 'CellMentor', dims= 1:optimal_params_meth$best_params$k, reduction.name = 'umap_cellmentor', verbose = F)

DimPlot(mouse_ref, reduction = 'umap_cellmentor', group.by = 'MyeloidAnnotation', cols = my_colors)
DimPlot(mouse_ref, reduction = 'umap_cellmentor', group.by = 'Batch') + NoLegend()


h_test <- qread(paste0(PREPRDATADIR, 'mouse_myel_CellMentor_bestmethod_v4_htest.qs'))
myeloid$CellMentor <- CreateDimReducObject(
  embeddings = t(as.matrix(h_test)),
  key = paste0('CellMentor', "_"),
  assay = DefaultAssay(myeloid),
  loadings = as.matrix(final_model@W)
)

myeloid <- myeloid %>%
  RunUMAP(reduction = 'CellMentor', dims= 1:optimal_params_meth$best_params$k, reduction.name = 'umap_cellmentor', verbose = F) %>% 
  FindNeighbors(dims = 1:5, reduction = 'CellMentor') %>%
  FindClusters(resolution = 0.3)

sce <- as.SingleCellExperiment(myeloid)

Idents(myeloid) <- 'RNA_snn_res.0.3'
markers <- FindAllMarkers(myeloid, only.pos = T, min.pct = 0.25)
top5 <- markers %>%
  dplyr::filter(!grepl("^Gm[0-9]", rownames(.))) %>% 
  dplyr::filter(!grepl("Rik$", rownames(.))) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC)

DotPlot(myeloid, dc_markers, cols = c('lightgrey', "#BA55D3FF")) +
  RotatedAxis()
# Monocytes
monocyte_markers <- c("Ly6c2", "Plac8", "Lyz2", "Ccl8", "Chil3", "Ifitm6")

# Macrophages
macrophage_markers <- c("F13a1", "Fcgr1", "Trem2", "Cd68", "Vcan", "Apoe", "C1qa", "Trem1")

# Neutrophils
neutrophil_markers <- c("S100a8", "S100a9", "Lcn2", "Chil1", "Cxcr2", "Retnlg")

# Dendritic Cells (DCs)
dc_markers <- c("H2-DMb1", "H2-Eb1", "H2-Aa", "Cd74", "Fscn1")

DotPlot(myeloid, top5$gene, cols = c('lightgrey', "#BA55D3FF")) +
  RotatedAxis()

Idents(object = myeloid) <- "RNA_snn_res.0.3"

new.cluster.ids <- c('0'='Macrophages',
                     '1'='DC',
                     '2'='Monocytes',
                     '3'='Macrophages',
                     '4'='Neutrophils'
)

myeloid <- RenameIdents(myeloid, new.cluster.ids)
myeloid$annotation <- Idents(myeloid)

DimPlot(myeloid, reduction = 'umap_cellmentor', group.by = 'annotation', cols = my_colors) + ggtitle('Annotation')

plot_grid(VlnPlot(myeloid, features = c('S100a8'), group.by = 'annotation', split.by = 'age', split.plot = T, cols = c("#FF69B4FF", "#6495EDFF")) + NoLegend(), 
          VlnPlot(myeloid, features = c('S100a9'), group.by = 'annotation', split.by = 'age', split.plot = T, cols = c("#FF69B4FF", "#6495EDFF")), nrow = 1, rel_widths = c(0.85, 1))

qsave(myeloid, paste0(PREPRDATADIR, 'mouse_myeloid_cellmentor.qs'))
