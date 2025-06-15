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
atlas_ref <- readRDS(paste0(DATADIR, 'NB_atlas/NB_atlas_v2.rds'))
atlas <- qread(paste0(PREPRDATADIR, 'ATLAS_object.qs'))

sample_metadata <- read.csv(paste0(DATADIR, 'samples_metadata_v5.csv'))
sample_metadata[sample_metadata==""] <- NA

rownames(sample_metadata) <- sample_metadata$Sample_dataset



clean_sample_id <- function(id) {
  # First remove _rep1 and _rep2
  id <- gsub("_rep[12]$", "", id)

  # Handle the special bonine cases
  id <- gsub("_sn$", "", id)      # Remove _sn suffix
  id <- gsub("_[12]$", "", id)    # Remove _1 or _2 suffix

  return(id)
}

atlas$less18M <- sample_metadata[clean_sample_id(atlas$SampleID),]$less18M
atlas$Sex <- sample_metadata[clean_sample_id(atlas$SampleID),]$Sex
atlas$MYCN_status <- sample_metadata[clean_sample_id(atlas$SampleID),]$MYCN
atlas$method <- sample_metadata[clean_sample_id(atlas$SampleID),]$Method

metadata_python <- read.csv("query_metadata_with_umap.csv", row.names = 1)

# Match to Seurat cells
metadata_python <- metadata_python[colnames(atlas), ]
atlas <- AddMetaData(atlas, metadata = metadata_python[c("Age", "celltypist_cell_label_fine", "Final_Annotation")])

# Add UMAP
atlas[['UMAP']] <- CreateDimReducObject(embeddings = as.matrix(metadata_python[, c("UMAP_1", "UMAP_2")]), key = "UMAP_", global = T, assay = "RNA")

ref <- subset(atlas_ref, Cell_type_wImmuneZoomAnnot %in% c('cDC1', 'cDC2/DC3', 'Classical monocyte', 'Macrophage',
                                                           'Migratory cDC', 'Neutrophil', 'Patrolling monocyte', 'pDC'))

# ref <- subset(ref, Study == 'Verhoeven2022')

ref@meta.data <- ref@meta.data %>%
  mutate(SimpleCellType = case_when(
    Cell_type_wImmuneZoomAnnot %in% c('cDC1', 'cDC2/DC3', 'Migratory cDC') ~ 'cDC',
    Cell_type_wImmuneZoomAnnot == 'pDC' ~ 'pDC',
    Cell_type_wImmuneZoomAnnot == 'Macrophage' ~ 'Macrophages',
    Cell_type_wImmuneZoomAnnot %in% c('Classical monocyte') ~ 'Classical monocytes',
    Cell_type_wImmuneZoomAnnot %in% c('Patrolling monocyte') ~ 'Patrolling monocytes',
    Cell_type_wImmuneZoomAnnot == 'Neutrophil' ~ 'Neutrophils',
    TRUE ~ 'Other'
  ))
ref <- subset(ref, SimpleCellType != 'Other')
query <- subset(atlas, Final_Annotation %in% c('DC', 'Macrophages', 'Monocytes', 'Neutrophils', 'pDCs'))

rm(atlas, atlas_ref)

#' Perform batch correction on reference and query datasets separately
#' 
#' @param ref_matrix Reference expression matrix
#' @param ref_metadata Reference metadata with batch information
#' @param query_matrix Query expression matrix
#' @param query_metadata Query metadata with batch information
#' @param var_features Number of variable features to select
#' @return List with batch-corrected reference and query matrices
correct_batches_separately <- function(ref_matrix, 
                                       ref_metadata, 
                                       query_matrix, 
                                       query_metadata,
                                       var_features = 3000) {
  # 1. Correct reference data
  ref_seurat <- CreateSeuratObject(counts = ref_matrix)
  ref_seurat$study <- ref_metadata$study_id
  ref_seurat$patient <- ref_metadata$patient_id
  ref_seurat$protocol <- ref_metadata$protocol
  
  # Process reference data
  ref_seurat <- NormalizeData(ref_seurat)
  ref_seurat <- FindVariableFeatures(ref_seurat, nfeatures = var_features)
  ref_seurat <- ScaleData(ref_seurat)
  ref_seurat <- RunPCA(ref_seurat, npcs = 50)
  
  # Apply Harmony to reference data (study level)
  ref_seurat <- RunHarmony(ref_seurat, "study")
  
  # 2. Correct query data separately using same workflow
  query_seurat <- CreateSeuratObject(counts = query_matrix)
  query_seurat$study <- query_metadata$study_id
  query_seurat$patient <- query_metadata$patient_id
  query_seurat$protocol <- query_metadata$protocol
  
  # Process query data
  query_seurat <- NormalizeData(query_seurat)
  query_seurat <- FindVariableFeatures(query_seurat, nfeatures = var_features)
  query_seurat <- ScaleData(query_seurat, features = rownames(ref_seurat))
  query_seurat <- RunPCA(query_seurat, npcs = 50, features = VariableFeatures(ref_seurat))
  
  # Apply Harmony to query data (study level)
  query_seurat <- RunHarmony(query_seurat, "study")
  
  # 3. Get corrected matrices from both objects
  # Get common features
  common_features <- intersect(rownames(ref_matrix), rownames(query_matrix))
  
  # Extract corrected expression
  ref_corrected <- back_project_pca(
    embeddings = Embeddings(ref_seurat, "harmony"),
    loadings = ref_seurat[["pca"]]@feature.loadings,
    features = common_features
  )
  
  query_corrected <- back_project_pca(
    embeddings = Embeddings(query_seurat, "harmony"),
    loadings = query_seurat[["pca"]]@feature.loadings,
    features = common_features
  )
  
  return(list(
    ref_matrix = ref_corrected,
    query_matrix = query_corrected,
    common_features = common_features
  ))
}

#' Back-project PCA/Harmony embeddings to expression space
#' 
#' @param embeddings Cell embeddings (cells × dimensions)
#' @param loadings PCA loadings (genes × dimensions)
#' @param features Features to include in output
#' @return Expression matrix in original gene space
back_project_pca <- function(embeddings, loadings, features) {
  # Make sure dimensions match
  loadings <- loadings[, 1:ncol(embeddings)]
  
  # Back project to gene expression space
  reconstructed <- t(embeddings %*% t(loadings))
  
  # Subset to requested features
  reconstructed <- reconstructed[features, ]
  
  return(as(reconstructed, "CsparseMatrix"))
}


print("CREATE CELLMENTOR OBJECT")
object = CreateCSFNMFobject(as.matrix(ref[['RNA']]@counts), ref$SimpleCellType,
                            query[['RNA']]@counts)
qsave(object, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_object_v6.qs'))


# object <- qread(paste0(PREPRDATADIR, 'atlas_myel_CellMentor_object_v6.qs'))


optimal_params_meth <- CellMentor(object, num_cores = 10, alpha_range = c(1), beta_range = c(1), 
                                  gamma_range = c(1), delta_range = c(1))
print('SAVE BEST METHOD')
qsave(optimal_params_meth, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_bestmethod_v6.qs'))

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
qsave(h_test, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_bestmethod_v5_htest.qs'))

print('DONE')