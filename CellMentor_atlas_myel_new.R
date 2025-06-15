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

CreateBatchRobustCSFNMFobject <- function(ref_matrix,
                                          ref_celltype,
                                          ref_batch,
                                          data_matrix,
                                          n_features = 5000,
                                          norm = TRUE,
                                          scale = TRUE,
                                          scale_by = "cells",
                                          verbose = TRUE,
                                          num_cores = 1) {
  
  if(!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required for batch-robust gene selection")
  }
  
  # Create reporter function for verbose output
  report <- function(message) {
    if(verbose) cat("[", format(Sys.time(), "%H:%M:%S"), "]", message, "\n")
  }
  
  report("Starting batch-robust CellMentor object creation")
  
  # Check for common genes between matrices
  report("Checking for common genes")
  ref_genes <- rownames(ref_matrix)
  data_genes <- rownames(data_matrix)
  
  if(is.null(ref_genes) || is.null(data_genes)) {
    stop("Both matrices must have gene names as rownames")
  }
  
  common_genes <- intersect(ref_genes, data_genes)
  report(paste("Found", length(common_genes), "common genes between reference and query"))
  
  if(length(common_genes) < 1000) {
    warning("Very few common genes found. Check gene name formatting.")
  }
  
  # Filter matrices to common genes
  ref_matrix <- ref_matrix[common_genes,]
  data_matrix <- data_matrix[common_genes,]
  
  # Create Seurat object for batch-robust gene selection
  report("Creating Seurat object for batch-robust gene selection")
  seurat_obj <- Seurat::CreateSeuratObject(counts = ref_matrix)
  seurat_obj$batch <- ref_batch
  
  # Split by batch
  report("Splitting by batch for integration feature selection")
  batch_list <- Seurat::SplitObject(seurat_obj, split.by = "batch")
  
  # Process each batch
  report(paste("Processing", length(batch_list), "batches"))
  for(i in seq_along(batch_list)) {
    batch_list[[i]] <- Seurat::NormalizeData(batch_list[[i]], verbose = FALSE)
    batch_list[[i]] <- Seurat::FindVariableFeatures(
      batch_list[[i]],
      selection.method = "vst",
      nfeatures = min(2000, ncol(batch_list[[i]])/2),
      verbose = FALSE
    )
  }
  
  # Find integration features
  report(paste("Selecting", n_features, "batch-robust integration features"))
  integration_features <- Seurat::SelectIntegrationFeatures(
    batch_list,
    nfeatures = min(n_features, length(common_genes)),
    verbose = FALSE
  )
  
  report(paste("Selected", length(integration_features), "batch-robust features"))
  
  # Create standard CSFNMF object with batch-robust genes
  report("Creating CellMentor object with batch-robust genes")
  object <- CreateCSFNMFobject(
    ref_matrix = ref_matrix[integration_features,],
    ref_celltype = ref_celltype,
    data_matrix = data_matrix[integration_features,],
    norm = norm,
    most.variable = FALSE,  # We already selected variable genes
    scale = scale,
    scale_by = scale_by,
    gene_list = NULL,
    verbose = verbose,
    num_cores = num_cores
  )
  
  # Store batch information and selected features
  object@annotation$batch <- ref_batch
  
  # Add attributes with batch-robust information
  attr(object, "batch_robust_info") <- list(
    n_original_genes = length(common_genes),
    n_selected_genes = length(integration_features),
    batch_robust_genes = integration_features,
    n_batches = length(unique(ref_batch))
  )
  
  report("Batch-robust CellMentor object creation complete")
  return(object)
}

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

ref <- subset(ref, Study == 'Verhoeven2022')

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

print("CREATE CELLMENTOR OBJECT")
object = CreateBatchRobustCSFNMFobject(
  ref_matrix = as.matrix(ref[['RNA']]@counts), 
  ref_celltype = ref$SimpleCellType,
  ref_batch = ref$Study,
  data_matrix = as.matrix(query[['RNA']]@counts),
  n_features = 5000,
  verbose = TRUE
)
qsave(object, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_object_new_v2.qs'))

# v1 - verhoeven only
# v2 - stupid subset - across batch genes
# v3 - smart subset
# object <- qread(paste0(PREPRDATADIR, 'atlas_myel_CellMentor_object_v3.qs'))

# v2 - CellMentor(object, num_cores = 10, alpha_range = c(1), beta_range = c(1), 
# gamma_range = c(1), delta_range = c(1)) - stupid subset
# v3 - CellMentor(object, num_cores = 10, alpha_range = c(1, 10, 20), beta_range = c(1, 10, 20), 
# gamma_range = c(0.001, 0.01), delta_range = c(0.1, 0.5)) - stupid subset
# v4 - only verhoeven - 1 set
# v5 - smart subset - 1 set

rank_result <- SelectRank(
  train_matrix = object@matrices@ref,
  max_p_value = 0.01,
  numCores = 10
)
full_object <- methods::new("traincsfnmf")
full_object@matrices <- object@matrices
full_object@annotation <- object@annotation
best_model <- RunCSFNMF(
  train_object = full_object,
  k = rank_result,
  const.alpha = 15,
  const.beta = 15,
  const.gamma = 0.001,
  const.delta = 0.5,
  max.iter = 100,
  verbose = T,
  num_cores = 10,
  whole_object = TRUE
)

print('SAVE BEST METHOD')
qsave(best_model, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_bestmethod_new_v2.qs'))

print('DONE')