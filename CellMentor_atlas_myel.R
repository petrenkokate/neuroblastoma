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
# ref <- subset(ref, Sample == 'Verhoeven2022_NB13')

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
ref$Study <- as.character(ref$Study)
ref <- subset(ref, Study %in% c('Kildisiute2021_10X', 'Slyper2020_cell', 'Slyper2020_nucleus', 'Costa2022'))
query <- subset(atlas, Final_Annotation %in% c('DC', 'Macrophages', 'Monocytes', 'Neutrophils', 'pDCs'))

rm(atlas, atlas_ref)

print("CREATE CELLMENTOR OBJECT")
object = CreateCSFNMFobject(as.matrix(ref[['RNA']]@counts), ref$SimpleCellType,
                            query[['RNA']]@counts)
qsave(object, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_object_v7.qs'))

# v1 - verhoeven only
# v2 - stupid subset - across batch genes
# v3 - smart subset
# v6 - only verhoeven - 1 sample Verhoeven2022_NB13 - 1 set
# v7 - only datasets that are not in our atlas
object <- qread(paste0(PREPRDATADIR, 'atlas_myel_CellMentor_object_v7.qs'))

# v2 - CellMentor(object, num_cores = 10, alpha_range = c(1), beta_range = c(1), 
# gamma_range = c(1), delta_range = c(1)) - stupid subset
# v3 - CellMentor(object, num_cores = 10, alpha_range = c(1, 10, 20), beta_range = c(1, 10, 20), 
# gamma_range = c(0.001, 0.01), delta_range = c(0.1, 0.5)) - stupid subset
# v4 - only verhoeven - 1 set
# v5 - smart subset - 1 set
# v6 - only verhoeven - 1 sample - 1 set
# v7 - only datasets that are not in our atlas
optimal_params_meth <- CellMentor(object, num_cores = 10, 
                                  gamma_range = c(0), delta_range = c(0))
print('SAVE BEST METHOD')
qsave(optimal_params_meth, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_bestmethod_7.qs'))

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