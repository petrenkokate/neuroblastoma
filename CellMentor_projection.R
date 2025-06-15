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

set.seed(100)
DATADIR = '~/neuroblastoma/data/'
WORKDIR = '~/neuroblastoma/'
PREPRDATADIR = '~/neuroblastoma/preprocessing_data/'

print("READ CELLMENTOR OBJECT")
optimal_params_meth <- qread(paste0(PREPRDATADIR, 'atlas_CellMentor_bestmethod_v2.qs'))

final_model <- optimal_params_meth$best_model
rm(optimal_params_meth)
h_test <- project_data(
  W = final_model@W,                    # Use learned W matrix
  X = final_model@matrices@data,        # Query/test data matrix
  seed = 1,
  num_cores = 5,
  chunk_size = 10000,
  verbose = TRUE
)

print('SAVE H TEST')
qsave(h_test, paste0(PREPRDATADIR, 'atlas_CellMentor_bestmethod_v2_htest.qs'))

print('DONE')