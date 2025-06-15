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


library(devtools)
setwd('~/CellMentor/')
load_all()

print("READ CELLMENTOR OBJECT")
object <- qread(paste0(PREPRDATADIR, 'atlas_CellMentor_object_v3.qs'))

optimal_params_meth <- select_optimal_parameters(object, num_cores = 10, alpha_range = c(1), beta_range = c(1),
                                                gamma_range = c(1), delta_range = c(1))
print('SAVE BEST METHOD')
qsave(optimal_params_meth, paste0(PREPRDATADIR, 'atlas_CellMentor_bestmethod_v3.qs'))

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
qsave(h_test, paste0(PREPRDATADIR, 'atlas_CellMentor_bestmethod_v3_htest.qs'))

print('DONE')