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

atlas_ref <- qread(paste0(PREPRDATADIR, 'atlas_ref_subset_v2.qs'))
atlas <- qread(paste0(PREPRDATADIR, 'ATLAS_object.qs'))

sce <- as.SingleCellExperiment(atlas)
rm(atlas)
print('SingleR')
pred_singler <- SingleR(
  test = sce, 
  ref = as.SingleCellExperiment(atlas_ref), 
  labels = atlas_ref$Cell_type, 
  BPPARAM = MulticoreParam(10)
)
print('SAVING')
qsave(pred_singler, paste0(PREPRDATADIR, 'ATLAS_singler.qs'))

print('DONE')