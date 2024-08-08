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

set.seed(100)
DATADIR = '~/neuroblastoma/data/'
WORKDIR = '~/neuroblastoma/'
PREPRDATADIR = '~/neuroblastoma/preprocessing_data/'

print('Read data')
seu_new <- qread(paste0(PREPRDATADIR, 'seu_list_mypipelineannot_v1.qs'))

print('Merge data')
seu_new <- merge(seu_new[[1]], seu_new[-1])

plan("multicore", workers = 20)
options(future.globals.maxSize = 5 * 1024^3)

print('Normalize data')
seu_new <- seu_new %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)

print('Integrate data')
seu_new <- seu_new %>% 
  IntegrateLayers(
    object = ., method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE) %>% 
  RunUMAP(reduction = "integrated.cca", dims = 1:20) 

seu_new <- JoinLayers(seu_new)

print('Save data')

qsave(seu_new, paste0(PREPRDATADIR, 'seu_integrated_mypipelineannot_v1.qs'))