library(data.table)
library(Seurat)
library(scDblFinder)
library(DoubletFinder)
library(dplyr)
library(remotes)
library(STACAS)
library(SingleR)
library(scGate)
library(ggplot2)
library(BiocParallel)
library(harmony)
library(RColorBrewer)
library(stringr)
library(cowplot)
library(scuttle)
library(infercna)
library(infercnv)
library(parallel)
library(future)
library(clustree)
library(dittoSeq)
library(ggh4x)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Mm.eg.db")
library(ReactomePA)
library(msigdbr)
library(qs)

set.seed(100)
DATADIR = '~/neuroblastoma/data/'
WORKDIR = '~/neuroblastoma/'
PREPRDATADIR = '~/neuroblastoma/preprocessing_data/'

malig_cells <- qread(paste0(PREPRDATADIR, 'malig_cells_mouse_gsvas3.qs'))
# print('Read data')
# seu <- qread(paste0(PREPRDATADIR, 'seu_mouse.qs'))
# 
# malig_cells <- subset(seu, annotation == 'malignant_cells')
# 
# print('Malignant cells analysis - step 1')
# 
# malig_cells <- malig_cells %>% 
#   NormalizeData() %>% 
#   FindVariableFeatures() %>% 
#   ScaleData() %>% 
#   RunPCA(verbose=FALSE) %>% 
#   RunHarmony(., 'age',
#              lambda = 1, verbose = FALSE) %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20, verbose=F) %>% 
#   FindNeighbors(dims = 1:20) %>%
#   FindClusters(resolution = 0.1)
# 
# malig_cells <- subset(malig_cells, RNA_snn_res.0.1 != 2)
# 
# print('Malignant cells analysis - step 2')
# 
# malig_cells <- malig_cells %>% 
#   NormalizeData() %>% 
#   FindVariableFeatures() %>% 
#   ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
#   RunPCA(verbose=FALSE) %>% 
#   RunHarmony(., 'age',
#              lambda = 1, verbose = FALSE) %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20, verbose=F) %>% 
#   FindNeighbors(dims = 1:20) %>%
#   FindClusters(resolution = seq(0.1, 1, 0.1))
# 
# # 
# malig_cells <- subset(malig_cells, RNA_snn_res.0.1 != 3)
# 
# Idents(malig_cells) <- malig_cells$RNA_snn_res.0.1

print('Gene sets downloading')
library(escape)
GS.hallmark <- getGeneSets(species = 'Mus musculus', library = "H")
# GS.reactome <- getGeneSets(species = 'Mus musculus', library = "C2", subcategory = 'CP:REACTOME')
GS.tf <- getGeneSets(species = 'Mus musculus', library = "C3", subcategory = 'GTRD')

print('ssGSEA - hallmarks')
malig_cells <- escape::runEscape(malig_cells,
                                 method = "ssGSEA",
                                 gene.sets = GS.hallmark,
                                 groups = 500,
                                 min.size = 5,
                                 new.assay.name = "hallmark.ssGSEA",
                                 BPPARAM = MulticoreParam(20))

print('GSVA - hallmarks - save')
qsave(malig_cells, paste0(PREPRDATADIR, 'malig_cells_mouse_gsvas1.qs'))

# print('ssGSEA - reactome')
# malig_cells <- escape::runEscape(malig_cells, 
#                                  method = "ssGSEA",
#                                  gene.sets = GS.reactome,
#                                  groups = 500,
#                                  min.size = 5,
#                                  new.assay.name = "reactome.ssGSEA",
#                                  BPPARAM = MulticoreParam(20))
# print('ssGSEA - reactome - save')
# qsave(malig_cells, paste0(PREPRDATADIR, 'malig_cells_mouse_gsvas2.qs'))

print('ssGSEA - TF')
malig_cells <- escape::runEscape(malig_cells, 
                                 method = "ssGSEA",
                                 gene.sets = GS.tf, 
                                 groups = 500, 
                                 min.size = 5,
                                 new.assay.name = "tf.ssGSEA",
                                 BPPARAM = MulticoreParam(20))
print('ssGSEA - TF - save')
qsave(malig_cells, paste0(PREPRDATADIR, 'malig_cells_mouse_gsvas3_ssGSEA.qs'))