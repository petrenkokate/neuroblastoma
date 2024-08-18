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

t_cells <- qread(paste0(PREPRDATADIR, 't_cells_mouse_v2.qs'))


print('Gene sets downloading')
library(escape)
GS.hallmark <- getGeneSets(species = 'Mus musculus', library = "H")
GS.reactome <- getGeneSets(species = 'Mus musculus', library = "C2", subcategory = 'CP:REACTOME')
GS.tf <- getGeneSets(species = 'Mus musculus', library = "C3", subcategory = 'GTRD')

print('ssGSEA - hallmarks')
t_cells <- escape::runEscape(t_cells,
                                 method = "ssGSEA",
                                 gene.sets = GS.hallmark,
                                 groups = 500,
                                 min.size = 5,
                                 new.assay.name = "hallmark.ssGSEA",
                                 BPPARAM = MulticoreParam(20))

print('GSVA - hallmarks - save')
qsave(t_cells, paste0(PREPRDATADIR, 't_cells_mouse_ssGSEAs1.qs'))

print('ssGSEA - reactome')
t_cells <- escape::runEscape(t_cells,
                                 method = "ssGSEA",
                                 gene.sets = GS.reactome,
                                 groups = 500,
                                 min.size = 5,
                                 new.assay.name = "reactome.ssGSEA",
                                 BPPARAM = MulticoreParam(20))
print('ssGSEA - reactome - save')
qsave(t_cells, paste0(PREPRDATADIR, 't_cells_mouse_ssGSEAs2.qs'))

print('ssGSEA - TF')
t_cells <- escape::runEscape(t_cells, 
                                 method = "ssGSEA",
                                 gene.sets = GS.tf, 
                                 groups = 500, 
                                 min.size = 5,
                                 new.assay.name = "tf.ssGSEA",
                                 BPPARAM = MulticoreParam(20))
print('ssGSEA - TF - save')
qsave(t_cells, paste0(PREPRDATADIR, 't_cells_mouse_ssGSEAs3.qs'))