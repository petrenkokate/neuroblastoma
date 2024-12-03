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
# library(infercnv)
library(parallel)
library(future)
library(clustree)
library(dittoSeq)
# install.packages("ggh4x")
library(ggh4x)
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(clusterProfiler)
# BiocManager::install('org.Mm.eg.db')
library("org.Mm.eg.db")
library(ReactomePA)
library(msigdbr)
# install.packages('qs')
library(qs)
# devtools::install_github("ncborcherding/escape")
# devtools::install_github("rcastelo/GSVA")
library(escape)
# install.packages('gprofiler2')
library(gprofiler2)
# remotes::install_github("carmonalab/ProjecTILs")
library(ProjecTILs)
library(gridExtra)
library(scCustomize)
library(enrichR)
# devtools::install_github("jinworks/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

set.seed(100)
DATADIR = '~/neuroblastoma/data/'
WORKDIR = '~/neuroblastoma/'
PREPRDATADIR = '~/neuroblastoma/preprocessing_data/'

my_colors <- c("#9ACD32FF", "#4682B4FF",  "#DDA0DDFF", "#FFA07AFF", "#8FBC8BFF",
               "#40E0D0FF",  "#F0E68CFF", "#5F9EA0FF", "#D2B48CFF",  
               "#6495EDFF", "#FF69B4FF", "#BA55D3FF",  "#F08080FF", "#32CD32FF",  
               "#FFDAB9FF",  "#87CEEBFF")

print("READ DATA")

seu <- qread(paste0(PREPRDATADIR, 'seu_mouse.qs'))
myeloid <- qread(paste0(PREPRDATADIR, 'myeloid_cells_mouse.qs'))
t_cells <- qread(paste0(PREPRDATADIR, 't_cells_mouse_v3.qs'))

Idents(object = seu) <- "RNA_snn_res.0.1"

new.cluster.ids <- c('0'='Malignant cells',
                     '1'='Malignant cells',
                     '2'='T cells',
                     '3'='B cells',
                     '4'='Myeloid cells',
                     '5'='Fibroblasts',
                     '6'='Malignant cells',
                     '6'='Malignant cells',
                     '7'='Endothelial cells',
                     '8'='Fibroblasts'
)

seu <- RenameIdents(seu, new.cluster.ids)
seu$annotation <- Idents(seu)

seu$deep_annotation <- seu$annotation
seu$deep_annotation <- as.character(seu$deep_annotation)
myeloid$annotation <- as.character(myeloid$annotation)
seu$deep_annotation[Cells(myeloid)] <- myeloid$annotation
seu$deep_annotation[Cells(t_cells)] <- t_cells$SingleR_imm_ref_v4

print("YOUNG DATA PROCESSING")
# Young
seu$samples <- seu$age
cellchat_young <- createCellChat(object = subset(seu, age == 'young'), group.by = "annotation", assay = "RNA")

CellChatDB <- CellChatDB.mouse 

CellChatDB.use <- subsetDB(CellChatDB) 
cellchat_young@DB <- CellChatDB.use
print("YOUNG DATA PROCESSING - identifyOverExpressedGenes")
cellchat_young <- subsetData(cellchat_young)
future::plan("multicore", workers = 8) # do parallel
cellchat_young <- identifyOverExpressedGenes(cellchat_young)
cellchat_young <- identifyOverExpressedInteractions(cellchat_young)

print("YOUNG DATA PROCESSING - computeCommunProb")
cellchat_young <- computeCommunProb(cellchat_young, type = "triMean")
cellchat_young <- filterCommunication(cellchat_young, min.cells = 10)

df.net.young <- subsetCommunication(cellchat_young)
print("YOUNG DATA PROCESSING - aggregateNet")
cellchat_young <- computeCommunProbPathway(cellchat_young)
cellchat_young <- aggregateNet(cellchat_young)
cellchat_young <- netAnalysis_computeCentrality(cellchat_young)


print("ADULT DATA PROCESSING")
cellchat_adult <- createCellChat(object = subset(seu, age == 'adult'), group.by = "annotation", assay = "RNA")

cellchat_adult@DB <- CellChatDB.use
print("ADULT DATA PROCESSING - identifyOverExpressedGenes")
cellchat_adult <- subsetData(cellchat_adult)
future::plan("multicore", workers = 8) # do parallel
cellchat_adult <- identifyOverExpressedGenes(cellchat_adult)
cellchat_adult <- identifyOverExpressedInteractions(cellchat_adult)
print("ADULT DATA PROCESSING - computeCommunProb")
cellchat_adult <- computeCommunProb(cellchat_adult, type = "triMean")
cellchat_adult <- filterCommunication(cellchat_adult, min.cells = 10)
df.net.adult <- subsetCommunication(cellchat_adult)
print("ADULT DATA PROCESSING - aggregateNet")
cellchat_adult <- computeCommunProbPathway(cellchat_adult)
cellchat_adult <- aggregateNet(cellchat_adult)
cellchat_adult <- netAnalysis_computeCentrality(cellchat_adult)

object.list <- list(young = cellchat_young, adult = cellchat_adult)
print("DATA SAVING")
qsave(object.list, paste0(PREPRDATADIR, 'cellchat_list_mouse.qs'))
print("DONE")
