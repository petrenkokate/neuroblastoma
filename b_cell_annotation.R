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


set.seed(100)
DATADIR = '~/neuroblastoma/data/'
WORKDIR = '~/neuroblastoma/'
PREPRDATADIR = '~/neuroblastoma/preprocessing_data/'

my_colors <- c("#9ACD32FF", "#4682B4FF",  "#DDA0DDFF", "#FFA07AFF", "#8FBC8BFF",
               "#40E0D0FF",  "#F0E68CFF", "#5F9EA0FF", "#D2B48CFF",  
               "#6495EDFF", "#FF69B4FF", "#BA55D3FF",  "#F08080FF", "#32CD32FF",  
               "#FFDAB9FF",  "#87CEEBFF")



bped <- celldex::BlueprintEncodeData()
hum_atlas_data <- celldex::HumanPrimaryCellAtlasData()
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes 
data(EnsemblGeneTable.Hs)
scGate_models_DB <- get_scGateDB()

print('Integration')

seu <- qread( paste0(PREPRDATADIR, 'seu_integrated_mypipelineannot_v1.qs'))

seu <- seu %>% 
  RunHarmony(., 'SampleID',
             lambda = 1, verbose = FALSE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose=F) 

seu$ANNOTATION_FINAL <- ifelse(seu$ANNOTATION_FINAL == "endothelial_cells", "stromal_cells", seu$ANNOTATION_FINAL)
seu <- subset(seu, ANNOTATION_FINAL %in% c("myeloid_cells", "B_cells", "B_cells",
                                           "malignanB_cells","stromal_cells"))

qsave(seu, paste0(PREPRDATADIR, 'seu_mypipelineannot_harmony.qs'))

b_cells_h <- subset(seu, ANNOTATION_FINAL == 'B_cells')

print('B_cells annotation')
b_cells_list <- lapply(b_cells_h$SampleID %>% unique, function(sample) {
  
  print('B_cells clustering')
  b_cell_sample <- sample
  try({
    b_cell_sample <- subset(b_cells_h, SampleID == sample) %>% 
      NormalizeData() %>%
      FindVariableFeatures()%>%
      ScaleData() %>%
      RunPCA() %>% 
      RunUMAP(dims = 1:20) %>%
      FindNeighbors(dims = 1:20) %>%
      FindClusters(resolution = 1.5)
    
    print('B_cells sce')
    sce <- as.SingleCellExperiment(b_cell_sample)
    colLabels(sce) <- b_cell_sample$RNA_snn_res.1.5
    Idents(b_cell_sample) <- 'RNA_snn_res.1.5'
    
    print('B_cells singler')
    # #cluster-wise annotation
    pred_bped_main <- SingleR(test = sce, ref = bped, clusters=colLabels(sce),
                              labels = bped$label.main, BPPARAM=MulticoreParam(15))
    # # Create a named vector for RenameIdents
    singleR_labels <- pred_bped_main$pruned.labels
    names(singleR_labels) <- rownames(pred_bped_main)
    
    # Update cluster identities in Seurat object based on SingleR results
    b_cell_sample <- RenameIdents(b_cell_sample, singleR_labels)
    b_cell_sample[['SingleR_main']] <- Idents(b_cell_sample)
    print('B_cells singler2')
    Idents(b_cell_sample) <- 'RNA_snn_res.1.5'
    # #cluster-wise annotation
    pred_bped_fine <- SingleR(test = sce, ref = bped, clusters=colLabels(sce),
                              labels = bped$label.fine, BPPARAM=MulticoreParam(15))
    # # Create a named vector for RenameIdents
    singleR_labels <- pred_bped_fine$pruned.labels
    names(singleR_labels) <- rownames(pred_bped_fine)
    
    # Update cluster identities in Seurat object based on SingleR results
    b_cell_sample <- RenameIdents(b_cell_sample, singleR_labels)
    b_cell_sample[['SingleR_fine']] <- Idents(b_cell_sample)
    
    # scGate
    print('B_cells scgate')
    b_cell_sample <- scGate(b_cell_sample, model = scGate_models_DB$human$TME_HiRes, ncores = 15)
  })
  
  return(b_cell_sample)
})

names(b_cells_list) <- b_cells_h$SampleID %>% unique
qsave(b_cells_list, paste0(PREPRDATADIR, 'b_cells_human_list_annotation.qs'))