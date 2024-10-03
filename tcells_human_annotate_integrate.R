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
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db")
library(ReactomePA)
library(msigdbr)
library(qs)

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

print('reading data')
seu <- qread(paste0(PREPRDATADIR, 'seu_naivetreatment.qs'))

t_cells_h <- subset(seu, ANNOTATION_FINAL == 'T_cells')

t_cells_list <- lapply(t_cells_h$SampleID %>% unique, function(sample) {
  
  try({
  print(sample)
  if (sample %in% c('M228AAA_T1_wienke','Tumor_27_dong','M277AAB_T_wienke')) {
    t_cell_sample <- subset(t_cells_h, SampleID == sample) %>% 
      NormalizeData() %>%
      FindVariableFeatures()%>%
      ScaleData() %>%
      RunPCA(npcs = 20) %>% 
      RunUMAP(dims = 1:5) %>%
      FindNeighbors(dims = 1:5) %>%
      FindClusters(resolution = 1.5)
  }
  else {
    t_cell_sample <- subset(t_cells_h, SampleID == sample) %>% 
      NormalizeData() %>%
      FindVariableFeatures()%>%
      ScaleData() %>%
      RunPCA() %>% 
      RunUMAP(dims = 1:20) %>%
      FindNeighbors(dims = 1:20) %>%
      FindClusters(resolution = 1.5)
  }
  sce <- as.SingleCellExperiment(t_cell_sample)
  colLabels(sce) <- t_cell_sample$RNA_snn_res.1.5
  Idents(t_cell_sample) <- 'RNA_snn_res.1.5'
  
  # #cluster-wise annotation
  pred_bped_main <- SingleR(test = sce, ref = bped, clusters=colLabels(sce),
                            labels = bped$label.main, BPPARAM=MulticoreParam(15))
  # # Create a named vector for RenameIdents
  singleR_labels <- pred_bped_main$pruned.labels
  names(singleR_labels) <- rownames(pred_bped_main)
  
  # Update cluster identities in Seurat object based on SingleR results
  t_cell_sample <- RenameIdents(t_cell_sample, singleR_labels)
  t_cell_sample[['SingleR_main']] <- Idents(t_cell_sample)
  
  Idents(t_cell_sample) <- 'RNA_snn_res.1.5'
  # #cluster-wise annotation
  pred_bped_fine <- SingleR(test = sce, ref = bped, clusters=colLabels(sce),
                            labels = bped$label.fine, BPPARAM=MulticoreParam(15))
  # # Create a named vector for RenameIdents
  singleR_labels <- pred_bped_fine$pruned.labels
  names(singleR_labels) <- rownames(pred_bped_fine)
  
  # Update cluster identities in Seurat object based on SingleR results
  t_cell_sample <- RenameIdents(t_cell_sample, singleR_labels)
  t_cell_sample[['SingleR_fine']] <- Idents(t_cell_sample)
  
  # scGate
  message('scGate is going...')
  t_cell_sample <- scGate(t_cell_sample, model = scGate_models_DB$human$TME_HiRes, ncores = 15)
  
  return(t_cell_sample)
  })
})

print('merging')
t_cells_h <- merge(x=t_cells_list[[1]], y=t_cells_list[2:length(t_cells_list)])

print('integration')
t_cells_h <- t_cells_h %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose=FALSE) %>%
  RunHarmony(., 'SampleID',
             lambda = 1, verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, verbose=F) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

print('join layers')
t_cells_h <- JoinLayers(t_cells_h)

print('saving')
qsave(t_cells_h, paste0(PREPRDATADIR, 't_cells_human_annotation_merged_v2.qs'))
