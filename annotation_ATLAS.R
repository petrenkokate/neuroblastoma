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


print("DATABASES")
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

print("FUNCTION")

annotate_neuroblastoma <- function(sample) {
  tryCatch({
    print(paste('RUN', sample, 'RUN'))
  # Label mapping dictionary between SingleR and final annotations
  singler_to_final <- c(
    "DC" = "myeloid_cells",
    "Smooth_muscle_cells" = "stromal_cells",
    "Epithelial_cells" = "stromal_cells",
    "B_cell" = "b_cells",
    "Neutrophils" = "myeloid_cells",
    "T_cells" = "t_cells",
    "Monocyte" = "myeloid_cells",
    "Erythroblast" = "other",
    "BM & Prog." = "other",
    "Endothelial_cells" = "endothelial_cells",
    "Gametocytes" = "other",
    "Neurons" = "malignant_cells",
    "Keratinocytes" = "other",
    "HSC_-G-CSF" = "myeloid_cells",
    "Macrophage" = "myeloid_cells",
    "NK_cell" = "t_cells",
    "Embryonic_stem_cells" = "other",
    "Tissue_stem_cells" = "stromal_cells",
    "Chondrocytes" = "other",
    "Osteoblasts" = "other",
    "BM" = "other",
    "Platelets" = "other",
    "Fibroblasts" = "stromal_cells",
    "iPS_cells" = "other",
    "Hepatocytes" = "other",
    "MSC" = "stromal_cells",
    "Neuroepithelial_cell" = "malignant_cells",
    "Astrocyte" = "malignant_cells",
    "HSC_CD34+" = "other",
    "CMP" = "myeloid_cells",
    "GMP" = "myeloid_cells",
    "MEP" = "other",
    "Myelocyte" = "myeloid_cells",
    "Pre-B_cell_CD34-" = "b_cells",
    "Pro-B_cell_CD34+" = "b_cells",
    "Pro-Myelocyte" = "myeloid_cells"
  )
  print(paste('PREPROCESSING', sample, 'PREPROCESSING'))
  # Pre-processing
  seu <- subset(seu, subset = SampleID == sample) %>% 
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1.5)
  print(paste('SINGLER', sample, 'SINGLER'))
  # SingleR annotation
  sce <- as.SingleCellExperiment(seu)
  colLabels(sce) <- seu$RNA_snn_res.1.5
  
  pred_singler <- SingleR(
    test = sce, 
    ref = hum_atlas_data, 
    clusters = colLabels(sce),
    labels = hum_atlas_data$label.main, 
    BPPARAM = MulticoreParam(10)
  )
  
  # Store SingleR labels
  
  # Store both raw and harmonized labels in metadata
  seu$singler_labels_raw <- plyr::mapvalues(
    seu$RNA_snn_res.1.5,
    from = rownames(pred_singler),
    to = pred_singler$labels)
  seu$singler_labels_raw[is.na(seu$singler_labels_raw)] <- "unknown"
  singler_labels <- singler_to_final[as.character(seu$singler_labels_raw)]
  names(singler_labels) <- colnames(seu)
  seu$singler_labels <- singler_labels
  seu$singler_labels[is.na(seu$singler_labels)] <- "unknown"
  print(paste('SCGATE', sample, 'SCGATE'))
  # Run scGate for major cell types
  run_scgate <- function(seu, model_name) {
    seu <- scGate(seu, model = model_name)
    return(seu$is.pure == "Pure")
  }
  
  # Score major cell types
  seu$immune_cells <- run_scgate(seu, scGate_models_DB$human$TME_broad$Immune)
  seu$stromal_cells <- run_scgate(seu, scGate_models_DB$human$generic$Stromal)
  seu$endothelial_cells <- run_scgate(seu, scGate_models_DB$human$generic$Endothelial)
  
  # Score immune subtypes
  seu$t_cells <- run_scgate(seu, scGate_models_DB$human$generic$Tcell) | 
    run_scgate(seu, scGate_models_DB$human$generic$NK)
  seu$b_cells <- run_scgate(seu, scGate_models_DB$human$generic$Bcell) | 
    run_scgate(seu, scGate_models_DB$human$generic$PlasmaCell)
  seu$myeloid_cells <- run_scgate(seu, scGate_models_DB$human$generic$Myeloid) | 
    run_scgate(seu, scGate_models_DB$human$generic$Monocyte)
  
  # Score malignant cells
  malignant_model <- gating_model(
    name = "Malignant",
    signature = c("MYCN", "MLLT11", "SLC6A2", "PHOX2B", "NXPH1", "SYT1")
  )
  seu$malignant_cells <- run_scgate(seu, malignant_model)
  print(paste('INTEGRATION', sample, 'INTEGRATION'))
  # Integrate annotations
  integrate_annotations <- function(seu) {
    seu$final_annotation <- "unknown"
    
    for (cluster in unique(seu$RNA_snn_res.1.5)) {
      cells_in_cluster <- WhichCells(seu, idents = cluster)
      cluster_data <- seu[, cells_in_cluster]
      
      # Get SingleR annotation
      singler_annot <- names(which.max(table(cluster_data$singler_labels)))
      
      # Calculate scGate scores for the cluster
      scgate_scores <- list(
        stromal_cells = mean(cluster_data$stromal_cells),
        endothelial_cells = mean(cluster_data$endothelial_cells),
        malignant_cells = mean(cluster_data$malignant_cells),
        t_cells = mean(cluster_data$t_cells),
        b_cells = mean(cluster_data$b_cells),
        myeloid_cells = mean(cluster_data$myeloid_cells)
      )
      
      # Set confidence thresholds
      high_conf <- 0.7
      med_conf <- 0.5
      
      # Get top scGate scores
      sorted_scores <- sort(unlist(scgate_scores), decreasing = TRUE)
      top_score <- sorted_scores[1]
      top_label <- names(sorted_scores)[1]
      
      # Determine final annotation
      if (top_score >= high_conf) {
        final_label <- top_label
      } else if (top_score >= med_conf && singler_annot == top_label) {
        final_label <- top_label
      } else if (singler_annot %in% names(scgate_scores) && scgate_scores[[singler_annot]] >= 0.3) {
        final_label <- singler_annot
      } else if (singler_annot == "malignant_cells" && scgate_scores$malignant_cells >= 0.3) {
        final_label <- "malignant_cells"
      } else {
        final_label <- "unknown"
      }
      
      seu$final_annotation[cells_in_cluster] <- final_label
    }
    
    return(seu)
  }
  
  # Run final annotation
  seu <- integrate_annotations(seu)
  print(paste('CONFIDENCE', sample, 'CONFIDENCE'))
  
  # Set final annotation as active identity
  Idents(seu) <- "final_annotation"
  print(seu$final_annotation %>% head())
  return(seu)
  }, error = function(e) {
    print(paste("Error in sample", sample, ":", e$message))
    return(NULL)
  })
}

seu <- qread(paste0(PREPRDATADIR, 'seu_list_preproc_ATLAS_v1.rds'))
annotated_samples <- mclapply(seu$SampleID %>% unique,
                            annotate_neuroblastoma, mc.cores = 5)

names(annotated_samples) <- seu$SampleID %>% unique
qsave(annotated_samples, paste0(PREPRDATADIR, 'seu_list_mypipelineannot_ATLAS.qs'))