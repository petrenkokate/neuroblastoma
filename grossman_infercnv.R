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

print('READ DATA')
seu <- qread(paste0(PREPRDATADIR, 'seu_list_preproc_grossman.qs'))

print('SET PARAMETERS FOR INFERCNV')

options(scipen = 100)
options(bitmapType="Xlib")
options(bitmapType='cairo')
plan("sequential")

calculate_cnv <- function(sample) {
  
  seu_obj <- subset(seu, SampleID == sample)
  print(paste(sample, 'is going on...'))
  
  # prepare annotation file
  cell_annotations <- seu_obj$celltype_humatlas_main_filt
  
  # Create a data frame with cell names and their annotations
  annotation_df <- data.frame(
    CellName = names(cell_annotations),
    CellType = cell_annotations,
    stringsAsFactors = FALSE) %>% 
    mutate(CellType = ifelse(CellType %in% c("Neurons", "Neuroepithelial_cell", "unknown"), 
                             "malignant", CellType))
  
  malignant_count <- sum(annotation_df$CellType == "malignant")
  
  if (malignant_count < 2) {
    print(paste(sample, 'SKIPPED'))
    return(NULL)
  }
  
  # Count the number of cells in each CellType
  cell_type_counts <- annotation_df %>%
    group_by(CellType) %>%
    summarise(Count = n())
  
  # Find cell types with fewer than 2 cells
  small_cell_types <- cell_type_counts %>%
    filter(Count < 10) %>%
    pull(CellType)
  
  # Update CellType to "unknown" for small cell types
  annotation_df <- annotation_df %>%
    mutate(CellType = ifelse(CellType %in% small_cell_types, "unknown", CellType))
  
  write.table(annotation_df, paste0(WORKDIR, "refs/", sample, "_annotation_file.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cell_annotations <- annotation_df$CellType %>% unique() 
  cell_annotations <- cell_annotations[!cell_annotations %in% c("malignant", "unknown")]
  
  # create the infercnv object
  infercnv_obj = try(CreateInfercnvObject(raw_counts_matrix=seu_obj[["RNA"]]$counts %>% as.matrix(),
                                          annotations_file=paste0(WORKDIR, "refs/", sample, "_annotation_file.txt"),
                                          delim="\t",
                                          gene_order_file=paste0(WORKDIR, 'refs/gencode_v19_gene_pos.txt'),
                                          ref_group_names=cell_annotations))
  
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = try(infercnv::run(infercnv_obj,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir=paste0(PREPRDATADIR, 'infercnv/', sample),  # dir is auto-created for storing outputs
                                   cluster_by_groups=F,   # cluster
                                   plot_steps = T,
                                   num_ref_groups = length(cell_annotations),
                                   analysis_mode='subclusters',
                                   denoise=T,
                                   HMM=T,
                                   num_threads = 40,
                                   png_res = 300,
                                   write_expr_matrix = T
  ) )
  
  return(infercnv_obj)
  
}
print('RUN RUN RUN RUN RUN RUN')
infercnv_list <- mclapply(seu$SampleID %>% unique, calculate_cnv, mc.cores = 2)

print('SAVING LIST')
qsave(infercnv_list, paste0(PREPRDATADIR, 'infercnv_list_grossman.qs'))

print('DONE SUCCESSFULLY! :)')

