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

set.seed(100)
DATADIR = '~/neuroblastoma/data/'
WORKDIR = '~/neuroblastoma/'
PREPRDATADIR = '~/neuroblastoma/preprocessing_data/'

my_colors <- c("#9ACD32FF", "#4682B4FF",  "#DDA0DDFF", "#FFA07AFF", "#8FBC8BFF",
               "#40E0D0FF",  "#F0E68CFF", "#5F9EA0FF", "#D2B48CFF",  
               "#6495EDFF", "#FF69B4FF", "#BA55D3FF",  "#F08080FF", "#32CD32FF",  
               "#FFDAB9FF",  "#87CEEBFF")

print('READ SEU OBJECT')
seu <- qread(paste0(PREPRDATADIR, 'seu_list_preproc_ATLAS_v1.rds'))

print('SET OPTIONS')
options(scipen = 100)
options(bitmapType="Xlib")
options(bitmapType='cairo')
plan("sequential")

calculate_cnv <- function(sample) {
  
  # Skip samples from the olsen dataset
  if (grepl("olsen", tolower(sample), ignore.case = TRUE)) {
    print(paste(sample, 'SKIPPED - Olsen sample'))
    return(NULL)
  }
  
  seu_obj <- subset(seu, SampleID == sample)
  print(paste(sample, 'is going on...'))
  
  # Determine genome version based on sample
  genome_version <- case_when(
    sample %in% c("Jansky2021") ~ "hg19",
    TRUE ~ "hg38"  # Default to hg38 for all other samples
  )
  
  # Select appropriate gene position file
  gene_pos_file <- ifelse(genome_version == "hg19",
                          paste0(WORKDIR, "refs_ATLAS/gencode_v47_hg19_gene_pos.txt"),
                          paste0(WORKDIR, "refs_ATLAS/gencode_v47_hg38_gene_pos.txt"))
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
  
  # Check for minimum number of malignant cells
  malignant_count <- sum(annotation_df$CellType == "malignant")
  if (malignant_count < 2) {
    print(paste(sample, 'SKIPPED - insufficient malignant cells'))
    return(NULL)
  }
  
  # Check for minimum number of reference cells
  ref_count <- sum(annotation_df$CellType != "malignant")
  if (malignant_count < 100) {
    print(paste(sample, 'SKIPPED - insufficient reference cells'))
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
  
  # Save annotation file
  annotation_file <- paste0(WORKDIR, "refs_ATLAS/", sample, "_annotation_file.txt")
  write.table(annotation_df, 
              annotation_file,
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE)
  
  cell_annotations <- annotation_df$CellType %>% unique()
  
  cell_annotations <- cell_annotations[!cell_annotations %in% c("malignant", "unknown")]
  
  # create the infercnv object
  infercnv_obj = try(CreateInfercnvObject(raw_counts_matrix=seu_obj[["RNA"]]$counts %>% as.matrix(),
                                          annotations_file=paste0(WORKDIR, "refs_ATLAS/", sample, "_annotation_file.txt"),
                                          delim="\t",
                                          gene_order_file=paste0(WORKDIR, 'refs_ATLAS/gencode_v19_gene_pos.txt'),
                                          ref_group_names=cell_annotations))
  
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = try(infercnv::run(infercnv_obj,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir=paste0(PREPRDATADIR, 'infercnv_ATLAS/', sample),  # dir is auto-created for storing outputs
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
  
  # Create and run inferCNV object with error handling
  tryCatch({
    # Create inferCNV object
    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = seu_obj[["RNA"]]$counts %>% as.matrix(),
      annotations_file = annotation_file,
      delim = "\t",
      gene_order_file = gene_pos_file,
      ref_group_names = cell_annotations
    )
    
    # Run inferCNV
    infercnv_obj <- infercnv::run(
      infercnv_obj,
      cutoff = 0.1,  # Use 0.1 for 10x-genomics data
      out_dir = out_dir,
      cluster_by_groups = FALSE,
      analysis_mode = 'subclusters',  # Better for malignant cell identification
      denoise = TRUE,
      HMM = TRUE,
      HMM_type = "i3",  # 3-state model (deletion, neutral, amplification)
      HMM_report_by = "cell",  # Get cell-level predictions
      tumor_subcluster_partition_method = "leiden",
      num_threads = 40,
      plot_steps = TRUE,
      plot_probabilities = TRUE,
      png_res = 300,
      write_expr_matrix = TRUE,
      window_length = 101,
      max_centered_threshold = 3,
      noise_logistic = TRUE,
      sd_amplifier = 1.5
    )
    
    # Extract HMM predictions and add to Seurat object
    if (!is.null(infercnv_obj)) {
      # Get HMM state predictions
      hmm_states <- infercnv_obj@hmm_states
      if (!is.null(hmm_states)) {
        # Calculate proportion of amplified or deleted states per cell
        cnv_scores <- sapply(colnames(seu_obj), function(cell) {
          if (cell %in% names(hmm_states)) {
            states <- hmm_states[[cell]]
            # Count proportion of non-neutral states
            mean(states != "neutral")
          } else {
            NA
          }
        })
        
        # Add scores to Seurat object
        seu_obj$cnv_score <- cnv_scores
        
        # Classify cells as malignant based on CNV score
        # Using mean + 2*SD threshold for conservative classification
        score_threshold <- mean(cnv_scores, na.rm = TRUE) + 
          2 * sd(cnv_scores, na.rm = TRUE)
        seu_obj$is_malignant <- seu_obj$cnv_score > score_threshold
        
        # Add subcluster information if available
        if (!is.null(infercnv_obj@tumor_subclusters)) {
          seu_obj$infercnv_subcluster <- infercnv_obj@tumor_subclusters[colnames(seu_obj)]
        }
      }
    }
    
    # Return results
    return(list(
      infercnv_obj = infercnv_obj,
      seurat_obj = seu_obj,
      cnv_scores = cnv_scores
    ))
    
  }, error = function(e) {
    print(paste("Error in sample", sample, ":", e$message))
    return(NULL)
  })
}
  
# Run the analysis
print('Starting inferCNV analysis...')
results <- lapply(seu$SampleID %>% unique(), function(sample) {
  calculate_cnv(sample)
})

# Save results
print('SAVING')
qsave(results, paste0(PREPRDATADIR, "infercnv_results_ATLAS.qs"))

print('SAVED SUCCESSFULLY :)')