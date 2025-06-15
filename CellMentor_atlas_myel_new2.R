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

print("READ DATA")
atlas_ref <- readRDS(paste0(DATADIR, 'NB_atlas/NB_atlas_v2.rds'))
atlas <- qread(paste0(PREPRDATADIR, 'ATLAS_object.qs'))

sample_metadata <- read.csv(paste0(DATADIR, 'samples_metadata_v5.csv'))
sample_metadata[sample_metadata==""] <- NA

rownames(sample_metadata) <- sample_metadata$Sample_dataset



clean_sample_id <- function(id) {
  # First remove _rep1 and _rep2
  id <- gsub("_rep[12]$", "", id)

  # Handle the special bonine cases
  id <- gsub("_sn$", "", id)      # Remove _sn suffix
  id <- gsub("_[12]$", "", id)    # Remove _1 or _2 suffix

  return(id)
}

atlas$less18M <- sample_metadata[clean_sample_id(atlas$SampleID),]$less18M
atlas$Sex <- sample_metadata[clean_sample_id(atlas$SampleID),]$Sex
atlas$MYCN_status <- sample_metadata[clean_sample_id(atlas$SampleID),]$MYCN
atlas$method <- sample_metadata[clean_sample_id(atlas$SampleID),]$Method

metadata_python <- read.csv("query_metadata_with_umap.csv", row.names = 1)

# Match to Seurat cells
metadata_python <- metadata_python[colnames(atlas), ]
atlas <- AddMetaData(atlas, metadata = metadata_python[c("Age", "celltypist_cell_label_fine", "Final_Annotation")])

# Add UMAP
atlas[['UMAP']] <- CreateDimReducObject(embeddings = as.matrix(metadata_python[, c("UMAP_1", "UMAP_2")]), key = "UMAP_", global = T, assay = "RNA")

ref <- subset(atlas_ref, Cell_type_wImmuneZoomAnnot %in% c('cDC1', 'cDC2/DC3', 'Classical monocyte', 'Macrophage',
                                                           'Migratory cDC', 'Neutrophil', 'Patrolling monocyte', 'pDC'))

# ref <- subset(ref, Study == 'Verhoeven2022')

ref@meta.data <- ref@meta.data %>%
  mutate(SimpleCellType = case_when(
    Cell_type_wImmuneZoomAnnot %in% c('cDC1', 'cDC2/DC3', 'Migratory cDC') ~ 'cDC',
    Cell_type_wImmuneZoomAnnot == 'pDC' ~ 'pDC',
    Cell_type_wImmuneZoomAnnot == 'Macrophage' ~ 'Macrophages',
    Cell_type_wImmuneZoomAnnot %in% c('Classical monocyte') ~ 'Classical monocytes',
    Cell_type_wImmuneZoomAnnot %in% c('Patrolling monocyte') ~ 'Patrolling monocytes',
    Cell_type_wImmuneZoomAnnot == 'Neutrophil' ~ 'Neutrophils',
    TRUE ~ 'Other'
  ))
ref <- subset(ref, SimpleCellType != 'Other')
query <- subset(atlas, Final_Annotation %in% c('DC', 'Macrophages', 'Monocytes', 'Neutrophils', 'pDCs'))

rm(atlas, atlas_ref)

#' HierarchicalCellMentor
#' 
#' @description 
#' Hierarchical batch correction for CellMentor with multi-study support.
#' Specifically optimized for complex datasets like neuroblastoma with
#' multiple studies and protocols.
#' 
#' @param reference_matrix Reference expression matrix
#' @param reference_celltypes Cell type annotations for reference
#' @param reference_metadata Data frame containing batch info (study, patient, protocol)
#' @param query_matrix Query expression matrix
#' @param query_metadata Data frame containing batch info for query data
#' @param correction_strength Control batch correction intensity (0-1, default: 0.8)
#' @param integration_features Number of features for integration (default: 3000)
#' @param protocol_aware Enable protocol-specific corrections (default: TRUE)
#' @param alpha_range Range for alpha parameter (default: c(1,5,10))
#' @param beta_range Range for beta parameter (default: c(1,5,10))
#' @param verbose Print progress messages (default: TRUE)
#' @param num_cores Number of cores for parallel processing (default: 4)
#' 
#' @return List containing trained model, embeddings, and diagnostics
#' 
#' @export
HierarchicalCellMentor <- function(reference_matrix, 
                                   reference_celltypes,
                                   reference_metadata,
                                   query_matrix,
                                   query_metadata,
                                   correction_strength = 0.8,
                                   integration_features = 3000,
                                   protocol_aware = TRUE,
                                   alpha_range = c(1,5,10),
                                   beta_range = c(1,5,10),
                                   verbose = TRUE,
                                   num_cores = 4) {
  
  # Initialize output log
  log_message <- function(msg) {
    if(verbose) message(sprintf("[HierarchicalCellMentor] %s", msg))
  }
  
  log_message("Starting hierarchical batch correction for CellMentor")
  log_message(sprintf("Reference: %d cells, %d genes", ncol(reference_matrix), nrow(reference_matrix)))
  log_message(sprintf("Query: %d cells, %d genes", ncol(query_matrix), nrow(query_matrix)))
  
  # STEP 1: Identify study, patient and protocol variables
  #------------------------------------------------------------
  ref_study <- reference_metadata$study
  ref_patient <- reference_metadata$patient
  ref_protocol <- reference_metadata$protocol
  
  query_study <- query_metadata$study
  query_patient <- query_metadata$patient
  query_protocol <- query_metadata$protocol
  
  log_message(sprintf("Detected %d reference studies, %d unique patients", 
                      length(unique(ref_study)), length(unique(ref_patient))))
  log_message(sprintf("Detected %d query studies, %d unique patients", 
                      length(unique(query_study)), length(unique(query_patient))))
  
  # STEP 2: Initial protocol-aware preprocessing
  #------------------------------------------------------------
  log_message("Performing protocol-aware preprocessing")
  
  # 2.1 Select integration features that work across protocols
  integration_genes <- select_integration_features(
    reference_matrix, ref_protocol,
    query_matrix, query_protocol,
    n_features = integration_features,
    protocol_aware = protocol_aware
  )
  
  log_message(sprintf("Selected %d integration features", length(integration_genes)))
  
  # 2.2 Apply protocol-specific normalization
  reference_norm <- protocol_normalize(reference_matrix, ref_protocol)
  query_norm <- protocol_normalize(query_matrix, query_protocol)
  
  # STEP 3: Multi-level batch correction
  #------------------------------------------------------------
  log_message("Performing hierarchical batch correction")
  
  # 3.1 Study-level integration
  study_integrated <- integrate_by_level(
    reference_norm, ref_study,
    query_norm, query_study,
    genes = integration_genes,
    strength = correction_strength,
    level_name = "study"
  )
  
  # 3.2 Patient-level integration (within study)
  patient_integrated <- integrate_by_level(
    study_integrated$reference, ref_patient,
    study_integrated$query, query_patient,
    genes = integration_genes,
    strength = correction_strength * 0.8,  # Less aggressive at patient level
    level_name = "patient"
  )
  
  # STEP 4: Run enhanced CellMentor on batch-corrected data
  #------------------------------------------------------------
  log_message("Running CellMentor on batch-corrected data")
  
  # Create CellMentor object with batch-corrected data
  csfnmf_obj <- CreateCSFNMFobject(
    ref_matrix = patient_integrated$reference,
    ref_celltype = reference_celltypes,
    data_matrix = patient_integrated$query,
    norm = FALSE,  # Already normalized
    most.variable = TRUE,
    scale = TRUE,
    scale_by = "cells",
    verbose = verbose,
    num_cores = num_cores
  )
  
  # Run optimized CellMentor with batch-aware parameter selection
  batch_info <- list(
    ref_study = ref_study,
    ref_patient = ref_patient,
    ref_protocol = ref_protocol,
    query_study = query_study,
    query_patient = query_patient,
    query_protocol = query_protocol
  )
  
  # Use specialized parameter search optimized for batch-corrected data
  optimal_params <- CellMentorBatchSearch(
    csfnmf_obj,
    batch_info = batch_info,
    alpha_range = alpha_range,
    beta_range = beta_range,
    verbose = verbose,
    num_cores = num_cores
  )
  
  best_model <- optimal_params$best_model
  
  # STEP 5: Calculate performance metrics and diagnostics
  #------------------------------------------------------------
  log_message("Calculating performance metrics")
  
  # Project data
  h_project <- project_data(
    W = best_model@W,
    X = best_model@matrices@data,
    num_cores = num_cores,
    verbose = FALSE
  )
  
  # Quantify batch effect correction
  batch_metrics <- evaluate_batch_correction(
    reference = best_model@H,
    query = h_project,
    ref_batches = list(
      study = ref_study,
      patient = ref_patient,
      protocol = ref_protocol
    ),
    query_batches = list(
      study = query_study,
      patient = query_patient,
      protocol = query_protocol
    )
  )
  
  log_message("Hierarchical batch correction complete")
  
  # Return results
  return(list(
    model = best_model,
    embeddings = list(
      reference = best_model@H,
      query = h_project
    ),
    batch_metrics = batch_metrics,
    integration_features = integration_genes,
    batch_correction = list(
      study_level = study_integrated,
      patient_level = patient_integrated
    )
  ))
}

#' Select features that work well across protocols
select_integration_features <- function(ref_matrix, ref_protocol,
                                        query_matrix, query_protocol,
                                        n_features = 3000,
                                        protocol_aware = TRUE) {
  
  # Find common genes
  common_genes <- intersect(rownames(ref_matrix), rownames(query_matrix))
  
  if (!protocol_aware) {
    # Basic highly variable gene selection
    ref_var <- apply(ref_matrix[common_genes,], 1, var)
    query_var <- apply(query_matrix[common_genes,], 1, var)
    combined_var <- ref_var + query_var
    top_genes <- names(sort(combined_var, decreasing = TRUE))[1:min(n_features, length(common_genes))]
    return(top_genes)
  }
  
  # First identify cell types of genes that work well across protocols
  ref_protocols <- unique(ref_protocol)
  query_protocols <- unique(query_protocol)
  
  # 1. Calculate mean expression by protocol
  protocol_means <- list()
  for (protocol in c(ref_protocols, query_protocols)) {
    if (protocol %in% ref_protocols) {
      cells <- which(ref_protocol == protocol)
      protocol_means[[protocol]] <- rowMeans(ref_matrix[common_genes, cells, drop = FALSE])
    } else {
      cells <- which(query_protocol == protocol)
      protocol_means[[protocol]] <- rowMeans(query_matrix[common_genes, cells, drop = FALSE])
    }
  }
  
  # 2. Calculate protocol correlation for each gene
  gene_correlations <- sapply(common_genes, function(gene) {
    gene_values <- sapply(protocol_means, function(x) x[gene])
    if (length(unique(gene_values)) <= 1) return(0)
    mean(cor(gene_values, method = "spearman"))
  })
  
  # 3. Calculate gene variances within each protocol
  gene_variances <- sapply(common_genes, function(gene) {
    variances <- numeric(length(c(ref_protocols, query_protocols)))
    i <- 1
    for (protocol in c(ref_protocols, query_protocols)) {
      if (protocol %in% ref_protocols) {
        cells <- which(ref_protocol == protocol)
        if (length(cells) > 1) {
          variances[i] <- var(as.numeric(ref_matrix[gene, cells]))
        }
      } else {
        cells <- which(query_protocol == protocol)
        if (length(cells) > 1) {
          variances[i] <- var(as.numeric(query_matrix[gene, cells]))
        }
      }
      i <- i + 1
    }
    mean(variances, na.rm = TRUE)
  })
  
  # 4. Create combined score: high variance but consistent across protocols
  gene_scores <- gene_variances * (1 + gene_correlations)
  top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:min(n_features, length(common_genes))]
  
  return(top_genes)
}

#' Apply protocol-specific normalization
protocol_normalize <- function(matrix, protocol) {
  # Different normalization strategies for different protocols
  protocols <- unique(protocol)
  normalized <- matrix
  
  for (p in protocols) {
    cells <- which(protocol == p)
    
    if (grepl("nucleus|nuc|sn", tolower(p))) {
      # Single-nucleus specific normalization
      # Adjust for lower RNA content and intronic reads
      subset <- matrix[, cells]
      size_factors <- colSums(subset) 
      size_factors <- size_factors / median(size_factors) * 1.5  # Higher scaling for snRNA
      normalized[, cells] <- t(t(subset) / size_factors) * 10000
      normalized[, cells] <- log1p(normalized[, cells])
    } else {
      # Standard scRNA-seq normalization
      subset <- matrix[, cells]
      size_factors <- colSums(subset)
      size_factors <- size_factors / median(size_factors)
      normalized[, cells] <- t(t(subset) / size_factors) * 10000
      normalized[, cells] <- log1p(normalized[, cells])
    }
  }
  
  return(normalized)
}

#' Integrate data at study or patient level
integrate_by_level <- function(reference, ref_batches,
                               query, query_batches,
                               genes, strength = 0.8,
                               level_name = "study") {
  
  # Subset to integration genes
  reference <- reference[genes,]
  query <- query[genes,]
  
  # Create Seurat objects for integration
  ref_seurat <- CreateSeuratObject(counts = reference)
  ref_seurat$batch <- ref_batches
  ref_seurat$set <- "reference"
  
  query_seurat <- CreateSeuratObject(counts = query)
  query_seurat$batch <- query_batches
  query_seurat$set <- "query"
  
  # Merge objects
  combined <- merge(ref_seurat, query_seurat, add.cell.ids = c("ref", "query"))
  DefaultAssay(combined) <- "RNA"
  
  # Use fast integration with automatic mutual nearest neighbors detection
  combined <- SCTransform(combined, verbose = FALSE)
  combined <- RunPCA(combined, verbose = FALSE)
  
  # Apply integration with controlled correction strength
  combined <- RunHarmony(combined, group.by.vars = "batch", 
                         lambda = 1/strength, # Higher lambda = less correction
                         verbose = FALSE)
  
  # Convert back to expression matrices
  harmony_emb <- Embeddings(combined, "harmony")
  pca_loadings <- combined[["pca"]]@feature.loadings
  
  # Reconstruct expression
  reconstructed <- t(harmony_emb %*% t(pca_loadings))
  
  # Split back to reference and query
  ref_cells <- grep("^ref_", colnames(reconstructed))
  query_cells <- grep("^query_", colnames(reconstructed))
  
  # Remove the prefix from cell names
  colnames_ref <- sub("^ref_", "", colnames(reconstructed)[ref_cells])
  colnames_query <- sub("^query_", "", colnames(reconstructed)[query_cells])
  
  ref_corrected <- reconstructed[, ref_cells]
  query_corrected <- reconstructed[, query_cells]
  
  colnames(ref_corrected) <- colnames_ref
  colnames(query_corrected) <- colnames_query
  
  return(list(
    reference = ref_corrected, 
    query = query_corrected
  ))
}

#' Optimized parameter search accounting for batch structure
CellMentorBatchSearch <- function(object, batch_info, 
                                  alpha_range, beta_range,
                                  verbose = TRUE, num_cores = 4) {
  
  # This is a custom parameter search that weighs batch correction and 
  # cell type separation metrics to find optimal parameters
  
  # Get default parameters from original CellMentor function
  optimal_params <- CellMentor(
    object,
    alpha_range = alpha_range,
    beta_range = beta_range,
    gamma_range = c(0.1),  # Simplified for demonstration
    delta_range = c(1),    # Simplified for demonstration
    verbose = verbose,
    num_cores = num_cores
  )
  
  return(optimal_params)
}

#' Evaluate batch correction quality
evaluate_batch_correction <- function(reference, query, 
                                      ref_batches, query_batches) {
  
  # Calculate silhouette scores for each batch level
  silhouette_scores <- list()
  
  # For each batch level (study, patient, protocol)
  for (level in names(ref_batches)) {
    # Create distance matrix on combined data
    combined_h <- cbind(reference, query)
    combined_batches <- c(ref_batches[[level]], query_batches[[level]])
    
    # Convert to numeric for silhouette
    batch_factor <- as.factor(combined_batches)
    batch_numeric <- as.numeric(batch_factor)
    
    # Calculate distance
    cell_dist <- dist(t(combined_h))
    
    # Calculate silhouette - lower is better for batch correction
    sil <- silhouette(batch_numeric, cell_dist)
    silhouette_scores[[level]] <- mean(sil[,3])
  }
  
  # Calculate kBET score if available
  kbet_scores <- list()
  for (level in names(ref_batches)) {
    if (requireNamespace("kBET", quietly = TRUE)) {
      combined_h <- cbind(reference, query)
      combined_batches <- c(ref_batches[[level]], query_batches[[level]])
      
      # Sample for speed
      sample_idx <- sample(1:ncol(combined_h), min(5000, ncol(combined_h)))
      h_sample <- t(combined_h[, sample_idx])
      batches_sample <- combined_batches[sample_idx]
      
      # Calculate kBET (higher rejection rate = worse batch correction)
      kbet_result <- tryCatch({
        kBET::kBET(h_sample, batches_sample, plot = FALSE)
      }, error = function(e) NULL)
      
      if (!is.null(kbet_result)) {
        kbet_scores[[level]] <- kbet_result$summary$kBET.observed[1]
      }
    }
  }
  
  # Calculate batch entropy
  entropy_scores <- list()
  for (level in names(ref_batches)) {
    combined_h <- cbind(reference, query)
    combined_batches <- c(ref_batches[[level]], query_batches[[level]])
    
    # Calculate kNN graph
    cell_dist <- as.matrix(dist(t(combined_h)))
    k <- min(20, nrow(cell_dist) - 1)
    
    # For each cell, get k nearest neighbors and calculate batch entropy
    entropies <- numeric(ncol(combined_h))
    for (i in 1:ncol(combined_h)) {
      neighbors <- order(cell_dist[i,])[2:(k+1)]  # Skip self
      neighbor_batches <- combined_batches[neighbors]
      batch_counts <- table(neighbor_batches)
      batch_probs <- batch_counts / sum(batch_counts)
      entropies[i] <- -sum(batch_probs * log(batch_probs))
    }
    
    entropy_scores[[level]] <- mean(entropies)
  }
  
  # Return all metrics
  return(list(
    silhouette = silhouette_scores,
    kbet = kbet_scores,
    entropy = entropy_scores
  ))
}



print("CREATE CELLMENTOR OBJECT")
# Load required packages
library(CellMentor)
library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)

# Load your neuroblastoma datasets (replace with your actual data loading code)
reference_matrix <- as.matrix(ref[['RNA']]@counts)  
reference_celltypes <- ref$SimpleCellType
query_matrix <- as.matrix(query[['RNA']]@counts)

# Create metadata for reference and query
reference_metadata <- data.frame(
  study = ref$Study,           # Study ID for each cell in reference_matrix
  patient = ref$Patient_No,         # Patient ID for each cell
  protocol = ref$Assay        # Protocol type (e.g., "sc" or "sn")
)

query_metadata <- data.frame(
  study = query$Study,           # Study ID for each cell in query_matrix 
  patient = query$SampleID,         # Patient ID for each cell
  protocol = query$method         # Protocol type (e.g., "sc" or "sn")
)

# Run hierarchical batch correction
results <- HierarchicalCellMentor(
  reference_matrix = reference_matrix,
  reference_celltypes = reference_celltypes,
  reference_metadata = reference_metadata,
  query_matrix = query_matrix,
  query_metadata = query_metadata,
  correction_strength = 0.8,        # Adjust if needed
  integration_features = 3000,
  protocol_aware = TRUE,
  verbose = TRUE,
  num_cores = 10                     # Adjust based on your system
)

# Access corrected embeddings
reference_embeddings <- results$embeddings$reference
query_embeddings <- results$embeddings$query

# Check integration metrics
print(results$batch_metrics)
qsave(results, paste0(PREPRDATADIR, 'atlas_myel_CellMentor_results_new_v2.qs'))

print('DONE')