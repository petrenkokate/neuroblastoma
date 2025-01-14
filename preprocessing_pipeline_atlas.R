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

bped <- celldex::BlueprintEncodeData()
hum_atlas_data <- celldex::HumanPrimaryCellAtlasData()
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes 
data(EnsemblGeneTable.Hs)
scGate_models_DB <- get_scGateDB()

#Dong2020
group1_dong <- ReadMtx(paste0(DATADIR, 'Dong2020/Group1/Exp_data_UMIcounts1.mtx'),
                       paste0(DATADIR,'Dong2020/Group1/Cells1.csv'), 
                       paste0(DATADIR,'Dong2020/Group1/Genes1.txt'), 
                       feature.column=1, skip.cell=1, cell.sep =',')
group1_meta <- fread(paste0(DATADIR, 'Dong2020/Group1/Cells1.csv'), sep = ',')

group2_dong <- ReadMtx(paste0(DATADIR, 'Dong2020/Group2/Exp_data_UMIcounts2.mtx'),
                       paste0(DATADIR,'Dong2020/Group2/Cells2.csv'), 
                       paste0(DATADIR,'Dong2020/Group2/Genes2.txt'), 
                       feature.column=1, skip.cell=1, cell.sep =',')
group2_meta <- fread(paste0(DATADIR, 'Dong2020/Group2/Cells2.csv'), sep = ',')

create_matrices_from_mtx<- function(meta_data, data_matrix, variable_name) {
  sample_names <- unique(meta_data$sample)
  invisible(lapply(sample_names, function(sample_name) {
    sample_cells <- meta_data[sample == sample_name, cell_name]
    sample_matrix <- data_matrix[, sample_cells, drop = FALSE]
    assign(paste0(sample_name, variable_name), sample_matrix, envir = .GlobalEnv)
  }))
}

print('DONG READ')
create_matrices_from_mtx(group1_meta, group1_dong, '_dong')
create_matrices_from_mtx(group2_meta, group2_dong, '_dong')
rm(group1_meta, group2_meta, group1_dong, group2_dong)

#Jansky2021
matrix_jansky <- ReadMtx(paste0(DATADIR, 'Jansky2021/Exp_data_UMIcounts.mtx'),
                         paste0(DATADIR,'Jansky2021/Cells.csv'), 
                         paste0(DATADIR,'Jansky2021/Genes.txt'), 
                         feature.column=1, skip.cell=1, cell.sep =',')
jansky_meta <- fread(paste0(DATADIR, 'Jansky2021/Cells.csv'), sep = ',')

print('JANSKY READ')
create_matrices_from_mtx(jansky_meta, matrix_jansky, '_jansky')
rm(jansky_meta, matrix_jansky)

#Verhoeven2022
load_matrices <- function(directory, pattern_files, sep, pattern_find, pattern_name) {
  file_names <- list.files(path = directory, pattern = pattern_files, full.names = TRUE)
  matrices <- lapply(file_names, function(file) {
    mat <- as.matrix(fread(file, sep = sep), rownames = 1)
    assign(gsub(pattern_find, pattern_name, basename(file)), mat, envir = .GlobalEnv)
  })
}

print('Verhoeven READ')
load_matrices(paste0(DATADIR, 'Verhoeven_GSE147766/'), "\\.csv$", sep=',', ".*_(NB\\d+).count.csv$", "\\1_verhoeven")

print('GROSSMANN READ')
# Get all HTA directories
dirs <- list.files(path = paste0(DATADIR, 'Grossman'), 
                   pattern = "^HTA", 
                   full.names = TRUE,
                   include.dirs = TRUE)

# Read the mapping table
sample_map <- read.csv(paste0(DATADIR, "Grossman/Grossman_str2.csv"))

# Function to create Seurat object for one directory
grossmann_read <- function(dir_path) {
  # Extract directory name as biospecimen_id
  biospecimen_id <- basename(dir_path)
  
  # Read the data
  counts <- Read10X(data.dir = dir_path,
                    gene.column = 1)
  #change name of genes to Symbol
  rownames(counts) <- stringr::str_remove(rownames(counts), "\\.[0-9]+$")
  rownames(counts) <- mapIds(org.Hs.eg.db,
                             keys = rownames(counts),
                             column = "SYMBOL",
                             keytype = "ENSEMBL")
  counts <- counts[!is.na(rownames(counts)), ]
  rownames(counts) <- make.unique(rownames(counts))
  # Add SampleID to metadata using biospecimen_id mapping
  SampleID <- paste0(
    sample_map$Patient_No[match(biospecimen_id, sample_map$HTAN_biospecimen_id)],
    "_grossmann"
  )
  invisible(assign(SampleID, counts, envir = .GlobalEnv))
}

lapply(dirs, grossmann_read)

print('WIENKE READ')

# wienke data
load_wienke_matrices <- function(directory, pattern_find="scR_CEL-Seq2(_\\d+)?", pattern_name="wienke") {
  # Read data and metadata
  data <- readRDS(paste0(directory, 'neuroblastoma_nectin2_tigit/objects/geneTxpTable.RDS'))
  meta <- read.table(paste0(directory, 'neuroblastoma_nectin2_tigit/objects/metadataTable.txt'),
                     header=TRUE, sep=" ", quote="\"", stringsAsFactors=FALSE)
  
  # Create mapping
  plate_mapping <- meta$file_id_GEO
  names(plate_mapping) <- rownames(meta)
  
  # Clean sample IDs
  sample_ids <- gsub(pattern_find, pattern_name, plate_mapping)
  sample_ids <- gsub("scR-CEL-Seq2(_\\d+)?", pattern_name, sample_ids)
  
  # Clean gene names
  rownames(data) <- make.unique(gsub(".*--", "", rownames(data)))
  
  # Split data by sample
  unique_samples <- unique(sample_ids)
  invisible(lapply(unique_samples, function(sample) {
    cells <- names(sample_ids)[sample_ids == sample]
    mat <- data[, cells, drop=FALSE]
    assign(sample, mat, envir=.GlobalEnv)
  }))
  
}

load_wienke_matrices(paste0(DATADIR, 'Wienke_data/'))

print('OLSEN READ')
# olsen dataset
load_olsen_matrices <- function(directory) {
  # Define file pattern and mapping
  file_pattern <- "GSM.*_NB\\d+\\.count\\.csv\\.gz$"
  sample_mapping <- c(
    "NB01" = "NB_01_olsen",
    "NB02" = "NB_02_olsen",
    "NB12" = "NB_12_olsen", 
    "NB16" = "NB_16_olsen",
    "NB17" = "NB_17_olsen",
    "NB19" = "NB_19_olsen",
    "NB21" = "NB_21_olsen",
    "NB26" = "NB_26_olsen"
  )
  
  # Get file paths
  files <- list.files(path=directory, pattern=file_pattern, full.names=TRUE)
  
  # Read and process each file
  invisible(lapply(files, function(file) {
    # Extract NB number from filename
    nb_id <- stringr::str_extract(basename(file), "NB\\d+")
    # Read gzipped CSV
    mat <- as.matrix(fread(file, sep=","), rownames=1)
    # Assign to global environment with mapped name
    assign(sample_mapping[nb_id], mat, envir=.GlobalEnv)
  }))
}

# Usage
load_olsen_matrices(paste0(DATADIR, "olsen2024_GSE147766/"))

print('BONINE READ')

# bonine dataset
load_bonine_matrices <- function(directory) {
  sample_mapping <- c(
    "Bonine2023_cell_KDP02" = "NB_001_bonine",
    "Bonine2023_cell_KDP06" = "NB_002_bonine",
    "Bonine2023_CITE_CS202" = "NB_003_bonine_sc",
    "Bonine2023_nucleus_CS210" = "NB_003_bonine_sn",
    "Bonine2023_nucleus_NBO001" = "NB_004_bonine_1",
    "Bonine2023_nucleus_NBO004" = "NB_004_bonine_2",
    "Bonine2023_nucleus_NBO006" = "NB_006_bonine"
  )
  
  lapply(names(sample_mapping), function(sample_id) {
    barcode_file <- list.files(directory, pattern=paste0(sample_id, "_barcodes.tsv.gz$"), full.names=TRUE)
    features_file <- list.files(directory, pattern=paste0(sample_id, "_features.tsv.gz$"), full.names=TRUE)
    matrix_file <- list.files(directory, pattern=paste0(sample_id, "_matrix.mtx.gz$"), full.names=TRUE)
    
    if(length(barcode_file) > 0 && length(features_file) > 0 && length(matrix_file) > 0) {
      # Create temporary directory
      temp_dir <- tempdir()
      sample_dir <- file.path(temp_dir, sample_id)
      dir.create(sample_dir, showWarnings=FALSE)
      
      # Copy files with standard 10x names
      file.copy(barcode_file, file.path(sample_dir, "barcodes.tsv.gz"))
      file.copy(features_file, file.path(sample_dir, "features.tsv.gz"))
      file.copy(matrix_file, file.path(sample_dir, "matrix.mtx.gz"))
      
      mat <- Read10X(data.dir=sample_dir, gene.column=2)
      invisible(assign(sample_mapping[sample_id], mat, envir=.GlobalEnv))
      
      # Cleanup
      unlink(sample_dir, recursive=TRUE)
    }
    return(NULL)
  })
}

load_bonine_matrices(paste0(DATADIR, "Bonine_GSE253865/"))

print('PATEL READ')
# patel dataset
load_patel_matrices <- function(directory) {
  # Get all sample directories
  sample_dirs <- list.dirs(directory, full.names=TRUE, recursive=FALSE)
  sample_dirs <- sample_dirs[!grepl("\\.(sh|txt)$", sample_dirs)]  # Exclude script files
  
  # Create sample name mapping
  create_sample_name <- function(dir_name) {
    htapp_num <- stringr::str_extract(basename(dir_name), "HTAPP-\\d+")
    htapp_num <- gsub("-", "_", htapp_num)
    return(paste0(htapp_num, "_patel"))
  }
  
  # Function to read data using temporary directory
  read_10x_with_temp <- function(data_dir, pattern_prefix) {
    # Find the files
    barcodes_file <- list.files(data_dir, pattern=paste0(pattern_prefix, ".*barcodes\\.tsv\\.gz$"), full.names=TRUE)[1]
    features_file <- list.files(data_dir, pattern=paste0(pattern_prefix, ".*features\\.tsv\\.gz$"), full.names=TRUE)[1]
    matrix_file <- list.files(data_dir, pattern=paste0(pattern_prefix, ".*matrix\\.mtx\\.gz$"), full.names=TRUE)[1]
    
    if(all(file.exists(barcodes_file, features_file, matrix_file))) {
      # Create temporary directory
      temp_dir <- file.path(tempdir(), basename(data_dir))
      dir.create(temp_dir, showWarnings=FALSE, recursive=TRUE)
      
      # Copy files with standard names
      file.copy(barcodes_file, file.path(temp_dir, "barcodes.tsv.gz"))
      file.copy(features_file, file.path(temp_dir, "features.tsv.gz"))
      file.copy(matrix_file, file.path(temp_dir, "matrix.mtx.gz"))
      
      # Read data
      mat <- Read10X(data.dir=temp_dir, gene.column=2)
      
      # Clean up
      unlink(temp_dir, recursive=TRUE)
      
      return(mat)
    }
    return(NULL)
  }
  
  
  for(sample_dir in sample_dirs) {
    base_name <- create_sample_name(sample_dir)
    sample_basename <- basename(sample_dir)
    
    # Check if channel2 exists to determine if we need rep suffixes
    has_channel2 <- dir.exists(file.path(sample_dir, "channel2"))
    
    # Process channel1 (always exists)
    channel1_dir <- file.path(sample_dir, "channel1")
    mat1 <- read_10x_with_temp(channel1_dir, paste0(sample_basename, "_channel1"))
    
    if(!is.null(mat1)) {
      if(has_channel2) {
        # If there's a channel2, add rep1 suffix
        sample_name <- paste0(base_name, "_rep1")
      } else {
        # If only channel1 exists, use base name
        sample_name <- base_name
      }
      invisible(assign(sample_name, mat1, envir=.GlobalEnv))
    }
    
    # Process channel2 if it exists
    if(has_channel2) {
      channel2_dir <- file.path(sample_dir, "channel2")
      mat2 <- read_10x_with_temp(channel2_dir, paste0(sample_basename, "_channel2"))
      
      if(!is.null(mat2)) {
        sample_name <- paste0(base_name, "_rep2")
        invisible(assign(sample_name, mat2, envir=.GlobalEnv))
      }
    }
  }
  
}


load_patel_matrices(paste0(DATADIR, "patel2024/"))

print('CREATE LIST')
all_samples <- ls()[!ls() %in% c("load_matrices", "create_matrices_from_mtx", "load_h5",
                                 "load_wienke_matrices", "dirs", "grossmann_read", 
                                 "my_colors", "sample_map", "load_olsen_matrices",
                                 "load_bonine_matrices", "NB_003_bonine_sc", "load_patel_matrices",
                                 "DATADIR", "WORKDIR", "PREPRDATADIR", 'bped', 'g2m.genes',
                                 'EnsemblGeneTable.Hs', "s.genes", "sample_metadata",
                                 "scGate_models_DB", 'hum_atlas_data')]
samples_list <- mget(all_samples)
rm(list = ls()[!ls() %in% c("my_colors","DATADIR", "WORKDIR", "PREPRDATADIR", 'bped', 'g2m.genes',
                            'EnsemblGeneTable.Hs', "s.genes", "sample_metadata",
                            "scGate_models_DB", 'hum_atlas_data', 'samples_list')])
gc()
seu_list <- lapply(names(samples_list), function(sample_name) {
  sample <- samples_list[[sample_name]]
  print(sample_name)
  seu <- CreateSeuratObject(counts = sample, min.cells = 3, min.features = 200)
  seu$SampleID <- sample_name
  return(seu)
})

names(seu_list) <- names(samples_list)

print('METADATA READ')

sample_metadata <- read.csv(paste0(DATADIR, 'samples_metadata_v5.csv'))
sample_metadata[sample_metadata==""] <- NA

rownames(sample_metadata) <- sample_metadata$Sample_dataset

print('DOUBLETS FINDING')

seu_list <- mclapply(seu_list, function(seu){
  print(unique(seu$SampleID))
  
  seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$percent_ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu[["complexity"]] <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
  message('Normalization is going on...')
  seu <- seu %>% 
    NormalizeData(verbose=F) %>% 
    FindVariableFeatures(verbose=F) %>% 
    ScaleData(vars.to.regress=c('nCount_RNA'), verbose=F) %>% 
    RunPCA(verbose=F) %>% 
    RunUMAP(dims = 1:20, verbose=F)
  # identify doublets
  # scDblFinder
  message('scDblFinder is going on...')
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, dbr.sd=1, verbose=F)
  seu[['scDblFinder.class']] <- sce$scDblFinder.class
  # DoubletFinder
  message('doubletFinder is going on...')
  sweep.list <- paramSweep(seu, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn |>
    arrange(desc(BCmetric))
  pK <- pK[1, 'pK']
  pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
  nExp <- round(ncol(seu) * 0.03)
  seu <- doubletFinder(seu, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:10)
  seu[['DoubletFinder.class']] <- seu@meta.data[, colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]
  seu[[colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]] <- NULL
  seu[[colnames(seu@meta.data)[grepl("pANN", colnames(seu@meta.data))]]] <- NULL
  seu[['singlet']] <- ifelse(seu$scDblFinder.class=='doublet' & seu$DoubletFinder.class=='Doublet', 'no', 'yes')
  
  
  return(seu)
}, mc.cores = 5)


print("PRELIMINARY ANNOTATION")

seu_list <- lapply(seu_list, function(seu) {
  print(unique(seu$SampleID))
  if(grepl("_rep", unique(seu$SampleID))) {
    current_sample <- strsplit(unique(seu$SampleID), "_rep")[[1]][1]  # Handle replicate cases
  }
  else if(grepl("_bonine_", unique(seu$SampleID))) {
    current_sample <- sub("_[^_]+$", "", unique(seu$SampleID))  # Handle bonine cases
  }
  else {current_sample = unique(seu$SampleID)}
  # filter low-quality cells and doublets
  method <- sample_metadata$Method[sample_metadata$Sample_dataset == current_sample]
  
  # Apply filtering based on method
  seu <- if(method == "sc") {
    subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & singlet == 'yes' & complexity > 0.8 & percent_mito < 25)
  } else if(method == "sn") {
    subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & singlet == 'yes' & complexity > 0.8 & percent_mito < 10)
  }
  seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
  seu <- seu%>%
    NormalizeData() |>
    FindVariableFeatures() |> 
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  # Preliminary annotation
  message('SingleR is going...')
  pred_bped_main <- SingleR(test = as.SingleCellExperiment(seu), ref = hum_atlas_data, labels = hum_atlas_data$label.main, BPPARAM=MulticoreParam(30))
  seu[['celltype_humatlas_main']] <- pred_bped_main$pruned.labels
  pred_bped_fine <- SingleR(test = as.SingleCellExperiment(seu), ref = hum_atlas_data, labels = hum_atlas_data$label.fine, BPPARAM=MulticoreParam(30))
  seu[['celltype_humatlas_fine']] <- pred_bped_fine$pruned.labels
  # scGate
  message('scGate is going...')
  seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
})

print('MERGING')
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)], add.cell.ids = names(seu_list))
rm(seu_list)
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
seu$celltype_humatlas_main[is.na(seu$celltype_humatlas_main)] <- 'unknown'
seu$celltype_humatlas_fine[is.na(seu$celltype_humatlas_fine)] <- 'unknown'

print('ANNOTATION FILT')
cells_names <- table(seu$celltype_humatlas_main) %>% 
  as.data.frame %>% 
  filter(Freq > 50) 

seu$celltype_humatlas_main_filt <- ifelse(seu$celltype_humatlas_main %in% cells_names[,1],
                                          seu$celltype_humatlas_main, "unknown")

print('SAVE DATA')
qsave(seu, paste0(PREPRDATADIR, 'seu_list_preproc_ATLAS_v1.rds'))

print('SAVED SUCCESFULLY :)')
