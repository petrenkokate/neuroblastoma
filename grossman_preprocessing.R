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

# Get all HTA directories
dirs <- list.files(path = paste0(DATADIR, 'Grossman'), 
                   pattern = "^HTA", 
                   full.names = TRUE,
                   include.dirs = TRUE)
print('READ DATA - Create Seurat objects')
# Function to create Seurat object for one directory
create_seurat <- function(dir_path) {
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
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # Add biospecimen_id to metadata
  seurat_obj$biospecimen_id <- biospecimen_id
  
  return(seurat_obj)
}

# Create list of Seurat objects
seurat_list <- lapply(dirs, create_seurat)

# Name the list elements with biospecimen IDs
names(seurat_list) <- basename(dirs)

print('READ DATA - Merge Seurat objects')
merged_seurat <- merge(seurat_list[[1]],
                       y = seurat_list[-1],
                       add.cell.ids = names(seurat_list))
merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])

sample_map <- read.csv(paste0(DATADIR, "Grossman/Grossman_str2.csv"))

# Add SampleID to metadata using biospecimen_id mapping
merged_seurat$SampleID <- paste0(
  sample_map$Patient_No[match(merged_seurat$biospecimen_id, sample_map$HTAN_biospecimen_id)],
  "_grossman"
)
print('READ DATA - CHECK SAVING')
qsave(merged_seurat, paste0(PREPRDATADIR, 'seu_grossman_CHECK.qs'))

print('READ DATA - Split Seurat objects')
seu_list <- SplitObject(merged_seurat, split.by = "SampleID")

print('DATA PREPROCESSING - Doublets detection')
rm(merged_seurat, seurat_list)
seu_list <- lapply(seu_list, function(seu){
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
  
  print(table(seu$singlet))
  print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "complexity"), ncol = 5) +
          ggtitle(unique(seu$SampleID)))
  
  return(seu)
})

print('DATA PREPROCESSING - Preliminary annotation')
seu_list <- lapply(seu_list, function(seu) {
  print(unique(seu$SampleID))
  # filter low-quality cells and doublets
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & singlet =='yes' & complexity > 0.85)
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
  seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 20)
})

print('DATA PREPROCESSING - Merge seurat objects')
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)], add.cell.ids = names(seu_list))
rm(seu_list)
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
seu$celltype_humatlas_main[is.na(seu$celltype_humatlas_main)] <- 'unknown'
seu$celltype_humatlas_fine[is.na(seu$celltype_humatlas_fine)] <- 'unknown'

print('DATA PREPROCESSING - Add metadata')
sample_metadata <- read.csv(paste0(DATADIR, 'samples_metadata_v4.csv'))
sample_metadata[sample_metadata==""] <- NA

rownames(sample_metadata) <- sample_metadata$Sample_dataset

sample_metadata <- sample_metadata %>% 
  filter(Sample_dataset %in% unique(seu$SampleID))

seu$sex <- sample_metadata[seu$SampleID,]$Sex
seu$less18M <- sample_metadata[seu$SampleID,]$less18M
seu$Survival.Death <- sample_metadata[seu$SampleID,]$Survival.Death
seu$Site <- sample_metadata[seu$SampleID,]$Site
seu$INSS_stage <- sample_metadata[seu$SampleID,]$INSS_stage
seu$Treatment <- sample_metadata[seu$SampleID,]$Treatment

seu$Study <- sapply(str_split(seu$SampleID, "_"), function(x) {
  if (length(x) == 2) {
    return(x[2])
  } else if (length(x) == 3) {
    return(x[3])
  } else {
    return(NA)
  }
})

cells_names <- table(seu$celltype_humatlas_main) %>% 
  as.data.frame %>% 
  filter(Freq > 50) 

seu$celltype_humatlas_main_filt <- ifelse(seu$celltype_humatlas_main %in% cells_names[,1],
                                          seu$celltype_humatlas_main, "unknown")

qsave(seu, paste0(PREPRDATADIR, 'seu_list_preproc_grossman.qs'))

print('DONE')
