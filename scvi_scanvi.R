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

print("READ DATA")

seu <- qread(paste0(PREPRDATADIR, 'seu_naivetreatment_v2.qs'))

path_to_python <- "/tmp/RtmpWlYYFy/rstudio/terminal/python"
use_python(path_to_python)
sc <- import('scanpy', convert = FALSE)
scvi <- import("scvi", convert = FALSE)
torch <- import("torch")
gpu_available <- torch$cuda$is_available()
cat("GPU available:", gpu_available, "\n")

print('PREPROCESSING')

seu[["RNA"]] <- as(seu[["RNA"]], "Assay")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
top2000 <- head(VariableFeatures(seu), 3000)
seu2000 <- seu[top2000]
seu2000[["RNA"]] <- as(seu2000[["RNA"]], "Assay")
adata <- convertFormat(seu2000, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)

scvi$settings$seed = 100L

print('SCVI')

adata$X = adata$X$tocsr()
# run setup_anndata, use column Study for batch
scvi$model$SCVI$setup_anndata(adata, batch_key = 'SampleID', continuous_covariate_keys = list('nFeature_RNA'))

# create the model
model = scvi$model$SCVI(adata, n_layers=2L, n_latent=30L, gene_likelihood="nb")

model$train()

print('SCVI DONE')

SCVI_LATENT_KEY = "X_scVI"
latent = model$get_latent_representation()

latent <- as.matrix(latent)
rownames(latent) = colnames(seu)
seu[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seu))

print('SCANVI')

scanvi_model = scvi$model$SCANVI$from_scvi_model(
  model,
  adata=adata,
  labels_key="ANNOTATION_FINAL",
  unlabeled_category="Unknown"
)

scanvi_model$train(max_epochs=20L, n_samples_per_label=100L)

print('SCANVI DONE')

latent_scanvi = scanvi_model$get_latent_representation()

latent_scanvi <- as.matrix(latent_scanvi)
rownames(latent_scanvi) = colnames(seu)
seu[["scanvi"]] <- CreateDimReducObject(embeddings = latent_scanvi, key = "scanvi_", assay = DefaultAssay(seu))

print('SAVING')

qsave(seu, paste0(PREPRDATADIR, 'seu_naivetreatment_v2_scvi.qs'))

print('DONE SUCCESSFULLY :)')