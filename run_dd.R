library("tibble")
library("gridExtra")
library("stringr")
library("depmap")
library("ExperimentHub")
library(matrixTests)
library(broom)
library(ggrepel)
library(tidyverse)
library(pheatmap)
library("gplots")
library(ggplot2)
library(hrbrthemes)
library(ggpubr)
library(tictoc)
library(dendextend)
library(viridis)
library(treemap)
library(EnhancedVolcano)
library(reshape2)
library(msigdbr)
library(dplyr)
library(httr)


read_EH_file <- function(file_name, containing_folder){
  file <- read.csv(paste(c(containing_folder, file_name, ".csv"), collapse=""))
  row.names(file) <- file[,1]
  file <- file[,!(names(file) %in% "X")]
  return(file)
}

DD_profile <- function(base_mat, mut_mat, cells, only_target = F, return_p_value_for_gene = c(F,NA)){
  library(parallel)
  library(dplyr)
  library(broom)
  
  # Determine which genes to process
  if(return_p_value_for_gene[1]){
    genes <- return_p_value_for_gene[2]
  } else {
    genes <- row.names(base_mat)
  }
  
  # Pre-filter and prepare the CRISPR data once, outside the loop
  gene_stat_base <- CRISPR %>% filter(depmap_id %in% cells)
  
  # Detect number of cores (use one less than available to prevent system overload)
  num_cores <- 10
  
  # Define the function to process a single gene
  process_gene <- function(gene) {
    gene_stat <- gene_stat_base
    gene_data <- data.frame("Sample" = colnames(mut_mat), "Status" = mut_mat[gene,])
    gene_stat <- full_join(gene_stat, gene_data, by = c('depmap_id' = 'Sample'))
    
    gene_stat <- gene_stat %>% mutate(Status=recode(Status,`0`="no",`1`="yes"))
    
    # Nested analysis
    df <- gene_stat %>%
      filter(!is.na(Status)) %>%
      filter(!is.na(gene_name)) %>%
      group_by(gene_name, Status) %>%
      nest() %>% 
      spread(key = Status, value = data) %>% 
      mutate(
        t_test = map2(no, yes, ~{t.test(.x$dependency, .y$dependency) %>%
            broom::tidy()}),
        no = map(no, nrow),
        yes = map(yes, nrow)
      ) %>% 
      unnest()
    
    df <- df %>% mutate(qvalue = p.adjust(p.value, method='fdr'))
    
    # If returning p-value for a specific gene
    if(return_p_value_for_gene[1]){
      return(df)
    }
    
    # The matrix is filtered according to thresholds
    tmp <- df[,c('estimate','qvalue')]
    
    # Apply thresholds
    tmp[-log10(tmp$qvalue) < 2,]$estimate <- 0  # Insignificant cells
    tmp[tmp$estimate < 0.2,]$estimate <- 0      # Small effects
    
    # Return the processed data for this gene
    result <- list(gene = gene, estimates = t(tmp$estimate))
    return(result)
  }
  
  # Process genes in parallel
  if(return_p_value_for_gene[1]) {
    # For single gene p-value requests, don't use parallelization
    return(process_gene(genes))
  } else {
    # Use parallel processing for multiple genes
    results <- mclapply(genes, process_gene, mc.cores = num_cores)
    
    # Combine results
    for(result in results) {
      base_mat[result$gene,] <- result$estimates
    }
    
    if(only_target){
      tmp <- list()
      tmp[[1]] <- base_mat
      tmp[[2]] <- base_mat[,intersect(colnames(base_mat),unique(SURF$UniProt.gene))]
      base_mat <- tmp
    }
    
    return(base_mat)
  }
}

box_profile <- function(profile_mat, mut_mat, driver_gene, type = "pan", cancer_type = NA, goi = NA,
                        dist = F, DD_thresh = 0.2, highlight_cells = NA, paper_version = F){
  tmp <- profile_mat[driver_gene,] != 0
  potential_partners <- names(tmp)[tmp]
  
  # If there are no potential parnters for the driver-gene, we have nothing to show...
  if(length(potential_partners) == 0){
    return(NA)
  }
  
  # Pcell - cells which are positive for a mutation in the driver-gene (MUT).
  # Ncell- cells which are negative for a mutation in the driver-gene (WT).
  tmp <- mut_mat[driver_gene,] == 1
  Pcell <- names(tmp)[tmp]
  tmp <- mut_mat[driver_gene,] == 0
  Ncell <- names(tmp)[tmp]
  
  Pcell <- gsub("\\.", "-", Pcell)
  Ncell <- gsub("\\.", "-", Ncell)
  
  # In a cancer-specific analysis, filter cell-lines that are not relevant to us.
  if(type == "specific"){
    Pcell <- c(Pcell[METADATA[Pcell,]$primary_disease == cancer_type])
    Ncell <- c(Ncell[METADATA[Ncell,]$primary_disease == cancer_type])
  }
  
  # If there are values in the goi (gene of interest) parameter, act as if only the goi(s) is/are the potential partners to the driver-gene.
  if(length(goi) >= 1 || !is.na(goi)){
    potential_partners <- goi
  }
  
  # Iterate through all the potential partners.
  # For each - calculate the dependency of g (the current potential partner) in each cell line.
  # If there are less than 5 cell-lines with/without a mutation in the driver-gene that also have dependency data for g, return NA.
  # Else, we calculate the mean dependency of each group.
  #   *If the calculation of the means failed for some reason, we send a ggplot object containing an error message for the shiny app.
  # If there are multiple genes inside goi, or if there is 1 gene and dist is set to TRUE, we check whether g is found in goi. If it does,
  #   we create pie-charts for the composition of cancer-types in the different parts of the boxplot.
  # Regardless of the last step, we create the dependency boxplot and save/return it (goi is empty - save, goi is not empty - return).
  for(g in potential_partners){ 
    Pcell_dependencies <- CRISPR %>% filter(depmap_id %in% Pcell) %>% filter(gene_name == g)
    Pcell_dependencies <- Pcell_dependencies[,c("depmap_id", "dependency")]
    
    Ncell_dependencies <- CRISPR %>% filter(depmap_id %in% Ncell) %>% filter(gene_name == g)
    Ncell_dependencies <- Ncell_dependencies[,c("depmap_id", "dependency")]
    
    # Check for sufficient number of cell-lines in the driver-gene-MUT and the driver-gene-WT groups after filtering for a specific partner.
    if(isTRUE(nrow(Pcell_dependencies) >= 5) & isTRUE(nrow(Ncell_dependencies) >= 5)){ #5? sure?
      Pcell_dependencies$group <- "MUT"
      P_mean <- mean(Pcell_dependencies$dependency)
      
      Ncell_dependencies$group <- "WT"
      N_mean <- mean(Ncell_dependencies$dependency)
      
      # Checking if the mean calculation succeeded.
      # In the past we used to also check whether the mean-difference was higher than the DD_thresh parameter.
      if(is.na(N_mean) || is.na(P_mean) ){#|| (N_mean-P_mean < DD_thresh)){
        print("is.na(N_mean) || is.na(P_mean)")
        graph <- ggplot(data.frame(), aes(1, 1)) +
          geom_text(aes(label = "Error: not enough instances in database to create graph"),
                    size = 10, hjust = 0.5, vjust = 0.5) +
          theme_void()  # Remove axis and background elements
        
        return(graph)
      }
      
      p_num <- nrow(Pcell_dependencies)
      n_num <- nrow(Ncell_dependencies)
      data <- full_join(Pcell_dependencies, Ncell_dependencies, by=colnames(Pcell_dependencies))
      if(type == "pan"){
        data$cancer_types <- METADATA[data$depmap_id,]$primary_disease
      }
      else{
        # This is not really "cancer-types", it's "cancer-subtypes".
        data$cancer_types <- METADATA[data$depmap_id,]$lineage_subtype
      }
      
      # If there are multiple genes in the goi parameter or if there is 1 gene in it and the dist parameter is set tot TRUE.
      if((length(goi) > 1 || !is.na(goi)) && dist){
        # If the current fene is in goi, we create pie-charts for it.
        if(g %in% goi){
          pies <- list()
          pies[[1]] <- table((data %>% filter(dependency < -0.5) %>% filter(group == "WT"))$cancer_types)
          pies[[1]] <- data.frame(pies[[1]])
          
          pies[[2]] <- table((data %>% filter(dependency > -0.5) %>% filter(group == "WT"))$cancer_types)
          pies[[2]] <- data.frame(pies[[2]])
          
          pies[[3]] <- table((data %>% filter(dependency < -0.5) %>% filter(group == "MUT"))$cancer_types)
          pies[[3]] <- data.frame(pies[[3]])
          
          pies[[4]] <- table((data %>% filter(dependency > -0.5) %>% filter(group == "MUT"))$cancer_types)
          pies[[4]] <- data.frame(pies[[4]])
          
          for(ind in 1:4){
            pies[[ind]] <- ggplot(pies[[ind]], aes(x="", y=Freq, fill=Var1)) +
              geom_bar(width=1, stat="identity", color="white") +
              coord_polar("y", start=0) +
              theme_void() +
              scale_fill_manual(name = "Cancer Types",
                                values = c('#FF0000', '#A98307', '#00FF00', '#00FFFF', '#0000FF', '#FF00FF', '#FF7D00', '#0082FF', '#7D00FF', '#FF0082', '#FF8080', '#666666', '#80FF80', '#FCFC7D','#390056', '#8080FF', '#000000', '#800000', '#808000', '#008000', '#008080', '#000080', '#800080', '#800040', '#00F6B3', '#002256', '#005655', '#564200', '#560000','#E55137','#D0D0D0'))
            plot(pies[[ind]])
          }
        }
      }
      
      data$name <- NA
      for(i in rownames(data)){
        data[i,]$name <- (METADATA %>% filter(depmap_id == data[i,]$depmap_id))$stripped_cell_line_name
      }
      
      # Creating the actual boxplot.
      DD_value <- format(round(N_mean-P_mean, 2), nsmall = 2)
      WT_percentage <- round(n_num/(p_num+n_num)*100, 2)
      MUT_percentage <-round(p_num/(p_num+n_num)*100, 2)
      box <- box_plot_creation(data, type, cancer_type, driver_gene, g, highlight_cells, paper_version,
                               DD_value, WT_percentage, MUT_percentage)
      
      # If there is no goi (we are iterating through all the potential partners), save the boxplot to a file.
      if(is.na(goi)){
        png(paste(c("/bigdata/zivcohen/graphs/boxplots_",type,"/",
                    ifelse(type=="specific",paste(c(make.names(cancer_type),": "),collapse=""),""),
                    driver_gene,"_",g,".png"), collapse=""), height = 480, width = 2*480)
        plot(box)
        dev.off()
      }
      else{
        # If there is a goi, return the boxplot.
        return(box)
      }
    }
    else{
      # This happens when there are not enough cell-lines.
      print(Pcell_dependencies)
      print(Ncell_dependencies)
      
      return(NA)
    }
  }
  
  
}

box_plot_creation <- function(data, type, cancer_type, driver_gene, partner_gene,
                              highlight_cells, paper_version, DD_value, WT_percentage, MUT_percentage){
  if(!paper_version){
    box <- ggplot(data, aes(x=group, y=dependency, fill=group)) +
      geom_boxplot(outlier.shape = NA, show.legend = TRUE) +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      geom_jitter(aes(color=cancer_types), size=2, alpha=0.9) +
      stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
      theme_ipsum() +
      theme(plot.title = element_text(size=11),
            legend.text = element_text(size=10),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.position = "right") +
      geom_signif(comparisons = list(c("WT","MUT")),test = "t.test") +
      ggtitle(label = paste(c(partner_gene, " dependency in ", driver_gene, "+/- ",
                              ifelse(type=="specific", cancer_type, ""),
                              " cells \nDD = ", DD_value), collapse=""),
              subtitle = paste(c("WT: ", WT_percentage,
                                 "%\nMUT: ", MUT_percentage,"%"),collapse="")) +
      geom_point(data = data[data$depmap_id %in% highlight_cells,], color = "#000080") +
      geom_text(aes(label = name),
                data = data[data$depmap_id %in% highlight_cells,],
                color = "#000080", hjust = -.1)
    
    return(box)
  }
  else{
    box <- ggplot(data, aes(x=group, y=dependency, fill=group)) +
      geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      geom_jitter(aes(color=cancer_types), size=2, alpha=0.9) +
      stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
      theme_ipsum() +
      theme(plot.title = element_text(size=15),
            legend.text = element_text(size=10),
            axis.title.x = element_text(size = 14, hjust = 0.5),
            axis.title.y = element_text(size = 14, hjust = 0.5),
            legend.position = "none") +
      xlab(paste(c(driver_gene, "'s Mutation Status"), collapse="")) +
      ylab(paste(c(partner_gene, "'s Dependency Score"), collapse="")) +
      geom_signif(comparisons = list(c("WT","MUT")),test = "t.test") +
      ggtitle(label = paste(c(partner_gene, " dependency in ", driver_gene, "+/- ",
                              ifelse(type=="specific", cancer_type, ""),
                              " cells \nDD = ", DD_value), collapse=""),
              subtitle = "") +
      geom_point(data = data[data$depmap_id %in% highlight_cells,], color = "#000080") +
      geom_text(aes(label = name),
                data = data[data$depmap_id %in% highlight_cells,],
                color = "#000080", hjust = -.1)
    
    return(box)
  }
}

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

enrichment_analysis <- function(GOIs, current_depmap_ids, approved_genes) {
  tmp <- funrar::stack_to_matrix(TPM %>% filter(depmap_id %in% current_depmap_ids),
                                 'gene_name', 'cell_line', 'rna_expression')
  tmp <- tmp[approved_genes,]
  
  ssgsea_scores <- ssgsea(tmp, gene_sets)
  enrichment_df <- melt(ssgsea_scores,
                        varnames = c("gene_set", "cell_line"),
                        value.name = "Enrichment")
  
  expanded_df <- expand.grid(gene_set = unique(enrichment_df$gene_set),
                             cell_line = unique(enrichment_df$cell_line),
                             gene_name = GOIs)
  
  expanded_df <- left_join(enrichment_df, expanded_df,
                           by = c("gene_set" = "gene_set", "cell_line" = "cell_line"))
  
  expanded_df <- left_join(expanded_df, CRISPR,
                           by = c("gene_name" = "gene_name", "cell_line" = "cell_line"))
  
  correlation_df <- expanded_df %>%
    group_by(gene_set, gene_name) %>%
    filter(complete.cases(Enrichment, dependency)) %>%
    summarize(Correlation = cor(Enrichment, dependency, method = "pearson"),
              P_Value = cor.test(Enrichment, dependency, method = "pearson")$p.value)
  
  return(correlation_df)
}

large_enrichment_analysis <- function(GOIs, approved_genes) {
  tmp <- funrar::stack_to_matrix(TPM, 'gene_name', 'cell_line', 'rna_expression')
  tmp <- tmp[approved_genes,]
  
  ssgsea_scores <- ssgsea(tmp, gene_sets)
  enrichment_df <- melt(ssgsea_scores,
                        varnames = c("gene_set", "cell_line"),
                        value.name = "Enrichment")
  
  print(length(unique(enrichment_df$gene_set)))
  print(length(unique(enrichment_df$cell_line)))
  print(length(GOIs))
  
  all_gene_sets <- unique(enrichment_df$gene_set)
  chunk_size <- 20  # Define a reasonable chunk size
  
  correlation_dfs <- list()  # List to store correlation data frames
  
  for (i in seq(1, length(all_gene_sets), chunk_size)) {
    print(paste(i, " to ", min(i + chunk_size - 1, length(all_gene_sets))))
    gene_set_chunk <- all_gene_sets[i:min(i + chunk_size - 1, length(all_gene_sets))]
    gene_set_df <- enrichment_df[enrichment_df$gene_set %in% gene_set_chunk, ]
    
    expanded_df <- expand.grid(
      gene_set = gene_set_chunk,
      cell_line = unique(gene_set_df$cell_line),
      gene_name = GOIs
    )
    
    expanded_df <- left_join(enrichment_df, expanded_df,
                             by = c("gene_set" = "gene_set", "cell_line" = "cell_line"))
    
    expanded_df <- left_join(expanded_df, CRISPR,
                             by = c("gene_name" = "gene_name", "cell_line" = "cell_line"))
    
    # Continue with your data manipulation and analysis steps for each chunk
    View(expanded_df)
    # Store the correlation data frame in the list
    correlation_df <- expanded_df %>%
      group_by(gene_set, gene_name) %>%
      filter(complete.cases(Enrichment, dependency)) %>%
      summarize(Correlation = cor(Enrichment, dependency, method = "pearson"),
                P_Value = cor.test(Enrichment, dependency, method = "pearson")$p.value)
    correlation_dfs[[length(correlation_dfs) + 1]] <- correlation_df
  }
  
  # Combine all correlation data frames from the list
  final_correlation_df <- do.call(rbind, correlation_dfs)
  
  return(final_correlation_df)
  
  expanded_df <- expand.grid(gene_set = unique(enrichment_df$gene_set),
                             cell_line = unique(enrichment_df$cell_line),
                             gene_name = GOIs)
  
  expanded_df <- left_join(enrichment_df, expanded_df,
                           by = c("gene_set" = "gene_set", "cell_line" = "cell_line"))
  
  expanded_df <- left_join(expanded_df, CRISPR,
                           by = c("gene_name" = "gene_name", "cell_line" = "cell_line"))
  
  correlation_df <- expanded_df %>%
    group_by(gene_set, gene_name) %>%
    filter(complete.cases(Enrichment, dependency)) %>%
    summarize(Correlation = cor(Enrichment, dependency, method = "pearson"),
              P_Value = cor.test(Enrichment, dependency, method = "pearson")$p.value)
  
  return(correlation_df)
}

# Set paths (same as in original script)
DEPMAP_FILES_FOLDER <- "./"
COMMON_GENES_FILE <- "ShinyData/common_genes.txt"
PAN_CANCER_DRIVER_GENE_ANALYSIS_PREFIX <- "ShinyData/MUT/Pan/"
CANCER_SPECIFIC_DRIVER_GENE_ANALYSIS_PREFIX <- "ShinyData/MUT/Specific/"
PATHWAY_ANALYSIS_PREFIX <- "ShinyData/Expression/"

# Load minimal required data
TPM <- read_EH_file("TPM", DEPMAP_FILES_FOLDER)
CRISPR <- read_EH_file("CRISPR", DEPMAP_FILES_FOLDER)
METADATA <- read_EH_file("METADATA", DEPMAP_FILES_FOLDER)
MUT <- read_EH_file("MUT", DEPMAP_FILES_FOLDER)
row.names(METADATA) <- METADATA$depmap_id

# Load MSigDB data 
oncogenic_gene_sets <- msigdbr(species = "Homo sapiens", category = "C6")
gene_sets <- as.list(unique(oncogenic_gene_sets$gs_name))
names(gene_sets) <- gene_sets
approved_genes <- intersect(unique(TPM$gene_name), unique(oncogenic_gene_sets$gene_symbol))
for (gs in gene_sets) {
  gene_sets[gs] <- as.list((oncogenic_gene_sets %>%
                              filter(gs_name == gs) %>%
                              filter(gene_symbol %in% approved_genes))
                           ["human_gene_symbol"])
}

# Load druggable genes
druggable_genes <- read.csv(paste(DEPMAP_FILES_FOLDER,"Corsello.csv",sep=""))
druggable_genes <- unique(unlist(strsplit(druggable_genes$target, ", ")))
druggable_genes <- na.omit(druggable_genes)

# Recreate pan_cell_mut_mat_filt
cat("Recreating mutation matrix...\n")
pan_cell_mut <- data.frame("cell_line" = MUT$depmap_id,"gene_name" = MUT$gene_name,"status" = MUT$var_annotation)
pan_cell_mut["status"][pan_cell_mut["status"] == "damaging"] <- 1
pan_cell_mut["status"][pan_cell_mut["status"] != 1] <- 0
pan_cell_mut$status <- as.numeric(pan_cell_mut$status)
pan_cell_mut_mat <- funrar::stack_to_matrix(pan_cell_mut,'gene_name','cell_line','status')
pan_cell_mut_mat[is.na(pan_cell_mut_mat)] = 0

# Filter the matrix
tmp <- colSums(pan_cell_mut_mat)
non_hyper_mut_cells <- names(tmp)[tmp<100]
pan_cell_mut_mat_filt <- pan_cell_mut_mat[,non_hyper_mut_cells]
tmp <- rowSums(pan_cell_mut_mat)
driver_genes <- names(tmp)[tmp>5]
pan_cell_mut_mat_filt <- pan_cell_mut_mat_filt[driver_genes,]
pan_cell_mut_mat_filt[pan_cell_mut_mat_filt>1] = 1

# Free some memory
rm(pan_cell_mut)
rm(tmp)
gc()

# Read common genes
common_genes <- readLines(COMMON_GENES_FILE)

# Recreate cancer_types_filt and filtered_primary_disease
cat("Identifying cancer types...\n")
valid_cells <- intersect(colnames(pan_cell_mut_mat_filt), unique(CRISPR$depmap_id))
valid_cells_df <- data.frame("cell_line" = valid_cells)
valid_cells_df$cancer_type <- METADATA[valid_cells_df$cell_line,]$primary_disease
row.names(valid_cells_df) <- valid_cells_df$cell_line
SURF <- read.csv(paste0('surfaceome.csv'))

cancer_types <- unique(valid_cells_df$cancer_type)
driver_count <- data.frame("cancer_type" = cancer_types, "cell_lines" = NA, "num_drivers" = NA, "drivers" = NA)
row.names(driver_count) <- driver_count$cancer_type

for(cancer in cancer_types){
  cells <- row.names(valid_cells_df[valid_cells_df$cancer_type == cancer,])
  driver_count[cancer, ]$cell_lines <- length(cells)
  
  if(length(cells) > 1){
    cancer_mut_mat <- pan_cell_mut_mat_filt[,cells]
    tmp <- rowSums(cancer_mut_mat)
    cancer_mut_mat <- cancer_mut_mat[names(tmp)[tmp >= 5], ]
    
    driver_count[cancer, ]$num_drivers <- length(row.names(cancer_mut_mat))
    driver_count[cancer, ]$drivers <- paste(sort(row.names(cancer_mut_mat)), collapse = ", ")
  }
}

cancer_types_filt <- driver_count[driver_count$num_drivers > 1 & !is.na(driver_count$num_drivers),]$cancer_type
filtered_primary_disease <- cancer_types_filt
plan(multicore, workers = 20)  
print('run DD')
tmp <- DD_profile(cross_cell_profile, pan_cell_mut_mat_filt, unique(CRISPR$depmap_id), T)
if(class(tmp) == "list"){
  cross_cell_profile <- tmp[[1]]
  cross_cell_profile_target <- tmp[[2]]
} else{
  cross_cell_profile <- tmp
}

# Filter and save
tmp1 <- rowSums(abs(cross_cell_profile))
tmp2 <- colSums(abs(cross_cell_profile))
cross_cell_profile_filt <- cross_cell_profile[names(tmp1)[tmp1>0], names(tmp2)[tmp2>0]]

tmp1 <- rowSums(abs(cross_cell_profile_target))
tmp2 <- colSums(abs(cross_cell_profile_target))
cross_cell_profile_target_filt <- cross_cell_profile_target[names(tmp1)[tmp1>0], names(tmp2)[tmp2>0]]

# Save the pan-cancer files
cat("Saving pan-cancer dependency profiles...\n")
write.table(cross_cell_profile_filt,
            file = paste(PAN_CANCER_DRIVER_GENE_ANALYSIS_PREFIX,"dependencyProfile_matrix.csv",sep=""),
            sep = ",", row.names = TRUE)
write.table(cross_cell_profile_target_filt,
            file = paste(PAN_CANCER_DRIVER_GENE_ANALYSIS_PREFIX,"dependencyProfile_targetable_matrix.csv",sep=""),
            sep = ",", row.names = TRUE)
write.table(pan_cell_mut_mat_filt,
            file = paste(PAN_CANCER_DRIVER_GENE_ANALYSIS_PREFIX,"mutation_matrix.csv",sep=""),
            sep = ",", row.names = TRUE)