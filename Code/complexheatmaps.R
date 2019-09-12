#************************************
# Build heatmaps at OTU/ASV and varying taxonomy levels
#************************************
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
# library(vegan)
library(reshape2)
# library(gplots)
# library(pheatmap)
library(grid)



# --------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

# library(devtools)
# install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap) # Make Complex Heatmaps
# --------------------
# install.packages("circlize")
library(circlize)  # circular visualization in R


####################################
# Define various colour palettes
# Various colour palettes
my_colour_palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
my_colour_palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
my_colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
my_colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
my_colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
my_colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
my_colour_palette_32_distinct <- c("#ea7e00","#ca0074","#d1c69b","#474007","#bb00ad","#9c80ff","#be3300","#542e72","#00b9f5","#09436b","#8b0036","#9ac8e6","#ff1059","#959eff","#154a11","#0290f4","#ff7762","#7dbf00","#ff8194","#834c00","#006e73","#f9bb5d","#d6c943","#017229","#00d3a8","#732427","#36e191","#6a8200","#efb3ea","#3227bb","#ff90e1","#e92a12")
# lesion_palette_7 <- c("#8558d6","#6ee268","#d247ad","#c9d743","#d7453e","#59a237","#d78f2a")
# patient_palette_45 <- c("#d64530","#585fb1","#795d97","#9e4773","#3f6921","#71692c","#a2b93c","#d571cc","#9b3e97","#33947a","#98ad66","#448a4e","#869ae0","#5ce7af","#e085a3","#dfdc87","#d19be2","#5cb735","#e38269","#3db6c0","#50b565","#50902c","#a98a2c","#dde84a","#db3d76","#5fe485","#7c8329","#b3e791","#6fe965","#5ebce9","#3c86c1","#2a6a45","#65b688","#6651d1","#af4ed3","#df872f","#56e4db","#737cea","#ac464b","#dd37b5","#995b2b","#daac6f","#92e2be","#a2e24b","#e0be3a")
my_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
my_colour_palette_10_soft <- c("#9E788F","#4C5B61","#678D58","#AD5233","#A0A083","#4D456A","#588578","#D0AC4C","#2A7BA0","#931621")
####################################
# Function
log_matrix <- function(mymat){
  out <- log(mymat, 10)
  out[is.infinite(out)] <- 0
  return(out)
}


gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

# Center log ratio transform
clr = function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}


filter_heatmap_matrix <- function(myheatmap, row_max = 0, prevalence = 0){
  internal_heatmap <- myheatmap
  internal_heatmap <- internal_heatmap[which(apply(internal_heatmap, 1, max) >= row_max), ]
  # keep only OTUs/taxa that are in more than this fraction of samples
  filter_fraction <- prevalence
  entry_prevalences <- apply(internal_heatmap, 1, function(x) {length(which(x > 0))})/dim(internal_heatmap)[2]
  entries_from_prevalences <- names(entry_prevalences)[entry_prevalences > filter_fraction]
  entries_from_prevalences <- entries_from_prevalences[!is.na(entries_from_prevalences)]
  return(internal_heatmap[entries_from_prevalences,])
}

####################################

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/combined_processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Load the OTU - taxonomy mapping file
# otu_taxonomy_map.df <- read.csv("Result_tables/combined/combined_otu_taxonomy_map.csv", header = T)

# Factorise discrete columns
metadata.df$Commodity <- factor(metadata.df$Commodity)
metadata.df$Sample_type <- factor(metadata.df$Sample_type)
metadata.df$Sample_treatment <- factor(metadata.df$Sample_treatment)

# Load relative abundance matrices
otu_class_rel.m <- as.matrix(read.table(file = "Result_tables/combined/combined_Class_relative_abundances.csv", sep = ",", header = T, row.names = 1))
otu_phylum_rel.m <- as.matrix(read.table(file = "Result_tables/combined/combined_Phylum_relative_abundances.csv", sep = ",", header = T, row.names = 1))

# Cleanup column (sample) names - Relative abundance matrices
colnames(otu_class_rel.m) <- gsub("_J.*", "", colnames(otu_class_rel.m))
colnames(otu_phylum_rel.m) <- gsub("_J.*", "", colnames(otu_phylum_rel.m))

# and correct metadata
rownames(metadata.df) <- gsub("_J.*", "", rownames(metadata.df))
metadata.df$Index <- gsub("_J.*", "", metadata.df$Index)


# Since we likely removed samples from the count matrix
# in the main script, remove them from the metadata.df here
# samples_removed <- metadata.df$Index[!metadata.df$Index %in% colnames(otu_genus_rare_rel.m)]
# metadata.df <- metadata.df[! metadata.df$Index %in% samples_removed,]

# Remove samples that are not in the metadata.
# otu_rare_log.m <- otu_rare_log.m[,colnames(otu_rare_log.m) %in% metadata.df$Index]
# otu_genus_rare_log.m <- otu_genus_rare_log.m[,colnames(otu_genus_rare_log.m) %in% metadata.df$Index]

# Remove samples that are not in the metadata.
# otu_rare_rel.m <- otu_rare_rel.m[,colnames(otu_rare_rel.m) %in% metadata.df$Index,drop=F]
# otu_genus_rel.m <- otu_genus_rare_rel.m[,colnames(otu_genus_rel.m) %in% metadata.df$Index,drop=F]
# otu_family_rel.m<- otu_family_rare_rel.m[,colnames(otu_family_rel.m) %in% metadata.df$Index,drop=F]

# ------------------------------------------------------------------------------------

# Function to create heatmap
make_heatmap <- function(myheatmap_matrix,
                         mymetadata,
                         filename,
                         my_row_labels = NULL,
                         height = 10,
                         width = 10,
                         heatmap_height = 10,
                         heatmap_width = 10,
                         plot_height =10,
                         plot_width =10,
                         variables = NULL, # Annotations
                         cluster_columns = T,
                         cluster_rows = T,
                         my_breaks = NULL,
                         my_palette = NULL,
                         palette_choice = NULL,
                         legend_title = NULL,
                         ...
                         ){
  
  # Assign internal objects
  internal_heatmap_matrix.m <- myheatmap_matrix
  internal_metadata.df <- mymetadata
  
  # Order/filter the heatmap matrix to order/entries of metadata
  internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,rownames(internal_metadata.df)]
  
  # Order the heatmap matrix by the variables
  internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,do.call(order, internal_metadata.df[variables])]
  
  # Order the metadata by the variables
  internal_metadata.df <- internal_metadata.df[do.call(order, internal_metadata.df[variables]),]

  # Create metadata just containing the variables
  metadata_just_variables <- internal_metadata.df[,variables, drop = F]
  
  # Check that rownames match colnames
  if (!all(rownames(internal_metadata.df) == colnames(internal_heatmap_matrix.m))){
    stop("Row names in metadata do not match column names in matrix")
  }
  
  # Create annotations
  
  colour_lists <- list()
  for (myvar in variables){
    var_colour_name <- paste0(myvar, "_colour")
    # Assumes there is a colour column for each variable in the metadata
    # If there is no colour column, create one and assign from palette
    # internal_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
    # internal_colour_palette_10_distinct <- my_colour_palette_20_distinct
    internal_colour_palette_10_distinct <- my_colour_palette_206_distinct
    
    if (!var_colour_name %in% names(internal_metadata.df)){
      myvar_values <- factor(as.character(sort(unique(internal_metadata.df[,myvar]))))
      myvar_colours <- setNames(internal_colour_palette_10_distinct[1:length(myvar_values)], myvar_values)
      all_variable_colours <- as.character(lapply(as.character(internal_metadata.df[,myvar]), function(x) myvar_colours[x]))
      internal_metadata.df[,paste0(myvar,"_colour")] <- all_variable_colours
    }
    
    metadata_subset <- unique(internal_metadata.df[,c(myvar, var_colour_name)])
    # Order by the variable column
    metadata_subset <- metadata_subset[order(metadata_subset[,myvar]),]
    # Factorise the variable column
    metadata_subset[,myvar] <- factor(metadata_subset[,myvar])
    metadata_subset <- metadata_subset[!is.na(metadata_subset[,myvar]),]
    named_colour_list <- setNames(as.character(metadata_subset[, var_colour_name]), as.character(metadata_subset[,myvar]))
    colour_lists[[myvar]] <- named_colour_list
  }
  
  ha <- HeatmapAnnotation(df = metadata_just_variables,
                          which = "column",
                          col = colour_lists,
                          gp = gpar(col = "black",lwd =.2),
                          gap = unit(.1,"cm"))
  
  if (is.null(my_palette)){
    if (is.null(palette_choice)) {palette_choice <- "blue"}
    if (!palette_choice %in% c("blue", "purple","red")) { palette_choice <- "blue"}
    if (palette_choice == "blue"){
      my_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))
    } 
    else if (palette_choice == "purple"){
      my_palette <- colorRampPalette(c("white", "#f9cdac","#f3aca2", "#ee8b97", "#e96a8d","#db5087","#b8428c", "#973490", "#742796","#5e1f88", "#4d1a70", "#3d1459","#2d0f41"))
    } else if (palette_choice == "red"){
      my_palette <- colorRampPalette(c("white", "#fded86","#fde86e", "#f9d063", "#f5b857","#f0a04b","#eb8a40", "#e77235","#e35b2c", "#c74e29","#9d4429","#753c2c","#4c3430"))
    } 
  } else{
    my_palette <- colorRampPalette(my_palette)
  }
  
  if (!is.null(my_breaks)){
    internal_breaks <- my_breaks
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = my_palette(length(internal_breaks)))
    
  } else{
    internal_breaks <- seq(min(internal_heatmap_matrix.m), max(internal_heatmap_matrix.m), length.out = 6)
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = my_palette(length(internal_breaks)))
  }
  
  my_row_labels.v = rownames(internal_heatmap_matrix.m)
  if (!is.null(my_row_labels)){
    my_row_labels.v <- as.character(lapply(my_row_labels.v, function(x) as.character(row_labels.df[row_labels.df[,1] == x,][,2])))
  }
  # Order the heatmap rows by the row labels names
  internal_heatmap_matrix.m <- internal_heatmap_matrix.m[order(my_row_labels.v),]
  my_row_labels.v <- my_row_labels.v[order(my_row_labels.v)]
  
  hm <- Heatmap(matrix = internal_heatmap_matrix.m,
                
                top_annotation = ha,
                
                # Colours
                col = col_fun,
                na_col = "grey",
                
                # Sizing
                show_heatmap_legend = F,
                row_names_max_width = unit(35,"cm"),
                row_labels = my_row_labels.v,
                # row_names_side = "left",
                # height = unit(height,"cm"),
                # width = unit(width,"cm"),
                # heatmap_height = unit(heatmap_height,"cm"),
                # heatmap_width = unit(heatmap_width,"cm"),
                # heatmap_width = unit(15,"cm"),
                
                # Titles
                column_title = "Sample",
                column_title_side = "bottom",
                row_title = "Taxa",
                row_title_side = "left",
                
                # Clustering
                cluster_columns = cluster_columns,
                cluster_rows = cluster_rows,
                clustering_method_columns = "average",
                clustering_method_rows = "average",
                show_column_dend = T, 
                show_row_dend = T,
                
                # Borders
                border = F,
                rect_gp = gpar(col = "white", lwd = 1),
                
                # Text appearance
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                ...
  )
  
  # Legend appearance
  
  hm_legend <- Legend(col_fun = col_fun,
                      at = internal_breaks,
                      title_position = "leftcenter-rot",
                      # grid_width= unit(.1, "cm"),
                      # title = "Relative abundance (%)",
                      title = legend_title,
                      direction = "vertical",
                      # title_position = "topcenter",
                      border = "black",
                      # legend_width = unit(7,"cm"),
                      )
  
  pdf(filename,height=plot_height,width=plot_width)
  draw(hm, annotation_legend_list = c(hm_legend))
  dev.off()
}

# metadata.df$Commodity <- factor(metadata.df$Commodity)
# metadata.df$Sample_type <- factor(metadata.df$Sample_type)
# metadata.df$Sample_treatment <- factor(metadata.df$Sample_treatment)

# Define the discrete variables
discrete_variables <- c("Commodity","Sample_type","Sample_treatment","study_accession")
# discrete_variables <- c("Commodity","study_accession")

# heatmap_class_rel.m <- filter_heatmap_matrix(otu_class_rel.m, row_max = 0.00, prevalence = 0.0)
# heatmap_family_rel.m <- filter_heatmap_matrix(otu_family_rel.m, row_max = 0.05, prevalence = 0.2)
# heatmap_genus_rel.m <- filter_heatmap_matrix(otu_genus_rel.m, row_max = 0.05, prevalence = 0.01)

# new_row_names <- unlist(lapply(rownames(heatmap_otu_rel.m), function(x) {paste0(x, "; ", otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$Genus)}))
# row_labels.df <- data.frame("Row_label" = rownames(heatmap_otu_rel.m), "Row_label_new" = new_row_names)


make_heatmap(otu_phylum_rel.m*100, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/combined/phylum_relative_abundance.pdf"),
             variables = discrete_variables,
             plot_height = 8,
             plot_width = 120,
             cluster_columns = F,
             cluster_rows = T,
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,1,.1))*100,
             legend_title = "Relative abundance %",
             palette_choice = 'purple'
)

make_heatmap(otu_class_rel.m*100, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/combined/class_relative_abundance.pdf"),
             variables = discrete_variables,
             plot_height = 20,
             plot_width = 120,
             cluster_columns = F,
             cluster_rows = T,
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,1,.1))*100,
             legend_title = "Relative abundance %",
             palette_choice = 'purple'
)

# Get the top 10 taxa per project and just show those in the plot




# Calculate the (min, max, mean, median, stdev, #samples) abundances of each taxa at each taxa level
generate_taxa_summary <- function(mydata, taxa_column, group_by_columns = NULL){
  # if (is.null(group_by_columns)){
  #   select_columns <- c(taxa_column, "Read_count", "Relative_abundance")
  # } else{
  select_columns <- c(taxa_column, group_by_columns, "Sample", "study_accession", "Read_count", "Relative_abundance")
  total_samples <- length(unique(mydata$Sample))
  total_projects <- length(unique(mydata$study_accession))
  
  taxa_group_summary <- 
    mydata %>% 
    dplyr::select_(.dots = select_columns) %>%
    dplyr::group_by_(.dots = c(taxa_column, group_by_columns)) %>%
    dplyr::mutate(N_samples = n_distinct(Sample), N_projects = n_distinct(study_accession)) %>% # number of unique samples/index
    dplyr::group_by_(.dots = c(group_by_columns)) %>%
    dplyr::mutate(N_total_samples_in_group = n_distinct(Sample),
                  N_total_projects_in_group = n_distinct(study_accession))  %>%
    dplyr::group_by_(.dots = c(group_by_columns, taxa_column)) %>%
    dplyr::select(-Sample, -study_accession) %>%
    dplyr::summarise(N_samples = max(N_samples),
                     N_total_samples_in_group = max(N_total_samples_in_group),
                     N_projects = max(N_projects),
                     N_total_projects_in_group = max(N_total_projects_in_group),
                     Percent_group_samples = round((max(N_samples) / max(N_total_samples_in_group))*100, 2),
                     Percent_total_samples = round((max(N_samples) / total_samples)*100, 2),
                     Percent_group_projects = round((max(N_projects) / max(N_total_projects_in_group))*100, 2),
                     Percent_total_projects = round((max(N_projects) / total_projects)*100, 2),
                     
                     Mean_read_count = round(mean(Read_count), 2),
                     Median_read_count = median(Read_count),
                     Min_read_count = min(Read_count),
                     Max_read_count = max(Read_count),
                     Summed_read_count = sum(Read_count),
                      
                     Mean_relative_abundance = round(mean(Relative_abundance), 5),
                     Median_relative_abundance = round(median(Relative_abundance), 5),
                     Min_relative_abundance = round(min(Relative_abundance),5),
                     Max_relative_abundance = round(max(Relative_abundance),5),
                     Summed_relative_abundance = round(sum(Relative_abundance),5),
    ) %>%
    as.data.frame()
  return(taxa_group_summary)
}


filter_summary_to_top_n <- function(taxa_summary, grouping_variables, abundance_column, my_top_n = 10){
  # Get the top N taxa as described in a provided taxa summary table.
  out <- 
    taxa_summary %>%
    dplyr::group_by_(.dots = c(grouping_variables)) %>%
    dplyr::arrange(dplyr::desc(get(abundance_column))) %>%
    dplyr::top_n(my_top_n, get(abundance_column)) %>% 
    # dplyr::arrange_(.dots = c(grouping_variables),abundance_column) %>%
    dplyr::arrange_(.dots = c(grouping_variables)) %>%
    as.data.frame()
  return(out)
}

genus_data.df <- read.csv("Result_tables/combined/combined_Genus_counts_abundances_and_metadata.csv",header = T)
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df,taxa_column = "taxonomy_genus",group_by_columns = c("Commodity", "study_accession"))
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)
write.csv(genus_taxa_summary_filtered.df, file = "Result_tables/combined/combined_top_10_genus_Commodity_study_accession.csv",row.names = F, quote = F)

genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df,taxa_column = "taxonomy_genus",group_by_columns = c("Commodity"))
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)
write.csv(genus_taxa_summary_filtered.df, file = "Result_tables/combined/combined_top_10_genus_Commodity.csv",row.names = F, quote = F)
