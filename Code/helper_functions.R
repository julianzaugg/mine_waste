# Uncomment and run to install the libraries that might be needed 
# install.packages("ggplot2")
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("RColorBrewer")
# install.packages("vegan")
# install.packages("reshape2")
# install.packages("gplots")
# install.packages("heatmap3")
# install.packages("ggfortify")
# install.packages("scales")
# install.packages("seqinr")

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("decontam")
# library(decontam)

# install.packages("devtools") #Installs devtools (if not already installed)
# devtools::install_github("donaldtmcknight/microDecon") #Installs microDecon
# library(microDecon)

# BiocManager::install("DESeq2")
library(DESeq2)

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(gplots)
library(grid)
library(phyloseq)
library(seqinr) # For writing fasta files




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

# Function that calculates the geometric mean with some error-protection bits. 
# DESeq2 does not appear to work (will throw an error) if every OTU (or genus or genome etc.) 
# contains at least one count of zero in every row of the count data.
# Specifically, the function "dds<-DESeq(dds, betaPrior = FALSE)" will fail
# One way to address this is to use the function below as input to DESeq2 to transform the data.
# Calculate the geometric means prior to estimating the size factors
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

df2matrix <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}
# Taxonomy-sample matrix to dataframe convertor
# Just converts a matrix where the row name is a taxonomy label (really can be anything)
m2df <- function(mymatrix, name_of_taxonomy_col = "taxonomy"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, name_of_taxonomy_col] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(name_of_taxonomy_col,cur_names)]
  return(mydf)
}

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
                     
                     Mean_read_count_non_zero = round(mean(Read_count), 2),
                     Mean_read_count = round(sum(Read_count)/N_total_samples_in_group, 2),
                     Median_read_count = median(Read_count),
                     Min_read_count = min(Read_count),
                     Max_read_count = max(Read_count),
                     Summed_read_count = sum(Read_count),
                     
                     Mean_relative_abundance_non_zero = round(mean(Relative_abundance), 5),
                     Mean_relative_abundance = round(sum(Relative_abundance)/N_total_samples_in_group, 5),
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


# ---------------------------------------------------------------------------------------------------------

generate_pca <- function(pca_object, mymetadata, variable_to_plot, colour_palette, limits = NULL, filename = NULL, include_legend = T, add_spider = F, add_ellipse = F,
                         point_alpha = 1, plot_width = 10, plot_height=10, point_size = 0.8, point_line_thickness = 1,
                         label_sites = F, label_species = F,
                         legend_x = NULL, legend_y = NULL, legend_x_offset = 0, legend_y_offset = 0,title_cex = 1,
                         legend_cex = 0.6,
                         plot_spiders = NULL, plot_ellipses = NULL,plot_hulls = NULL, legend_cols = 2, legend_title = NULL,
                         label_ellipse = F, ellipse_label_size = 0.5, ellipse_border_width = 1,variable_colours_available = F, 
                         plot_title = NULL, use_shapes = F, my_levels = NULL,
                         plot_arrows = F, arrow_colour = "black", arrow_alpha = 1,
                         label_arrows=T,arrow_label_size = .5, num_top_species = 5, arrow_scalar = 1,
                         arrow_label_colour = "black", arrow_thickness = .2,arrow_label_font_type = 1,
                         specie_labeller_function = NULL, arrow_label_offset = 0, file_type = "pdf"){
  pca.scores <- try(scores(pca_object, choices=c(1,2,3)))
  if(inherits(pca.scores, "try-error")) {
    return()
  }
  
  # Get component x,y coordinates
  pca_site_scores <- scores(pca_object, display = "sites")
  pca_specie_scores <- scores(pca_object, display = "species")
  pca_percentages <- (pca_object$CA$eig/sum(pca_object$CA$eig)) * 100
  
  # Remove NA entries from the metadata and from the PCA
  internal_metadata <- mymetadata[!is.na(mymetadata[[variable_to_plot]]),]
  pca_site_scores <- pca_site_scores[rownames(pca_site_scores) %in% rownames(internal_metadata),]
  
  if (!is.null(limits)){
    x_min <- limits[1]
    x_max <- limits[2]
    y_min <- limits[3]
    y_max <- limits[4]
  }
  else {
    x_min <- round(lapply(min(pca_site_scores[,1]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    x_max <- round(lapply(max(pca_site_scores[,1]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    y_min <- round(lapply(min(pca_site_scores[,2]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    y_max <- round(lapply(max(pca_site_scores[,2]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    
  }
  my_xlab = paste("PC1 (", round(pca_percentages[1],1), "%)", sep = "")
  my_ylab = paste("PC2 (", round(pca_percentages[2],1), "%)", sep = "")
  metadata_ordered.df <- internal_metadata[order(rownames(internal_metadata)),]
  metadata_ordered.df <- metadata_ordered.df[order(metadata_ordered.df[[variable_to_plot]]),]
  
  # ------------------------------------------------------------------------------------
  # Ensure outcome variable is factored
  # Refactor the variable column so that the levels are consistent
  if (!is.null(my_levels)){
    metadata_ordered.df[,variable_to_plot] <- factor(metadata_ordered.df[,variable_to_plot], levels = my_levels)
  } else{
    # Uncomment to factorise and order levels alphabetically
    metadata_ordered.df[,variable_to_plot] <- factor(metadata_ordered.df[,variable_to_plot], levels = sort(unique(as.character(metadata_ordered.df[,variable_to_plot]))))
    # Uncomment to factorise and inherit the current ordering
    # metadata_ordered.df[,variable_to_plot] <- factor(metadata_ordered.df[,variable_to_plot])
  }
  # ------------------------------------------------------------------------------------
  
  if (!is.null(filename)){
    pdf(filename, height=plot_height,width=plot_width)  
    
  }
  if (!is.null(filename)){
    if (file_type == "pdf"){
      pdf(filename = filename,height = plot_height, width = plot_width)
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      svg(filename = filename,height = plot_height, width = plot_width)
      # svglite(file = filename,height = plot_height, width = plot_width)
    }
  }
  
  plot(0,
       type='n',
       # x = 0, y=0,
       xlim = c(x_min,x_max),
       ylim = c(y_min,y_max),
       xlab = my_xlab,
       ylab = my_ylab)
  
  # Make grid
  grid(NULL,NULL, lty = 2, col = "grey80")
  
  # Assign (unique) colours and shapes for each grouping variable
  variable_values <- levels(metadata_ordered.df[[variable_to_plot]])
  # variable_values <- unique(as.character(metadata_ordered.df[[variable_to_plot]]))
  
  
  # If variable colour column "variable_colour" in metadata, use colours from there
  if (variable_colours_available == T){
    colour_col_name <- paste0(variable_to_plot, "_colour")
    variable_colours <- setNames(as.character(unique(metadata_ordered.df[[colour_col_name]])), as.character(unique(metadata_ordered.df[[variable_to_plot]])))
  } else{
    variable_colours <- setNames(colour_palette[1:length(variable_values)], variable_values)  
  }
  if (use_shapes == T){
    # variable_shapes <- setNames(c(25,24,23,22,21,8,6,5,4,3,2,1)[1:length(variable_values)],variable_values)
    variable_shapes <- setNames(rep(c(25,24,23,22,21),length(variable_values))[1:length(variable_values)],variable_values)
  } else{
    variable_shapes <- setNames(rep(c(21),length(variable_values))[1:length(variable_values)],variable_values)  
  }
  annotation_dataframe <- data.frame(variable_colours, variable_shapes,stringsAsFactors = F)
  annotation_dataframe$variable_outline_colours <- as.character(annotation_dataframe$variable_colours)
  if (point_line_thickness != 0){
    annotation_dataframe[annotation_dataframe$variable_shapes > 15,"variable_outline_colours"] <- "black"  
  }
  
  if (!is.null(my_levels)){
    annotation_dataframe <- annotation_dataframe[my_levels,]
  }
  
  # Order the site scores by the order of the rows in the metadata
  # print(dim(pca_site_scores))
  # print(dim(metadata_ordered.df))
  # print(head(pca_site_scores))
  # print(head(pca_site_scores))
  pca_site_scores <- pca_site_scores[rownames(metadata_ordered.df),]
  
  all_sample_colours <- as.character(
    lapply(
      as.character(metadata_ordered.df[rownames(pca_site_scores),variable_to_plot]), 
      function(x) as.character(annotation_dataframe[x,"variable_colours"])
    )
  )
  # all_sample_colours <- as.character(
  #   lapply(
  #     as.character(metadata_ordered.df[rownames(pca_site_scores),variable_to_plot]), 
  #     function(x) variable_colours[x]
  #   )
  # )
  
  all_sample_shapes <- as.numeric(
    lapply(
      as.character(metadata_ordered.df[rownames(pca_site_scores),variable_to_plot]), 
      function(x) annotation_dataframe[x,"variable_shapes"][[1]]
    )
  )
  
  # all_sample_shapes <- as.numeric(
  #   lapply(
  #     as.character(sort(metadata_ordered.df[rownames(pca_site_scores),variable_to_plot])), 
  #     function(x) variable_shapes[x][[1]]
  #   )
  # )
  
  # Set the outline colours for all samples based on the sample colours and refering to the annotation dataframe created above
  all_sample_outline_colours <- as.character(unlist(lapply(all_sample_colours, function(x) annotation_dataframe[annotation_dataframe$variable_colours == x, "variable_outline_colours"])))
  # Need to construct the legend outline colour vector.
  # legend_point_outline_colours <- annotation_dataframe$variable_outline_colours
  points(pca_site_scores, 
         cex = point_size,
         lwd = point_line_thickness,
         pch = all_sample_shapes,
         col = alpha(all_sample_outline_colours,point_alpha),
         # col = alpha("black",point_alpha),
         bg = alpha(all_sample_colours, point_alpha),
  )
  
  # Plot arrows for species / variables
  plot_arrows_func <- function(){
    
    left_pc1.v <- rownames(pca_specie_scores[order(pca_specie_scores[,1]),][1:num_top_species,])
    right_pc1.v <- rownames(pca_specie_scores[order(pca_specie_scores[,1]),][(length(pca_specie_scores[,1]) - num_top_species):length(pca_specie_scores[,1]),])
    
    left_pc2.v <- rownames(pca_specie_scores[order(pca_specie_scores[,2]),][1:num_top_species,])
    right_pc2.v <- rownames(pca_specie_scores[order(pca_specie_scores[,2]),][(length(pca_specie_scores[,2]) - num_top_species):length(pca_specie_scores[,2]),])
    
    top_vars.v <- unique(c(left_pc1.v, right_pc1.v, left_pc2.v, right_pc2.v))
    arrows(0,0, 
           arrow_scalar * pca_specie_scores[top_vars.v,1], 
           arrow_scalar * pca_specie_scores[top_vars.v,2], 
           length =0.05, 
           col = alpha(arrow_colour, arrow_alpha),
           lwd = arrow_thickness,
           lty = 1)
    
    if (label_arrows){
      if (!is.null(specie_labeller_function)){
        # text(x = pca_specie_scores[top_vars.v,1],
        #      y = pca_specie_scores[top_vars.v,2],
        #      labels = specie_labeller_function(top_vars.v),
        #      cex = .5,
        #      pos = sample(c(1,2,3,4),1))
        for (tv in top_vars.v){
          text(x = pca_specie_scores[tv,1]* arrow_scalar,
               y = pca_specie_scores[tv,2]* arrow_scalar,
               labels = specie_labeller_function(tv),
               cex = arrow_label_size,
               # Values of 1, 2, 3 and 4, respectively indicate positions below, 
               # to the left of, above and to the right of the specified (x,y) coordinates.
               pos = sample(c(1,2,3,4),1),
               # pos = 2,
               offset = arrow_label_offset,
               col = alpha(arrow_label_colour,1),
               font = arrow_label_font_type)
        }
      } else{
        for (tv in top_vars.v){
          text(x = pca_specie_scores[tv,1] * arrow_scalar,
               y = pca_specie_scores[tv,2] * arrow_scalar,
               labels = tv,
               cex = arrow_label_size,
               # Values of 1, 2, 3 and 4, respectively indicate positions below, 
               # to the left of, above and to the right of the specified (x,y) coordinates.
               pos = sample(c(1,2,3,4),1),
               # pos = 2,
               offset = arrow_label_offset,
               col = alpha(arrow_label_colour,1),
               font = arrow_label_font_type)
        }
      }
    }
    
    
  }
  if (plot_arrows == T){
    plot_arrows_func()
  }
  
  # Plot ellipses that are filled
  plot_ellipses_func <- function () {
    for (member in variable_values) {
      if (nrow(metadata_ordered.df[metadata_ordered.df[[variable_to_plot]] == member,]) > 2){ # if too few samples, skip plotting ellipse
        ordiellipse(pca_site_scores,
                    groups = metadata_ordered.df[[variable_to_plot]],
                    kind = "ehull",
                    lwd = ellipse_border_width,
                    # border = variable_colours[member][[1]],
                    border = alpha(variable_colours[member][[1]],point_alpha),
                    # col = variable_colours[member][[1]],
                    col = alpha(variable_colours[member][[1]],point_alpha),
                    show.groups = member,
                    alpha = .05,
                    draw = "polygon",
                    label = F,
                    cex = .5)
      }
    }
  }
  
  # Plot hulls that are filled
  plot_hulls_func <- function () {
    for (member in variable_values){
      if (nrow(metadata_ordered.df[metadata_ordered.df[[variable_to_plot]] == member,]) > 2){ # if too few samples, skip plotting ellipse}
        ordihull(pca_site_scores,
                 groups = metadata_ordered.df[[variable_to_plot]],
                 lwd = ellipse_border_width,
                 # border = variable_colours[member][[1]],
                 border = alpha(variable_colours[member][[1]],point_alpha),
                 # col = variable_colours[member][[1]],
                 col = alpha(variable_colours[member][[1]],point_alpha),
                 show.groups = member,
                 alpha = .05,
                 draw = "polygon",
                 label = F,
                 cex = .5)
      }
    }
  }
  if (hasArg(plot_hulls)){
    if (plot_hulls == T){
      plot_hulls_func()    
    }
  }
  
  plot_ellipses_labels_func <- function(label_ellipse = F){
    # Repeat to have labels clearly on top of all ellipses
    for (member in variable_values){
      if (nrow(metadata_ordered.df[metadata_ordered.df[[variable_to_plot]] == member,]) > 2){ # if too few samples, skip plotting ellipse
        ordiellipse(pca_site_scores,
                    groups = metadata_ordered.df[[variable_to_plot]],
                    kind = "ehull",
                    # border = variable_colours[member][[1]],
                    border = NA,
                    # col = variable_colours[member][[1]],
                    col = NA,
                    show.groups = member,
                    alpha = 0,
                    draw = "polygon",
                    label = label_ellipse,
                    cex = ellipse_label_size)
      }
    }
  }
  
  if (hasArg(plot_ellipses) | hasArg(label_ellipse)){
    if (plot_ellipses == T){
      plot_ellipses_func()    
      plot_ellipses_labels_func(label_ellipse = label_ellipse)
    } else if (label_ellipse == T){
      plot_ellipses_labels_func(label_ellipse = label_ellipse)
    }
  } 
  
  #Plot spiders
  plot_spiders_func <- function (label_spider = F) {
    for (member in variable_values){
      if (nrow(metadata_ordered.df[metadata_ordered.df[[variable_to_plot]] == member,]) > 2){ # if too few samples, skip plotting ellipse
        ordispider(pca_site_scores,
                   groups = metadata_ordered.df[[variable_to_plot]],
                   # col = variable_colours[member][[1]],
                   col = alpha(variable_colours[member][[1]],point_alpha),
                   show.groups = member,
                   #alpha = .05,
                   label = label_spider)
      }
    }
  }
  if (hasArg(plot_spiders)){
    if (plot_spiders == T){
      plot_spiders_func(F)    
    }
  }
  
  # points(pca_specie_scores, 
  #        cex = 0.8,
  #        col = "red",
  #        bg = "red",
  #        pch = 3,
  # )
  
  if (label_sites == T){
    text(x = pca_site_scores[,1],
         y = pca_site_scores[,2],
         labels = rownames(pca_site_scores),
         cex = .5,
         pos = 2)
  }
  if (label_species == T){
    text(x = pca_specie_scores[,1],
         y = pca_specie_scores[,2],
         labels = rownames(pca_specie_scores),
         cex = .5,
         pos = 2)
    # arrows(x = pca_specie_scores[,1],
    #        y = pca_specie_scores[,2])
  }
  
  if (is.null(legend_x) || is.null(legend_y)){
    legend_x <- x_min + legend_x_offset
    legend_y <- y_max + legend_y_offset
  }
  if (is.null(legend_title)){
    legend_title <- variable_to_plot
  }
  
  if (!is.null(plot_title)){
    title(main = plot_title, cex.main = title_cex)
  }
  
  if (include_legend){
    legend(
      # title = bold(variable_to_plot),
      title = as.expression(bquote(bold(.(legend_title)))),
      title.col="black",
      # x = x_min-4,
      # y = y_max-6,
      x = legend_x,
      y = legend_y,
      # legend= variable_values,
      # pch= unique(all_sample_shapes),
      # col= legend_point_outline_colours,
      # pt.bg = unique(all_sample_colours),
      legend= rownames(annotation_dataframe),
      pch= annotation_dataframe$variable_shapes,
      col= as.character(annotation_dataframe$variable_outline_colours),
      pt.bg = as.character(annotation_dataframe$variable_colours),
      #bg = "white",
      bty = "n",
      ncol = legend_cols,
      cex = legend_cex,
      # pt.cex = 0.6,
      pt.lwd = point_line_thickness,
      y.intersp =1,
      x.intersp =1,
    )  
  }
  
  if (!is.null(filename)){
    dev.off()
  }
}


calculate_PC_taxa_contributions <- function(pca_object){
  # Calculate the percentage contribution from each taxa for PC1-3. Requires unscaled values that are squared
  pc1_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,1]^2, 3),value.name = "PC1_contribution_percentage")
  pc2_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC2_contribution_percentage")
  pc3_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC3_contribution_percentage")
  
  data.frame(pc1_contribution, pc2_contribution, pc3_contribution)
}

calculate_PC_abundance_correlations <- function(pca_object, mydata.df, taxa_column, variables = NULL){
  # pca_object <- genus_pca
  # mydata.df <- genus_data.df
  # taxa_column <-  "taxonomy_genus"
  # variables <- discrete_variables
  
  # Assume "Sample" and taxa_column are in mydata.df
  # mydata.df should be the combined table from the specific taxa level
  taxa_contributions <- calculate_PC_taxa_contributions(pca_object)
  
  # PC scores for each sample (site)
  pca_site_scores <- m2df(scores(pca_object, display = "sites"),"Sample") # These will be scaled scores!
  pca_species_scores <- m2df(scores(pca_object, display = "species"),taxa_column) # These will be scaled scores!
  
  # Filter to taxonomy in pca object (because we have filtered during the PCA calculations)
  abundance_pc_scores.df <- subset(mydata.df, get(taxa_column) %in% rownames(pca_object$CA$v))
  abundance_pc_scores.df <- left_join(abundance_pc_scores.df, pca_site_scores, by = "Sample")
  abundance_pc_scores.df <- abundance_pc_scores.df[,c("Sample","Relative_abundance", variables,taxa_column, "PC1", "PC2")]
  
  # Filter out entries where there are no PC scores
  abundance_pc_scores.df <- abundance_pc_scores.df[!is.na(abundance_pc_scores.df$PC1),]
  abundance_pc_scores.df$Relative_abundance <- abundance_pc_scores.df$Relative_abundance*100
  
  # Calculate the correlation between the abundances for each taxa and the PC1 and PC2 scores
  abundance_pc_correlations.df <- 
    abundance_pc_scores.df %>% 
    dplyr::group_by_(taxa_column) %>% dplyr::summarise(Pearson_PC1 = cor(PC1, Relative_abundance, method = "pearson"),
                                                       Pearson_PC2 = cor(PC2, Relative_abundance, method = "pearson"),
                                                       Spearman_PC1 = cor(PC1, Relative_abundance, method = "spearman"),
                                                       Spearman_PC2 = cor(PC2, Relative_abundance, method = "spearman"), 
                                                       N_Samples = n_distinct(Sample)) %>% 
    # filter(N_Samples >= 5) %>% 
    arrange(desc(abs(Pearson_PC1))) %>%
    as.data.frame()
  
  # Combine the correlation table with the PCA specie scores
  abundance_pc_correlations.df <- left_join(abundance_pc_correlations.df, pca_species_scores, by = taxa_column)
  
  
  # Add percentage contributions
  abundance_pc_correlations.df$PC1_contribution_percentage <- as.numeric(lapply(abundance_pc_correlations.df[,taxa_column],  
                                                                                function(x) taxa_contributions[x,]$PC1_contribution_percentage))
  abundance_pc_correlations.df$PC2_contribution_percentage <- as.numeric(lapply(abundance_pc_correlations.df[,taxa_column],  
                                                                                function(x) taxa_contributions[x,]$PC2_contribution_percentage))
  abundance_pc_correlations.df$PC3_contribution_percentage <- as.numeric(lapply(abundance_pc_correlations.df[,taxa_column], 
                                                                                function(x) taxa_contributions[x,]$PC3_contribution_percentage))
  abundance_pc_correlations.df
}
# ---------------------------------------------------------------------------------------------------------

# Function to create heatmap
make_heatmap <- function(myheatmap_matrix,
                         mymetadata = NULL,
                         filename= NULL,
                         #...,
                         # Dataframe with two columns. First must match row entry, second the new label
                         my_row_labels = NULL, 
                         my_col_labels = NULL, # Same as my_row_labels, though with columns
                         height = 10, # Not currently used
                         width = 10, # Not currently used
                         heatmap_height = 10, # Not currently used
                         heatmap_width = 10, # Not currently used
                         plot_height =10,
                         plot_width =10,
                         column_title_size = 10,
                         row_title_size = 10,
                         annotation_bar_name_size = 10,
                         variables = NULL, # Annotations
                         cluster_columns = T,
                         cluster_rows = T,
                         my_breaks = NULL,
                         my_palette = NULL,
                         palette_choice = NULL,
                         legend_title = NULL,
                         column_title = "Sample",
                         row_title = "Taxa",
                         legend_labels = NULL,
                         my_annotation_palette = NULL,
                         discrete_legend = FALSE, # Whether or not to display continuous legend as discrete
                         simple_anno_size = unit(.5, "cm"), # size of annotations
                         show_column_dend = F,
                         show_row_dend = F,
                         do_not_order = F,
                         show_cell_values = F,
                         # If show_cell_values =T, cells less than this will have a black font colour
                         # and above white
                         cell_fun_value_col_threshold = 15,
                         my_cell_fun = NULL,
                         show_legend = T,
                         show_top_annotation = T,
                         row_name_size = 6,
                         col_name_size = 6,
                         col_annotation_label_size = 6,
                         col_annotation_title_size = 6,
                         col_annotation_legend_grid_height = .2,
                         col_annotation_legend_grid_width = .2,
                         grid_thickness = 1,
                         ...){
  # print(list(...))
  argList<-list(...) # argument list for checking unspecified optional parameters
  # print(argList$cell_fun)
  # return(1)
  
  # Assign internal objects
  internal_heatmap_matrix.m <- myheatmap_matrix
  internal_metadata.df <- mymetadata
  
  if (!is.null(internal_metadata.df)){
    # Order/filter the heatmap matrix to order/entries of metadata
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,rownames(internal_metadata.df),drop = F]
    
    if (!is.null(variables)){
      # Order the heatmap matrix by the variables
      internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,do.call(order, internal_metadata.df[,variables,drop=F]),drop =F]
      # Order the metadata by the variables
      internal_metadata.df <- internal_metadata.df[do.call(order, internal_metadata.df[,variables,drop=F]),,drop=F]
      # Create metadata just containing the variables
      metadata_just_variables <- internal_metadata.df[,variables, drop = F]    
    }
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
      if (is.null(my_annotation_palette)){
        internal_colour_palette <- my_colour_palette_206_distinct
      } else{
        internal_colour_palette <- my_annotation_palette
      }
      if (!var_colour_name %in% names(internal_metadata.df)){
        myvar_values <- factor(as.character(sort(unique(internal_metadata.df[,myvar]))))
        myvar_colours <- setNames(internal_colour_palette[1:length(myvar_values)], myvar_values)
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
    
    # Appearance of the column annotations
    #HeatmapAnnotation
    if (show_top_annotation == T){
      if(is.null(variables)){
        print("No variables specified, cannot add annotations")
        show_top_annotation <- F
      } else{
        ha <- columnAnnotation(df = metadata_just_variables,
                               # which = "column",
                               col = colour_lists,
                               gp = gpar(col = "black",lwd =.2),
                               gap = unit(.1,"cm"),
                               show_annotation_name = T,
                               # annotation_legend_param, # ?color_mapping_legend for options
                               show_legend = show_legend,
                               simple_anno_size = simple_anno_size,
                               annotation_legend_param = list(labels_gp = gpar(fontsize = col_annotation_label_size),
                                                              title_gp = gpar(fontsize = col_annotation_title_size,fontface = "bold"),
                                                              grid_height = unit(col_annotation_legend_grid_height,"cm"),
                                                              grid_width = unit(col_annotation_legend_grid_width,"cm")),
                               annotation_name_gp = gpar(fontsize = annotation_bar_name_size))      
      }
      # HeatmapAnnotation(annotation_)
      # ?color_mapping_legend
      
    }
  }
  
  # TODO - add option for row annotation
  if (is.null(my_palette)){
    if (is.null(palette_choice)) {palette_choice <- "blue"}
    if (!palette_choice %in% c("blue", "purple","red","dark_bluered","bluered")) { palette_choice <- "blue"}
    if (palette_choice == "blue"){
      my_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))
    } 
    else if (palette_choice == "purple"){
      my_palette <- colorRampPalette(c("white", "#f9cdac","#f3aca2", "#ee8b97", "#e96a8d","#db5087","#b8428c", "#973490", "#742796","#5e1f88", "#4d1a70", "#3d1459","#2d0f41"))
    } else if (palette_choice == "red"){
      my_palette <- colorRampPalette(c("white", "#fded86","#fde86e", "#f9d063", "#f5b857","#f0a04b","#eb8a40", "#e77235","#e35b2c", "#c74e29","#9d4429","#753c2c","#4c3430"))
    } else if (palette_choice == "dark_bluered"){
      my_palette <- colorRampPalette(c("#08306B","#FFD92F","#67001F"))
    } else if (palette_choice == "bluered"){
      my_palette <- colorRampPalette(c("#17468a","#ffdd47","#99113a"))
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
    my_row_labels.v <- as.character(lapply(my_row_labels.v, function(x) as.character(my_row_labels[my_row_labels[,1] == x,][,2])))
  }
  my_col_labels.v = colnames(internal_heatmap_matrix.m)
  if (!is.null(my_col_labels)){
    my_col_labels.v <- as.character(lapply(my_col_labels.v, function(x) as.character(my_col_labels[my_col_labels[,1] == x,][,2])))
  }
  
  if (do_not_order != T){
    # Order the heatmap rows by the row labels names
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[order(my_row_labels.v),]
    my_row_labels.v <- my_row_labels.v[order(my_row_labels.v)]    
  }
  
  # if show values and no function provided
  if (show_cell_values == T & is.null(my_cell_fun)){ 
    my_cell_fun <- function(j, i, x, y, width, height, fill) {
      # if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold & internal_heatmap_matrix.m[i, j] != 0){
      if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold){
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "black"))}
      else if(internal_heatmap_matrix.m[i, j] >= cell_fun_value_col_threshold ) {
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "white"))
      }
    }
  }
  
  # Legend appearance
  if (is.null(legend_labels)){
    my_labels <- internal_breaks
  } else{
    my_labels <- legend_labels
  }
  if (discrete_legend == TRUE){
    hm_legend <- Legend(
      labels = rev(my_labels),
      at = internal_breaks,
      labels_gp = gpar(fontsize = 6),
      legend_gp = gpar(fill = rev(col_fun(internal_breaks))), # For discrete
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = 6,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  } else{
    hm_legend <- Legend(
      col_fun = col_fun, # For continuous
      labels = my_labels,
      at = internal_breaks,
      labels_gp = gpar(fontsize = 6),
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = 6,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  }
  
  hm <- Heatmap(matrix = internal_heatmap_matrix.m,
                
                # top_annotation = ha,
                
                # Colours
                col = col_fun,
                na_col = "grey",
                
                # Sizing
                show_heatmap_legend = F,
                heatmap_legend_param = list(hm_legend),
                row_names_max_width = unit(35,"cm"),
                row_labels = my_row_labels.v,
                column_labels = my_col_labels.v,
                # row_names_side = "left",
                # height = unit(height,"cm"),
                # width = unit(width,"cm"),
                # heatmap_height = unit(heatmap_height,"cm"),
                # heatmap_width = unit(heatmap_width,"cm"),
                # heatmap_width = unit(15,"cm"),
                
                # Titles
                column_title = column_title,
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize = column_title_size),
                row_title = row_title,
                row_title_side = "left",
                row_title_gp = gpar(fontsize = row_title_size),
                
                # Clustering
                cluster_columns = cluster_columns,
                cluster_rows = cluster_rows,
                clustering_method_columns = "average",
                clustering_method_rows = "average",
                show_column_dend = show_column_dend, 
                show_row_dend = show_row_dend,
                # column_dend_height = unit(2, "cm"),
                # row_dend_width = unit(3, "cm"),
                
                # Borders
                border = F,
                rect_gp = gpar(col = "white", lwd = grid_thickness),
                
                # Text appearance
                row_names_gp = gpar(fontsize = row_name_size),
                column_names_gp = gpar(fontsize = col_name_size),
                cell_fun = my_cell_fun,
                ...
  )
  
  if (show_top_annotation == T & exists("ha")){
    hm <- ha %v% hm
  }
  if (!is.null(filename)){
    pdf(filename,height=plot_height,width=plot_width)
    draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T)
    dev.off()    
  }
  draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T)
  return(list("heatmap" = hm, "legend" = hm_legend))
  
}


# ******************************************************************************************
# Taken from : https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
# ******************************************************************************************


calculate_alpha_diversity_significance <- function(mydata, variable){
  # Assumes there are Shannon, Chao1 and Simpson columns
  # This is run assuming unpaired data.
  results.df <- data.frame("Group_1" = character(),
                           "Group_2" = character(),
                           "Shannon_MannW_pvalue" = character(),
                           "Simpson_MannW_pvalue" = character(),
                           "Chao1_MannW_pvalue" = character(),
                           "Shannon_KrusW_pvalue" = character(),
                           "Simpson_KrusW_pvalue" = character(),
                           "Chao1_KrusW_pvalue" = character())
  group_combinations <- combn(as.character(unique(mydata[,variable])), 2)
  
  for (i in 1:ncol(group_combinations)) {
    group_1 <- group_combinations[1,i]
    group_2 <- group_combinations[2,i]
    group_1_meta <- subset(mydata, get(variable) == group_1)
    group_2_meta <- subset(mydata, get(variable) == group_2)
    
    # Mann-Whitney test on the Shannon diversity
    wilcox_shannon_test <- wilcox.test(group_1_meta$Shannon, group_2_meta$Shannon, exact = F)
    # Mann-Whitney test on the Simpson diversity
    wilcox_simpson_test <- wilcox.test(group_1_meta$Simpson, group_2_meta$Simpson, exact = F)
    # Mann-Whitney test on the Chao1 diversity
    wilcox_chao1_test <- wilcox.test(group_1_meta$Chao1, group_2_meta$Chao1, exact = F)
    
    # Kruskal-Wallis (pairwise) test on the Shannon diversity
    kruskal_shannon_test <- kruskal.test(Shannon~get(variable), data = subset(mydata, get(variable) %in% c(group_1, group_2)))
    # Kruskal-Wallis (pairwise) test on the Simpson diversity
    kruskal_simpson_test <- kruskal.test(Simpson~get(variable), data = subset(mydata, get(variable) %in% c(group_1, group_2)))
    # Kruskal-Wallis (pairwise) test on the Chao1 diversity
    kruskal_chao1_test <- kruskal.test(Chao1~get(variable), data = subset(mydata, get(variable) %in% c(group_1, group_2)))
    
    results.df <- rbind(results.df, data.frame("Group_1" = group_1, 
                                               "Group_2" = group_2, 
                                               "Shannon_MannW_pvalue" = round(wilcox_shannon_test$p.value,6),
                                               "Simpson_MannW_pvalue" = round(wilcox_simpson_test$p.value,6),
                                               "Chao1_MannW_pvalue" = round(wilcox_chao1_test$p.value,6),
                                               "Shannon_KrusW_pvalue" = round(kruskal_shannon_test$p.value,6),
                                               "Simpson_KrusW_pvalue" = round(kruskal_simpson_test$p.value,6),
                                               "Chao1_KrusW_pvalue" = round(kruskal_chao1_test$p.value,6)
    ))
  }
  results.df$Shannon_MannW_padj <- round(p.adjust(results.df$Shannon_MannW_pvalue,method = "BH"),6)
  results.df$Simpson_MannW_padj <- round(p.adjust(results.df$Simpson_MannW_pvalue,method = "BH"),6)
  results.df$Chao1_MannW_padj <- round(p.adjust(results.df$Chao1_MannW_pvalue,method = "BH"),6)
  results.df$Shannon_KrusW_padj <- round(p.adjust(results.df$Shannon_KrusW_pvalue,method = "BH"),6)
  results.df$Simpson_KrusW_padj <- round(p.adjust(results.df$Simpson_KrusW_pvalue,method = "BH"),6)
  results.df$Chao1_KrusW_padj <- round(p.adjust(results.df$Chao1_KrusW_pvalue,method = "BH"),6)
  results.df
}


summarise_alpha_diversities <- function(mydata, group_by_columns){
  summary.df <- mydata %>% 
    dplyr::group_by_(.dots = c(group_by_columns)) %>%
    dplyr::summarise(Shannon_Mean =mean(Shannon),
              Shannon_Stdev=sd(Shannon),
              Shannon_Max=max(Shannon), 
              Shannon_Min=min(Shannon), 
              Shannon_Median=median(Shannon), 
              
              Simpson_Mean=mean(Simpson), 
              Simpson_Stdev=sd(Simpson),
              Simpson_Max=max(Simpson), 
              Simpson_Min=min(Simpson), 
              Simpson_Median=median(Simpson), 
              
              Chao1_Mean=mean(Chao1), 
              Chao1_Stdev=sd(Chao1),
              Chao1_Max=max(Chao1), 
              Chao1_Min=min(Chao1), 
              Chao1_Median=median(Chao1),
              
              N_samples = n_distinct(Index),
              N_projects = n_distinct(study_accession)) %>% as.data.frame()
  summary.df
}


run_permanova_custom <- function(my_metadata, my_formula, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  result <- adonis(my_formula,data = my_metadata, permu=permutations,method= my_method)
  # result <- adonis(my_formula,data = my_metadata, permu=999,method="bray")
  for (r in rownames(result$aov.tab)){
    variable <- r
    Degress_of_freedom <- result$aov.tab[r,]$Df[1]
    SumOfSqs <- round(result$aov.tab[r,]$SumsOfSqs[1], 3)
    meanSqs <- round(result$aov.tab[r,]$MeanSqs[1], 3)
    F.model <- round(result$aov.tab[r,]$F.Model[1], 3)
    R2 <- round(result$aov.tab[r,]$R2[1], 3)
    p_value <- round(result$aov.tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degress_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       R2,
                                                       p_value))
  }
  # my_formula_string <- paste0(as.character(my_formula)[2], as.character(my_formula)[1], as.character(my_formula)[3])
  my_formula_string <- paste0(as.character(my_formula)[1], as.character(my_formula)[3])
  print(paste0("FORMULA: ", my_formula_string))
  print(result)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Formula <- my_formula_string
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}



run_permdisp_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  dist_matrix <- vegdist(t(my_data), method = my_method)
  betadisper_object <- with(my_metadata, betadisper(dist_matrix, group = get(my_group)))
  permutest_results <- permutest(betadisper_object, permutations = permutations, parallel = 2)
  
  for (r in rownames(permutest_results$tab)){
    variable <- r
    Degrees_of_freedom <- permutest_results$tab[r,]$Df[1]
    SumOfSqs <- round(permutest_results$tab[r,]$`Sum Sq`[1],3)
    meanSqs <- round(permutest_results$tab[r,]$`Mean Sq`[1], 3)
    F.model <- round(permutest_results$tab[r,]$F[1], 3)
    N_permutations <- permutest_results$tab[r,]$N.Perm[1]
    p_value <- round(permutest_results$tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degrees_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       N_permutations,
                                                       p_value))
  }
  print(permutest_results)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","Permutations","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Group <- my_group
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}


run_anosim_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  anosim_object <- with(my_metadata, anosim(x = t(my_data), 
                                            grouping = get(my_group),
                                            permutations = permutations,
                                            distance = my_method,
                                            parallel = 2))
  
  stat_sig_table <- rbind(stat_sig_table, data.frame(my_group,
                                                     anosim_object$statistic,
                                                     anosim_object$signif,
                                                     anosim_object$permutations))
  names(stat_sig_table) <- c("Variable","R_statistic", "Significance","Permutations")
  stat_sig_table$Method <- my_method
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}
