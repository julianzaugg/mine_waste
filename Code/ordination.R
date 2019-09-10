library(vegan)
library(ggplot2)
library(ggfortify)

############################################################
# Various colour colour_palettes
my_colour_palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
my_colour_palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
my_colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
my_colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
my_colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
my_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
my_colour_palette_10_soft <- c("#9E788F","#4C5B61","#678D58","#AD5233","#A0A083","#4D456A","#588578","#D0AC4C","#2A7BA0","#931621")
my_colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
############################################################


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

# For each rowname (OTU), get the corresponding taxonomy_species
# Assumes "OTU.ID" and "taxonomy_species" columns in the provided map dataframe
assign_taxonomy_to_otu <- function(otutable, taxon_map){
  taxonomies <- c()
  for (otuid in rownames(otutable)){
    taxonomies <- c(taxonomies, as.character(taxon_map[taxon_map$OTU.ID == otuid,]$taxonomy_species))
  }
  return(taxonomies)
}

# Function that takes the metadata and a list of variables (column names) and returns those samples (rownames) with NA entries
get_samples_missing_data <- function(my_metadata, variables){
  samples_missing_data <- c()
  for (name in variables) {
    samples_missing_data <- c(samples_missing_data, rownames(my_metadata[is.na(my_metadata[[name]]),]))
  }
  return(unique(samples_missing_data))
}

df2matrix <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")


# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/combined_processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Load count table at the OTU level. These are the counts for OTUs that were above our abundance thresholds
# otu_rare.df <- read.table("Result_tables/combined_OTU_counts_rarefied.csv", sep =",", header =T)
otu.df <- read.table("Result_tables/combined/combined_OTU_counts.csv", sep =",", header =T)
otu_genus.df <- read.table("Result_tables/combined/combined_Genus_counts.csv", sep =",", header =T)

# Create matrices
# otu_rare.m <- otu_rare.df
# rownames(otu_rare.m) <- otu_rare.df$OTU.ID
# otu_rare.m$OTU.ID <- NULL
# otu_rare.m <- as.matrix(otu_rare.m)

otu.m <- df2matrix(otu.df)
otu_genus.m <- df2matrix(otu_genus.df)

# Filter by reads per sample if you don't want to use the existing filtering
minimum_reads <- 0
# temp <- otu_rare.m[,colSums(otu_rare.m) == minimum_reads]
# otu_rare.m <- otu_rare.m[,colSums(otu_rare.m) >= minimum_reads]
otu.m <- otu.m[,colSums(otu.m) >= minimum_reads]
# metadata.df <- metadata.df[rownames(metadata.df) %in% colnames(otu_rare.m),]
metadata.df <- metadata.df[rownames(metadata.df) %in% colnames(otu.m),]


# Order the matrices and metadata to be the same order
metadata.df <- metadata.df[order(rownames(metadata.df)),]
# otu_rare.m <- otu_rare.m[,order(rownames(metadata.df))]
otu.m <- otu.m[,order(rownames(metadata.df))]
otu_genus.m <- otu_genus.m[,order(rownames(metadata.df))]

# summary((apply(otu_rare.m, 1, function(x) {length(which(x > 0))}) / length(colnames(otu_rare.m))) > 0.01)
# otu_rare.m <- otu_rare.m[apply(otu_rare.m, 1, function(x) {length(which(x > 0))}) / length(colnames(otu_rare.m)) > 0.01,]
# otu_rare_filtered.m <- otu_rare.m[apply(otu_rare.m,1,max) >= 25,]
# otu_rare_clr_filtered.m <- clr(otu_rare_filtered.m)

# Filter by prevalance. Feature that are in at least %N percent of samples.
# TODO - features that are also in % of projects? Or % sample per project?
prevalence_fraction <- 0.05
summary((apply(otu.m, 1, function(x) {length(which(x > 0))}) / length(colnames(otu.m))) > prevalence_fraction)

otu_filtered.m <- otu.m[apply(otu.m, 1, function(x) {length(which(x > 0))}) / length(colnames(otu.m)) > prevalence_fraction,]
# Filter to those features that have a maximum read count of # in at least one sample
dim(otu_filtered.m)
otu_filtered.m <- otu_filtered.m[apply(otu_filtered.m,1,max) >= 25,]
otu_clr_filtered.m <- clr(otu_filtered.m)

otu_genus_filtered.m <- otu_genus.m[apply(otu_genus.m, 1, function(x) {length(which(x > 0))}) / length(colnames(otu_genus.m)) > prevalence_fraction,]
otu_genus_filtered.m <- otu_genus_filtered.m[apply(otu_genus_filtered.m,1,max) >= 25,]
otu_genus_clr.m <- clr(otu_genus_filtered.m)

# If there are negative values, assign them a value of zero
# otu_rare_clr_filtered.m[which(otu_rare_clr_filtered.m < 0)] <- 0
otu_clr_filtered.m[which(otu_clr_filtered.m < 0)] <- 0
# otu_rare_clr_rel_filtered.m[which(otu_rare_clr_rel_filtered.m < 0)] <- 0
metadata_filtered.df <- metadata.df[rownames(metadata.df) %in% colnames(otu_clr_filtered.m),]
# rownames(metadata_filtered.df) == colnames(otu_clr_filtered.m)

# --------------------------------------------------------------------------------
# Ordination analysis
generate_pca <- function(pca_object, mymetadata, variable_to_plot, colour_palette, limits = NULL, filename = NULL, include_legend = T, add_spider = F, add_ellipse = F,
                         point_alpha = 1, plot_width = 10, plot_height=10, point_size = 0.8, point_line_thickness = 1,
                         label_sites = F, label_species = F,
                         legend_x = NULL, legend_y = NULL, legend_x_offset = 0, legend_y_offset = 0,
                         legend_cex = 0.6,
                         plot_spiders = NULL, plot_ellipses = NULL,plot_hulls = NULL, legend_cols = 2, legend_title = NULL,
                         label_ellipse = F, ellipse_label_size = 0.5, ellipse_border_width = 1,variable_colours_available = F, 
                         plot_title = NULL, use_shapes = F, my_levels = NULL){
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
  
  plot(pca_object,
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
    title(main = plot_title)
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

# Sample_treatment
# Commodity
# Sample_type

# samples <- rownames(subset(metadata.df, Commodity == "Arsenic"))
# temp <- rda(t(otu_genus_clr.m[,samples]), data = metadata.df, scale = T) # ~1 makes it unconstrained
temp <- rda(t(otu_genus_clr.m), data = metadata.df, scale = T) # ~1 makes it unconstrained

# temp <- rda(t(clr(otu_genus.m)), scale = T)
# temp <- rda(t(otu_genus_clr.m), data = metadata.df) # ~1 makes it unconstrained
# temp <- rda(t(otu_clr_filtered.m), data = metadata.df) # ~1 makes it unconstrained
# temp <- rda(t(otu_clr_filtered.m), data = metadata.df) # ~1 makes it unconstrained
# temp <- prcomp(t(otu_clr_filtered.m), center = T, scale = T)
# capscale(t(otu_filtered.m)~1, distance = "bray")
# plot(temp)
# plot(temp, scaling = 3)
# pca.res <- prcomp(t(otu_clr_filtered.m), center = T, scale. = T)
# temp <- rda(t(otu_clr_filtered.m)~Commodity, data = metadata.df) # ~1 makes it unconstrained
# temp <- capscale(t(otu_clr_filtered.m) ~ Commodity, data = metadata_filtered.df)
# temp <- metaMDS(t(otu_genus_clr.m),distance = "euclidean",trace=0)

generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -5, legend_y = 4,
             point_size = .6, point_line_thickness = 0.3,point_alpha =.8,
             legend_title = "Sample type",
             legend_cex = .5,
             plot_title = "",
             limits = c(-5,4,0,2),
             plot_spiders = F,
             plot_ellipses = F,
             plot_hulls = F,
             use_shapes = T,
             ellipse_border_width = .5,
             include_legend = T,
             label_ellipse = F, ellipse_label_size = .5,
             colour_palette = rev(my_colour_palette_15),
             variable_to_plot = "Sample_type", legend_cols = 1,
             variable_colours_available = F,
             # my_levels = c(""),
             filename = paste0("Result_figures/combined/combined_Sample_type_pca.pdf"))

generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -5, legend_y = 4,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.7,
             legend_title = "Commodity",
             legend_cex = .5,
             plot_title = "",
             limits = c(-5,4,0,2),
             plot_spiders = F,
             plot_ellipses = F,
             plot_hulls = F,
             use_shapes = T,
             ellipse_border_width = .5,
             include_legend = T,
             label_ellipse = F, ellipse_label_size = .5,
             colour_palette = my_colour_palette_15,
             variable_to_plot = "Commodity", legend_cols = 1,
             variable_colours_available = F,
             # my_levels = c(""),
             filename = paste0("Result_figures/combined/combined_Commodity_pca.pdf"))



generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -5, legend_y = 4,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.7,
             legend_title = "Sample treatment",
             legend_cex = .5,
             plot_title = "",
             limits = c(-5,4,0,2),
             plot_spiders = F,
             plot_ellipses = F,
             plot_hulls = F,
             use_shapes = T,
             ellipse_border_width = .5,
             include_legend = T,
             label_ellipse = F, ellipse_label_size = .5,
             colour_palette = my_colour_palette_15,
             variable_to_plot = "Sample_treatment", legend_cols = 1,
             variable_colours_available = F,
             # my_levels = c(""),
             filename = paste0("Result_figures/combined/combined_Sample_treatment_pca.pdf"))



generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -5, legend_y = 4,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.7,
             legend_title = "Study accession",
             legend_cex = .4,
             plot_title = "",
             limits = c(-5,4,0,2),
             plot_spiders = F,
             plot_ellipses = F,
             plot_hulls = F,
             use_shapes = T,
             ellipse_border_width = .5,
             include_legend = T,
             label_ellipse = F, ellipse_label_size = .5,
             colour_palette = my_colour_palette_206_distinct,
             variable_to_plot = "study_accession", legend_cols = 2,
             variable_colours_available = F,
             # my_levels = c(""),
             filename = paste0("Result_figures/combined/combined_study_accession_pca.pdf"))


