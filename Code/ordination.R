library(vegan)
library(ggplot2)
library(ggfortify)
source("Code/helper_functions.R")

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


common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
  panel.background = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white", size = 1),
  legend.key=element_blank(),
  legend.direction="vertical",
  legend.background = element_rect(colour ="white", size = .3),
  legend.text.align = 0,
  legend.title = element_text(size=10, face="bold"),
  legend.title.align = 0.5,
  legend.margin = margin(c(2,2,2,2)),
  legend.key.height=unit(.4,"cm"),
  legend.text = element_text(size = 8),
  axis.text = element_text(size = 9, colour = "black"),
  axis.title = element_text(size = 10,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))

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

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")


# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/other/combined_processed_metadata.csv", sep =",", header = T)

# Remove unknown commodity samples
metadata.df <- subset(metadata.df, Commodity != "Unknown")

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Order the metadata
metadata.df <- metadata.df[order(rownames(metadata.df)),]

# Load count table at the OTU level. These are the counts for OTUs that were above our abundance thresholds
# otu_rare.df <- read.table("Result_tables/combined_OTU_counts_rarefied.csv", sep =",", header =T)
# otu.df <- read.table("Result_tables/combined/combined_OTU_counts.csv", sep =",", header =T)
# otu.m <- df2matrix(otu.df)
# Filter by reads per sample if you don't want to use the existing filtering
# minimum_reads <- 0
# otu.m <- otu.m[,colSums(otu.m) >= minimum_reads]
# metadata.df <- metadata.df[rownames(metadata.df) %in% colnames(otu.m),]
# otu_rare.m <- t(rrarefy(t(otu.m),sample = 5000))
# Filter by prevalance. Feature that are in at least %N percent of samples.
# TODO - features that are also in % of projects? Or % sample per project?
# prevalence_fraction <- 0.05
# summary((apply(otu.m, 1, function(x) {length(which(x > 0))}) / length(colnames(otu.m))) > prevalence_fraction)

# otu_filtered.m <- otu.m[apply(otu.m, 1, function(x) {length(which(x > 0))}) / length(colnames(otu.m)) > prevalence_fraction,]
# Filter to those features that have a maximum read count of # in at least one sample
# dim(otu_filtered.m)
# otu_filtered.m <- otu_filtered.m[apply(otu_filtered.m,1,max) >= 25,]
# otu_clr_filtered.m <- clr(otu_filtered.m)
# If there are negative values, assign them a value of zero
# otu_clr_filtered.m[which(otu_clr_filtered.m < 0)] <- 0

# Load full count table
# otu.df <- read.table("Result_tables/combined/combined_OTU_counts.csv", sep =",", header =T)
otu_genus.df <- read.table("Result_tables/combined/count_tables/combined_Genus_counts.csv", sep =",", header =T)
otu_class.df <- read.table("Result_tables/combined/count_tables/combined_Class_counts.csv", sep =",", header =T)

# Load data with metadata
# otu_data.df <- read.csv("Result_tables/combined/combined_OTU_counts_abundances_and_metadata.csv",header = T)
genus_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Genus_counts_abundances_and_metadata.csv",header = T)
class_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Class_counts_abundances_and_metadata.csv",header = T)

# Add colours (since these are not generated for individual projects)
# dim(genus_data.df)
genus_data.df <- left_join(genus_data.df, metadata.df[c("Index",grep("colour", names(metadata.df), value =T))],
          by = c("Sample" = "Index"))
class_data.df <- left_join(class_data.df, metadata.df[c("Index",grep("colour", names(metadata.df), value =T))],
                           by = c("Sample" = "Index"))
# dim(genus_data.df)
# Generate summary for each commodity + project
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df,
                                               taxa_column = "taxonomy_genus",
                                               group_by_columns = c("Commodity", "study_accession"))
class_taxa_summary.df <- generate_taxa_summary(mydata = class_data.df,
                                               taxa_column = "taxonomy_class",
                                               group_by_columns = c("Commodity", "study_accession"))

# Filter to taxa that are in at least 10 percent of samples for each group
length(unique(genus_taxa_summary.df$taxonomy_genus))
genus_taxa_summary_filtered.df <- genus_taxa_summary.df %>% filter(Percent_group_samples > 10)
length(unique(genus_taxa_summary_filtered.df$taxonomy_genus))

# length(unique(class_taxa_summary.df$taxonomy_class))
class_taxa_summary_filtered.df <- class_taxa_summary.df %>% filter(Percent_group_samples > 10)
# length(unique(class_taxa_summary_filtered.df$taxonomy_class))

# Create matrix (also filter to rows of interest)
otu_genus.m <- df2matrix(otu_genus.df)[unique(genus_taxa_summary_filtered.df$taxonomy_genus),]
otu_class.m <- df2matrix(otu_class.df)[unique(class_taxa_summary_filtered.df$taxonomy_class),]

# Order the matrices and metadata to be the same order
# otu.m <- otu.m[,rownames(metadata.df)]
otu_genus.m <- otu_genus.m[,rownames(metadata.df)]
otu_class.m <- otu_class.m[,rownames(metadata.df)]

# Filter to those entries that have at least # counts in at least one sample
# dim(otu_genus.m)
# otu_genus_filtered.m <- otu_genus.m[apply(otu_genus.m,1,max) >= 50,]
# dim(otu_genus_filtered.m)
# otu_class_filtered.m <- otu_class.m[apply(otu_class.m,1,max) >= 50,]

# CLR transform
otu_genus_clr.m <- clr(otu_genus_filtered.m)
otu_class_clr.m <- clr(otu_class_filtered.m)
# otu_genus_rare_clr.m <- clr(otu_genus_rare.m)

# Correct negative values if present
otu_genus_clr.m[which(otu_genus_clr.m < 0)] <- 0
otu_class_clr.m[which(otu_class_clr.m < 0)] <- 0
# otu_genus_rare_clr.m[which(otu_genus_rare_clr.m < 0)] <- 0

# --------------------------------------------------------------------------------
# Ordination analysis

# Sample_treatment
# Commodity
# Sample_type
# temp <- decostand(t(otu_genus_rare_clr.m), method = "hellinger")
# m.pcoa <- capscale(t(otu_genus_rare.m)~1, data = metadata.df, dist = "bray")
# m.pcoa <- capscale(t(otu_genus_clr.m)~1, data = metadata.df, dist = "bray")
# 
# generate_pca(m.pcoa, mymetadata = metadata.df,
#              plot_height = 5, plot_width = 5,
#              legend_x = -3, legend_y = 1,
#              point_size = .6, point_line_thickness = 0.3,point_alpha =.8,
#              legend_title = "Commodity",
#              legend_cex = .5,
#              plot_title = "",
#              # limits = c(-5,4,0,2),
#              plot_spiders = F,
#              plot_ellipses = F,
#              plot_hulls = F,
#              use_shapes = T,
#              ellipse_border_width = .5,
#              include_legend = T,
#              label_ellipse = F, ellipse_label_size = .5,
#              colour_palette = my_colour_palette_15,
#              variable_to_plot = "Commodity", legend_cols = 1,
#              variable_colours_available = F,
#              # my_levels = c(""),
#              filename = paste0("Result_figures/combined/combined_Commodity_pca__bray.pdf"))


# ------------------------------------------------------------------
#     TESTING
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


# pca_object <- temp
# mymetadata <- metadata.df
# variable_to_plot <- "Commodity"
# pca_site_scores <- scores(pca_object, display = "sites")
# pca_specie_scores <- scores(pca_object, display = "species")
# pca_percentages <- (pca_object$CA$eig/sum(pca_object$CA$eig)) * 100
# # Remove NA entries from the metadata and from the PCA
# internal_metadata <- mymetadata[!is.na(mymetadata[[variable_to_plot]]),]
# pca_site_scores <- pca_site_scores[rownames(pca_site_scores) %in% rownames(internal_metadata),]
# metadata_ordered.df <- internal_metadata[order(rownames(internal_metadata)),]
# metadata_ordered.df <- metadata_ordered.df[order(metadata_ordered.df[[variable_to_plot]]),]
# pca_site_scores <- pca_site_scores[rownames(metadata_ordered.df),]
# ------------------------------------------------------------------
# GENUS LEVEL

# Generate ordination object
# temp <- rda(t(otu_genus_clr.m), data = metadata.df, scale = T) # ~1 makes it unconstrained
temp <- rda(t(otu_genus_clr.m), data = metadata.df)

# ------------------------------------
# Calculate correlations between abundances for each taxa from each sample and the PC1 and PC2
# PC scores for each sample (site)
pca_site_scores <- m2df(scores(temp, display = "sites"), "Sample")

# Filter to taxonomy in pca object (because we have filtered)
genus_abundance_pc_scores.df <- subset(genus_data.df, taxonomy_genus %in% rownames(otu_genus_clr.m))

# Combine with existing metadata
genus_abundance_pc_scores.df <- left_join(genus_abundance_pc_scores.df, pca_site_scores, by = "Sample")

# Filter to columns of interest
genus_abundance_pc_scores.df <- genus_abundance_pc_scores.df[,c("Sample", "study_accession", "Commodity", "Sample_type", "Sample_treatment","taxonomy_genus","Domain","Phylum","Class","Order","Family", "Genus","Relative_abundance", "PC1", "PC2", grep("colour", names(metadata.df),value = T))]

# Filter out entries where there is no PC score
genus_abundance_pc_scores.df <- genus_abundance_pc_scores.df[!is.na(genus_abundance_pc_scores.df$PC1),]

genus_abundance_pc_scores.df$Relative_abundance <- genus_abundance_pc_scores.df$Relative_abundance*100

# Calculate the correlation between the abundances for each taxa and the PC1 and PC2 scores
genus_abundance_pc_correlations.df <- 
  genus_abundance_pc_scores.df %>% 
  group_by(taxonomy_genus) %>% summarise(Pearson_PC1 = cor(PC1, Relative_abundance, method = "pearson"),
                                         Pearson_PC2 = cor(PC2, Relative_abundance, method = "pearson"),
                                         Spearman_PC1 = cor(PC1, Relative_abundance, method = "spearman"),
                                         Spearman_PC2 = cor(PC2, Relative_abundance, method = "spearman"), 
                                         N_Samples = n_distinct(Sample),
                                         N_Projects = n_distinct(study_accession)) %>% 
  # filter(N_Samples > 5, N_Projects > 1) %>% 
  filter(N_Samples >= 5) %>% 
  arrange(desc(abs(Pearson_PC1))) %>%
  as.data.frame()

# Determine if any of the genera are in the top 10 for the commodity or a study
commodity_genus_top_10_taxa <- unique(filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)$taxonomy_genus)
commodity_study_genus_top_10_taxa <- unique(filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                              grouping_variables = c("Commodity","study_accession"),
                                                              abundance_column = "Mean_relative_abundance",
                                                              my_top_n = 10)$taxonomy_genus)
genus_abundance_pc_correlations.df$In_top_10_for_Commodity <- "no"
genus_abundance_pc_correlations.df$In_top_10_for_Commodity_accession <- "no"
genus_abundance_pc_correlations.df[genus_abundance_pc_correlations.df$taxonomy_genus %in% commodity_genus_top_10_taxa,]$In_top_10_for_Commodity <- "yes"
genus_abundance_pc_correlations.df[genus_abundance_pc_correlations.df$taxonomy_genus %in% commodity_study_genus_top_10_taxa,]$In_top_10_for_Commodity_accession <- "yes"

# Write to file
write.csv(genus_abundance_pc_correlations.df,"Result_tables/combined/ordination/Genus_relative_abundances_PC1_PC2_correlations.csv", quote = F, row.names = F)

# Plot abundances against PC1 and PC2
genus_abundance_pc_correlations_filtered.df <- genus_abundance_pc_correlations.df %>% filter(N_Samples >= 10, abs(Pearson_PC1) >= .7 | abs(Spearman_PC1) >= .7, In_top_10_for_Commodity_accession == "yes")
# genus_abundance_pc_correlations_filtered.df <- genus_abundance_pc_correlations.df %>% filter(N_Samples >= 10, In_top_10_for_Commodity_accession == "yes")

for (genus in unique(genus_abundance_pc_correlations_filtered.df$taxonomy_genus)){
  data_subset <- subset(genus_abundance_pc_scores.df, taxonomy_genus == genus)
  commodities.v <- setNames(as.character(unique(data_subset[c("Commodity","Commodity_colour")])$Commodity_colour),  as.character(unique(data_subset[c("Commodity","Commodity_colour")])$Commodity))
  myplot <- ggplot(data_subset, aes(x = PC1, y = Relative_abundance, fill = Commodity, shape = Commodity)) +
    geom_point(alpha = .8, size = 2) +
    ggtitle(genus) + 
    ylab("Relative abundance") +
    scale_fill_manual(values = commodities.v) +
    scale_shape_manual(values = rep(c(25,24,23,22,21),4)) +
    # scale_y_continuous(breaks = seq(0,100,5), limits = c(0,100)) + 
    scale_y_continuous(breaks = seq(0,100,5)) + 
    scale_x_continuous(breaks = seq(-3.5,3.5,1), limits = c(-3.5,3.5)) + 
    common_theme +
    theme(plot.title = element_text(size = 4))
  genus_splitted <- unlist(strsplit(genus,";"))
  # out_name <- paste0(gsub(" ", "_", paste(genus_splitted[2], genus_splitted[3],genus_splitted[6], sep = "___")), "_PC1.pdf")
  out_name <- paste0(gsub(";", "-", genus), "__PC1.pdf")
  print(out_name)
  ggsave(paste0("Result_figures/combined/ordination/abundances_vs_PC/", out_name), plot = myplot, width = 12, height = 12, units = "cm")
}

genus_abundance_pc_correlations_filtered.df <- genus_abundance_pc_correlations.df %>% filter(N_Samples >= 10, abs(Pearson_PC2) >= .7 | abs(Spearman_PC2) >= .7, In_top_10_for_Commodity_accession == "yes")
for (genus in unique(genus_abundance_pc_correlations_filtered.df$taxonomy_genus)){
  data_subset <- subset(genus_abundance_pc_scores.df, taxonomy_genus == genus)
  commodities.v <- setNames(as.character(unique(data_subset[c("Commodity","Commodity_colour")])$Commodity_colour),  as.character(unique(data_subset[c("Commodity","Commodity_colour")])$Commodity))
  myplot <- ggplot(data_subset, aes(x = PC2, y = Relative_abundance, fill = Commodity, shape = Commodity)) +
    geom_point(alpha = .8, size = 2) +
    ggtitle(genus) + 
    ylab("Relative abundance") +
    scale_fill_manual(values = commodities.v) +
    scale_shape_manual(values = rep(c(25,24,23,22,21),4)) +
    # scale_y_continuous(breaks = seq(0,100,5), limits = c(0,100)) + 
    scale_y_continuous(breaks = seq(0,100,5)) + 
    scale_x_continuous(breaks = seq(-3.5,3.5,1), limits = c(-3.5,3.5)) + 
    common_theme +
    theme(plot.title = element_text(size = 4))
  genus_splitted <- unlist(strsplit(genus,";"))
  # out_name <- paste0(gsub(" ", "_", paste(genus_splitted[2], genus_splitted[3],genus_splitted[6], sep = "___")), "_PC1.pdf")
  out_name <- paste0(gsub(";", "-", genus), "__PC2.pdf")
  print(out_name)
  ggsave(paste0("Result_figures/combined/ordination/abundances_vs_PC/", out_name), plot = myplot, width = 12, height = 12, units = "cm")
}

# ------------------------------------
# Generate ordination plots
generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -6, legend_y = 5,
             # legend_x = -2, legend_y = 2,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.7,
             legend_title = "Commodity",
             legend_cex = .5,
             plot_title = "",
             limits = c(-6,4,-1,2),
             # limits = c(-2,2,-2,2),
             plot_spiders = F,
             plot_ellipses = F,
             plot_hulls = F,
             use_shapes = T,
             ellipse_border_width = .5,
             include_legend = T,
             label_ellipse = F, ellipse_label_size = .3,
             colour_palette = my_colour_palette_15,
             variable_to_plot = "Commodity", legend_cols = 1,
             variable_colours_available = T,
             # my_levels = c(""),
             filename = paste0("Result_figures/combined/ordination/Commodity_genus_pca.pdf"))


generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -6, legend_y = 5,
             point_size = .6, point_line_thickness = 0.3,point_alpha =.8,
             legend_title = "Sample type",
             legend_cex = .5,
             plot_title = "",
             limits = c(-6,4,-1,2),
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
             filename = paste0("Result_figures/combined/ordination/Sample_type_genus_pca.pdf"))


generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -6, legend_y = 5,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.7,
             legend_title = "Sample treatment",
             legend_cex = .5,
             plot_title = "",
             limits = c(-6,4,-1,2),
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
             filename = paste0("Result_figures/combined/ordination/Sample_treatment_genus_pca.pdf"))



generate_pca(temp, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -6, legend_y = 5,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.7,
             legend_title = "Study accession",
             legend_cex = .4,
             plot_title = "",
             limits = c(-6,4,-1,2),
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
             filename = paste0("Result_figures/combined/ordination/Study_accession_genus_pca.pdf"))




run_permanova_custom <- function(my_metadata, my_formula){
  stat_sig_table <- NULL
  result <- adonis(my_formula,data = my_metadata, permu=999,method="euclidean")
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
  print(result)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table
}

permanova_results <- run_permanova_custom(my_metadata = metadata.df,
                                          my_formula = as.formula(t(otu_genus_clr.m)~Commodity+study_accession+Commodity:study_accession))
