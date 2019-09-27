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

source("Code/helper_functions.R")


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

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/other/combined_processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Load the OTU - taxonomy mapping file
# otu_taxonomy_map.df <- read.csv("Result_tables/combined/combined_otu_taxonomy_map.csv", header = T)

# Factorise discrete columns
metadata.df$Commodity <- factor(metadata.df$Commodity)
metadata.df$Sample_type <- factor(metadata.df$Sample_type)
metadata.df$Sample_treatment <- factor(metadata.df$Sample_treatment)

# Load relative abundance matrices
otu_genus_rel.m <- as.matrix(read.table(file = "Result_tables/combined/relative_abundance_tables/combined_Genus_relative_abundances.csv", sep = ",", header = T, row.names = 1))
otu_class_rel.m <- as.matrix(read.table(file = "Result_tables/combined/relative_abundance_tables/combined_Class_relative_abundances.csv", sep = ",", header = T, row.names = 1))
otu_phylum_rel.m <- as.matrix(read.table(file = "Result_tables/combined/relative_abundance_tables/combined_Phylum_relative_abundances.csv", sep = ",", header = T, row.names = 1))

# Cleanup column (sample) names - Relative abundance matrices
colnames(otu_genus_rel.m) <- gsub("_J.*", "", colnames(otu_genus_rel.m))
colnames(otu_class_rel.m) <- gsub("_J.*", "", colnames(otu_class_rel.m))
colnames(otu_phylum_rel.m) <- gsub("_J.*", "", colnames(otu_phylum_rel.m))

# And correct metadata
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

# FULL HEATMAPS
make_heatmap(otu_phylum_rel.m*100, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/phylum_relative_abundance_heatmap.pdf"),
             variables = discrete_variables,
             plot_height = 8,
             plot_width = 120,
             cluster_columns = F,
             cluster_rows = T,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             palette_choice = 'purple'
)

make_heatmap(otu_class_rel.m*100, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/class_relative_abundance_heatmap.pdf"),
             variables = discrete_variables,
             plot_height = 20,
             plot_width = 120,
             cluster_columns = F,
             cluster_rows = T,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             palette_choice = 'purple'
)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# CLASS LEVEL
class_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Class_counts_abundances_and_metadata.csv",header = T)

# Remove unknown commodities
class_data.df <- subset(class_data.df, Commodity != "Unknown")

# Generate taxonomy summary for commodity and study
class_taxa_summary.df <- generate_taxa_summary(mydata = class_data.df,taxa_column = "taxonomy_class",group_by_columns = c("Commodity", "study_accession"))

# Get top taxa by mean abundance
class_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = class_taxa_summary.df, 
                                                          grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)
write.csv(class_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_study_accession_class_top_10.csv",row.names = F, quote = F)

# Generate matrix for heatmap
heatmap.m <- class_taxa_summary.df[c("study_accession", "taxonomy_class","Mean_relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_class %in% class_taxa_summary_filtered.df$taxonomy_class,]
heatmap.m <- heatmap.m %>% spread(study_accession, Mean_relative_abundance,fill = 0)

heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- unique(metadata.df[,c("Commodity", "study_accession","Sample_type","Sample_treatment", grep("colour", names(metadata.df), value =T)), drop = F])
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$study_accession

source("Code/helper_functions.R")
make_heatmap(heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/Study_accession_class_top_10_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity","Sample_type","Sample_treatment"),
             column_title = "Study accession",
             row_title = "Class",
             plot_height = 7,
             plot_width = 10,
             cluster_columns = T,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm")
             
)

# temp2 <- log(temp*100, 2)
# temp2[is.infinite(temp2)] <- 0
# make_heatmap(temp2, 
#              mymetadata = temp_metadata.df,
#              filename = paste0("Result_figures/combined/top_10_class_per_accession_mean_relative_abundance_log_scale.pdf"),
#              variables = c("Commodity","Sample_type","Sample_treatment"),
#              column_title = "Study accession",
#              plot_height = 7,
#              plot_width = 10,
#              cluster_columns = F,
#              cluster_rows = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_name_size = 6,
#              my_annotation_palette = my_colour_palette_15,
#              #legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = round(c(0,log(c(0.001, 0.005,0.05, seq(.1,.6,.1))*100,2)),2),
#              legend_title = "Mean relative abundance %",
#              palette_choice = 'purple'
# )


class_taxa_summary.df <- generate_taxa_summary(mydata = class_data.df,taxa_column = "taxonomy_class",group_by_columns = c("Commodity"))
class_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = class_taxa_summary.df, grouping_variables = c("Commodity"),abundance_column = "Mean_relative_abundance",my_top_n = 10)
write.csv(class_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_class_top_10.csv",row.names = F, quote = F)

heatmap.m <- class_taxa_summary.df[c("Commodity", "taxonomy_class","Mean_relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_class %in% class_taxa_summary_filtered.df$taxonomy_class,]
heatmap.m <- heatmap.m %>% spread(Commodity, Mean_relative_abundance,fill = 0)

heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- unique(metadata.df[,c("Commodity","Commodity_colour"), drop = F])
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$Commodity

make_heatmap(heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/Commodity_class_top_10_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity"),
             column_title = "Commodity",
             row_title = "Class",
             plot_height = 5,
             plot_width = 7,
             cluster_columns = T,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 0,
             my_annotation_palette = my_colour_palette_15,
             #my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,1,.1))*100,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm")
)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# GENUS LEVEL

# Get the top 10 taxa per project and just show those in the plot
genus_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Genus_counts_abundances_and_metadata.csv",header = T)

# Remove unknown commodities
genus_data.df <- subset(genus_data.df, Commodity != "Unknown")

genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df,
                                               taxa_column = "taxonomy_genus",
                                               group_by_columns = c("Commodity", "study_accession"))
write.csv(genus_taxa_summary.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_study_accession_genus.csv",row.names = F, quote = F)
# genus_taxa_summary.df <- genus_taxa_summary.df[genus_taxa_summary.df$Percent_group_samples >= 10,]
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)
write.csv(genus_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_study_accession_genus_top_10.csv",row.names = F, quote = F)

heatmap.m <- genus_taxa_summary.df[c("study_accession", "taxonomy_genus","Mean_relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_genus %in% genus_taxa_summary_filtered.df$taxonomy_genus,]
heatmap.m <- heatmap.m %>% spread(study_accession, Mean_relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)

heatmap_metadata.df <- unique(metadata.df[,c("Commodity", "study_accession","Sample_type","Sample_treatment", grep("colour", names(metadata.df), value =T)), drop = F])
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$study_accession

make_heatmap(heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/Study_accession_genus_top_10_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity","Sample_type","Sample_treatment"),
             column_title = "Study accession",
             row_title = "Genus",
             plot_height = 30,
             plot_width = 24,
             cluster_columns = T,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             legend_title = "Mean relative abundance %",
             discrete_legend = T,
             palette_choice = 'purple',
             row_dend_width = unit(25, "cm")
)

# -----------
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df,taxa_column = "taxonomy_genus",group_by_columns = c("Commodity"))
write.csv(genus_taxa_summary.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_genus.csv",row.names = F, quote = F)
# genus_taxa_summary.df <- genus_taxa_summary.df[genus_taxa_summary.df$Percent_group_samples >= 10,]
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)
write.csv(genus_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_genus_top_10.csv",row.names = F, quote = F)

heatmap.m <- genus_taxa_summary.df[c("Commodity", "taxonomy_genus","Mean_relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_genus %in% genus_taxa_summary_filtered.df$taxonomy_genus,]
heatmap.m <- heatmap.m %>% spread(Commodity, Mean_relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- unique(metadata.df[,c("Commodity", "Commodity_colour"), drop = F])
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$Commodity

make_heatmap(heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/Commodity_genus_top_10_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity"),
             column_title = "Commodity",
             row_title = "Genus",
             plot_height = 12,
             plot_width = 9.5,
             cluster_columns = T,
             cluster_rows = T,
             column_title_size = 10,
             annotation_name_size = 0,
             row_title_size = 10,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm")
)

# ------------------------------------------------------------------------
# Just Genera in c__Gammaproteobacteria and c__Alphaproteobacteria
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df,taxa_column = "taxonomy_genus",group_by_columns = c("Commodity","study_accession"))
genus_taxa_summary.df <- genus_taxa_summary.df %>% filter(grepl("c__Gammaproteobacteria|c__Alphaproteobacteria", taxonomy_genus))

genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)

heatmap.m <- genus_taxa_summary.df[c("study_accession", "taxonomy_genus","Mean_relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_genus %in% genus_taxa_summary_filtered.df$taxonomy_genus,]
heatmap.m <- heatmap.m %>% spread(study_accession, Mean_relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)

heatmap_metadata.df <- unique(metadata.df[,c("Commodity", "study_accession","Sample_type","Sample_treatment", grep("colour", names(metadata.df), value =T)), drop = F])
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$study_accession

make_heatmap(heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/Commodity_study_accession_top_10_Gammaproteobacteria_and_Alphaproteobacteria_genus_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity","Sample_type","Sample_treatment"),
             column_title = "Study accession",
             row_title = "Genus",
             plot_height = 18,
             plot_width = 15,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             annotation_name_size = 10,
             row_title_size = 10,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm")
)


# ---
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df,taxa_column = "taxonomy_genus",group_by_columns = c("Commodity"))
genus_taxa_summary.df <- genus_taxa_summary.df %>% filter(grepl("c__Gammaproteobacteria|c__Alphaproteobacteria", taxonomy_genus))
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)

heatmap.m <- genus_taxa_summary.df[c("Commodity", "taxonomy_genus","Mean_relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_genus %in% genus_taxa_summary_filtered.df$taxonomy_genus,]
heatmap.m <- heatmap.m %>% spread(Commodity, Mean_relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- unique(metadata.df[,c("Commodity", "Commodity_colour"), drop = F])
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$Commodity

make_heatmap(heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/Commodity_top_10_Gammaproteobacteria_and_Alphaproteobacteria_genus_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity"),
             column_title = "Commodity",
             row_title = "Genus",
             plot_height = 12,
             plot_width = 10,
             cluster_columns = F,
             cluster_rows = F,
             column_title_size = 10,
             annotation_name_size = 0,
             row_title_size = 10,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm")
)
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# FAMILY LEVEL
family_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Family_counts_abundances_and_metadata.csv",header = T)

# Remove unknown commodities
family_data.df <- subset(family_data.df, Commodity != "Unknown")

family_taxa_summary.df <- generate_taxa_summary(mydata = family_data.df,taxa_column = "taxonomy_family",group_by_columns = c("Commodity", "study_accession"))
family_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = family_taxa_summary.df, grouping_variables = c("Commodity", "study_accession"),abundance_column = "Mean_relative_abundance",my_top_n = 10)
write.csv(family_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_study_accession_family_top_10.csv",row.names = F, quote = F)

family_taxa_summary.df <- generate_taxa_summary(mydata = family_data.df,taxa_column = "taxonomy_family",group_by_columns = c("Commodity"))
family_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = family_taxa_summary.df, grouping_variables = c("Commodity"),abundance_column = "Mean_relative_abundance",my_top_n = 10)
write.csv(family_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_family_top_10.csv",row.names = F, quote = F)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# PHYLUM LEVEL
phylum_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Phylum_counts_abundances_and_metadata.csv",header = T)

# Remove unknown commodities
phylum_data.df <- subset(phylum_data.df, Commodity != "Unknown")

phylum_taxa_summary.df <- generate_taxa_summary(mydata = phylum_data.df,taxa_column = "taxonomy_phylum",group_by_columns = c("Commodity", "study_accession"))
phylum_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = phylum_taxa_summary.df, grouping_variables = c("Commodity", "study_accession"),abundance_column = "Mean_relative_abundance",my_top_n = 10)
write.csv(phylum_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_study_accession_phylum_top_10.csv",row.names = F, quote = F)

phylum_taxa_summary.df <- generate_taxa_summary(mydata = phylum_data.df,taxa_column = "taxonomy_phylum",group_by_columns = c("Commodity"))
phylum_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = phylum_taxa_summary.df, grouping_variables = c("Commodity"),abundance_column = "Mean_relative_abundance",my_top_n = 10)
write.csv(phylum_taxa_summary_filtered.df, file = "Result_tables/combined/taxa_summary_tables/Commodity_phylum_top_10.csv",row.names = F, quote = F)
