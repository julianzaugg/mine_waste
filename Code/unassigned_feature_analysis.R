# This script analyses the most abundant features that did not recieve a taxonomic assignment
# from QIIME, i.e. "Unassigned". These features have been BLASTed against the NCBI database
# and the top hit extracted from the results. The taxonomy for the top hit has also
# been extracted.


detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap) # Make Complex Heatmaps
library(circlize)  # circular visualization in R


source("Code/helper_functions.R")

# Summary function specific to the unassigned feature analysis
generate_taxa_summary_unassigned <- function(mydata, taxa_column, group_by_columns = NULL){
  select_columns <- c(taxa_column, group_by_columns, "Sample", "study_accession", "Relative_abundance")
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
                     
                     Mean_relative_abundance = round(mean(Relative_abundance), 5),
                     Median_relative_abundance = round(median(Relative_abundance), 5),
                     Min_relative_abundance = round(min(Relative_abundance),5),
                     Max_relative_abundance = round(max(Relative_abundance),5),
                     Summed_relative_abundance = round(sum(Relative_abundance),5),
    ) %>%
    as.data.frame()
  return(taxa_group_summary)
}


# ---------------------------------------------------------------------------------------------------------
# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/other/combined_processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Load abundance data for unassigned features. These abundances are calculated prior to filtering taxa or features.
# feature_abundances.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_OTU_counts_abundances_and_metadata.csv", header = T)
feature_abundances.df <- read.csv("Result_tables/combined/other/combined_most_abundant_unassigned.csv", header = T)
feature_abundances.df <- feature_abundances.df[feature_abundances.df$Sample %in% metadata.df$Index,]

# Load blast results for the most abundant features
# May need to manually fix file to remove bad characters, e.g. '#'
top_blast_hits.df <- read.table("Additional_results/unassigned_blast_final.tsv", sep = "\t", header = T) 
top_blast_hits.df <- top_blast_hits.df[c("query_id","subject_id","subject_scientific_name","Taxonomy")]

# strain
# species	name of a species (coincide with organism name for species-level nodes)
# genus	genus name when available
# family	family name when available
# order	order name when available
# class	class name when available
# phylum	phylum name when available
# kingdom	kingdom name when available
# superkingdom	superkingdom (domain) name when available
top_blast_hits.df <- separate(top_blast_hits.df, "Taxonomy", into = c("Super_kingdom","Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species","Strain"), remove =F, sep = ";")
top_blast_hits.df$Taxonomy <- NULL

# top_blast_hits.df[which(top_blast_hits.df$Species ==""),]$Species <- top_blast_hits.df[which(top_blast_hits.df$Species ==""),]$Strain
top_blast_hits.df[is.na(top_blast_hits.df)] <- "Unassigned"
top_blast_hits.df[top_blast_hits.df == ""] <- "Unassigned"


# top_blast_hits.df$Phylum <- as.character(lapply(top_blast_hits.df$Phylum, FUN = function(x) gsub("^", "p_", x)))

# Recreate the full taxonomy string with the 'prettier' taxa labels
top_blast_hits.df$taxonomy_kingdom <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, sep =";"))
top_blast_hits.df$taxonomy_phylum <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, sep =";"))
top_blast_hits.df$taxonomy_class <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, sep =";"))
top_blast_hits.df$taxonomy_order <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, sep =";"))
top_blast_hits.df$taxonomy_family <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, Family, sep =";"))
top_blast_hits.df$taxonomy_genus <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, Family,Genus, sep =";"))
top_blast_hits.df$taxonomy_species <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, Family,Genus, Species, sep =";"))
top_blast_hits.df$taxonomy_strain <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, Family,Genus, Species,Strain, sep =";"))


full_table.df <- left_join(feature_abundances.df, metadata.df,by = c("Sample" = "Index"))
full_table.df <- left_join(full_table.df, top_blast_hits.df, by = c("OTU.ID" = "query_id"))

full_table.df$study_accession.y <- NULL
names(full_table.df)[names(full_table.df) == "study_accession.x"] <- "study_accession"
full_table.df <- full_table.df[!is.na(full_table.df$taxonomy_strain),]
full_table.df <- full_table.df[,c(grep("RepSeq", names(full_table.df), value = T, invert = T), "RepSeq")]
write.csv(full_table.df, file = "Result_tables/combined/other/unassigned_blast_results_processed.csv", quote = F, row.names = F)

# Remove unknown commodity samples
full_table.df <- subset(full_table.df, Commodity != "Unknown")


# Generate taxa summaries
strain_taxa_summary.df <- generate_taxa_summary_unassigned(mydata = full_table.df,taxa_column = "taxonomy_strain",group_by_columns = c("Commodity", "study_accession"))
species_taxa_summary.df <- generate_taxa_summary_unassigned(mydata = full_table.df,taxa_column = "taxonomy_species",group_by_columns = c("Commodity", "study_accession"))
genus_taxa_summary.df <- generate_taxa_summary_unassigned(mydata = full_table.df,taxa_column = "taxonomy_genus",group_by_columns = c("Commodity", "study_accession"))


strain_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = strain_taxa_summary.df, 
                                                          grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)

species_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = species_taxa_summary.df, 
                                                          grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, 
                                                          grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance",
                                                          my_top_n = 10)


heatmap_strain.m <- strain_taxa_summary.df[c("study_accession", "taxonomy_strain","Mean_relative_abundance")]
# heatmap_strain.m <- heatmap_strain.m[heatmap_strain.m$taxonomy_strain %in% strain_taxa_summary_filtered.df$taxonomy_strain,]
heatmap_strain.m <- heatmap_strain.m %>% spread(study_accession, Mean_relative_abundance,fill = 0)
heatmap_strain.m <- df2matrix(heatmap_strain.m)

heatmap_genus.m <- genus_taxa_summary.df[c("study_accession", "taxonomy_genus","Mean_relative_abundance")]
heatmap_genus.m <- heatmap_genus.m[heatmap_genus.m$taxonomy_genus %in% genus_taxa_summary_filtered.df$taxonomy_genus,]
heatmap_genus.m <- heatmap_genus.m %>% spread(study_accession, Mean_relative_abundance,fill = 0)
heatmap_genus.m <- df2matrix(heatmap_genus.m)

heatmap_metadata.df <- unique(metadata.df[,c("Commodity", "study_accession","Sample_type","Sample_treatment", grep("colour", names(metadata.df), value =T)), drop = F])
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$study_accession


make_heatmap(heatmap_genus.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/combined/heatmaps/Study_accession_unassigned_genus_top_10_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity","Sample_type","Sample_treatment"),
             column_title = "Study accession",
             row_title = "Genus",
             plot_height = 15,
             plot_width = 12,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 10,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             legend_title = "Mean relative abundance %",
             discrete_legend = T,
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             # row_dend_width = unit(5, "cm")
)


genus_taxa_summary.df <- generate_taxa_summary_unassigned(mydata = full_table.df,taxa_column = "taxonomy_genus",group_by_columns = c("Commodity"))
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
             filename = paste0("Result_figures/combined/heatmaps/Commodity_unassigned_genus_top_10_mean_relative_abundance_heatmap.pdf"),
             variables = c("Commodity"),
             column_title = "Commodity",
             row_title = "Genus",
             plot_height = 10,
             plot_width = 9,
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

