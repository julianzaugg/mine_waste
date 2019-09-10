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


# ---------------------------------------------------------------------------------------------------------

# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/combined_processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Load abundance data
feature_abundances.df <- read.csv("Result_tables/combined/combined_most_abundant_unassigned.csv", header = T)
feature_abundances.df <- feature_abundances.df[feature_abundances.df$Sample %in% metadata.df$Index,]

# Load blast results
top_blast_hits.df <- read.table("Data/unassigned_blast_final.tsv", sep = "\t", header = T)
top_blast_hits.df <- top_blast_hits.df[c("query_id","subject_id","subject_scientific_name","Taxonomy")]

# species	name of a species (coincide with organism name for species-level nodes)
# genus	genus name when available
# family	family name when available
# order	order name when available
# class	class name when available
# phylum	phylum name when available
# kingdom	kingdom name when available
# superkingdom	superkingdom (domain) name when available
top_blast_hits.df <- separate(top_blast_hits.df, "Taxonomy", into = c("Super_kingdom","Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species","Strain"), remove =F, sep = ";")
# top_blast_hits.df[which(top_blast_hits.df$Species ==""),]$Species <- top_blast_hits.df[which(top_blast_hits.df$Species ==""),]$Strain
top_blast_hits.df[is.na(top_blast_hits.df)] <- "Unassigned"


top_blast_hits.df$Phylum <- as.character(lapply(top_blast_hits.df$Phylum, FUN = function(x) gsub("^", "p_", x)))

# Recreate the full taxonomy string with the 'prettier' taxa labels
top_blast_hits.df$taxonomy_species <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, Family,Genus, Species, sep =";"))
top_blast_hits.df$taxonomy_genus <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, Family,Genus, sep =";"))
top_blast_hits.df$taxonomy_family <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, Family, sep =";"))
top_blast_hits.df$taxonomy_order <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, Order, sep =";"))
top_blast_hits.df$taxonomy_class <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, Class, sep =";"))
top_blast_hits.df$taxonomy_phylum <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, Phylum, sep =";"))
top_blast_hits.df$taxonomy_kingdom <- with(top_blast_hits.df, paste(Super_kingdom, Kingdom, sep =";"))



full_table.df <- left_join(feature_abundances.df, metadata.df,by = c("Sample" = "Index"))
full_table.df <- left_join(full_table.df, top_blast_hits.df, by = c("OTU.ID" = "query_id"))

full_table.df$study_accession.y <- NULL
names(full_table.df)[names(full_table.df) == "study_accession.x"] <- "study_accession"
full_table.df <- full_table.df[!is.na(full_table.df$Taxonomy),]
# full_table.df %>% select(-instrument_platform)
# OTU.ID
# Sample
names(full_table.df)
# with(full_table.df, paste0(Super_kingdom, Kingdom, Phylum, Class, Genus))

ggplot(full_table.df, aes(x = Sample, y = Relative_abundance, fill = taxonomy_class)) +
  geom_bar(stat = "identity",show.legend = F) +
  theme(axis.text.x = element_text(angle = 90, size = 4))
  # scale_fill_manual()
