# Processing required for SILVA based tree

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

# library("knitr")
# library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
# .inst <- .cran_packages %in% installed.packages()
# if(any(!.inst)) {
#   install.packages(.cran_packages[!.inst])
# }
# .inst <- .bioc_packages %in% installed.packages()
# if(any(!.inst)) {
#   source("http://bioconductor.org/biocLite.R")
#   biocLite(.bioc_packages[!.inst], ask = F)
# }
# 
# sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
library(phyloseq)
library(phangorn)
library(DECIPHER)
library(dada2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# devtools::install_github("GuangchuangYu/treeio")
# BiocManager::install("treeio")
library(treeio); packageVersion("treeio")

# devtools::install_github("GuangchuangYu/ggtree")
# BiocManager::install("ggtree")
library(ggtree); packageVersion("ggtree")

# install.packages("ips")
library(ips)

# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")
source("Code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/other/combined_processed_metadata.csv", sep =",", header = T)

# Remove unknown commodity samples from metadata
metadata.df <- subset(metadata.df, Commodity != "Unknown")

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/combined/other/combined_otu_taxonomy_map.csv", header = T)

# Load the genus data.
genus_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Genus_counts_abundances_and_metadata.csv",header = T)

# Remove unknown commodities
genus_data.df <- subset(genus_data.df, Commodity != "Unknown")


# Summarise the genus data for each study_accession
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df, taxa_column = "taxonomy_genus", group_by_columns = c("Commodity", "study_accession"))

# Get the top 10 genera for each study_accession
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance", my_top_n = 10)
top_10_genera.df <- melt(unique(genus_taxa_summary_filtered.df$taxonomy_genus))
names(top_10_genera.df) <- "Genus"
top_10_genera.df$Genus_silva_format <- top_10_genera.df$Genus
top_10_genera.df$Genus_silva_format <- gsub("d__", "D_0__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("p__", "D_1__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("c__", "D_2__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("o__", "D_3__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("f__", "D_4__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("g__", "D_5__", top_10_genera.df$Genus_silva_format)

# Write the lsit of top 10 genera to file
write.csv(top_10_genera.df, 
          file = "Result_tables/combined/other/combined_study_accession_top_10_genera.csv",
          row.names = F,
          quote = F)

# The above list was used to randomly select 5 representatives for each genera from the SILVA database,
# specifically using the following files:
# taxonomy/16S_only/97/majority_taxonomy_7_levels.txt
# rep_set_aligned/97/97_alignment.fna

# Ladderize tree
mytree <- read_tree("Additional_results/SILVA_extract/SILVA_genera_1_representatives_cleaned_fasttree.newick")
write.tree(ladderize(mytree),file = "Additional_results/SILVA_extract/SILVA_genera_1_representatives_cleaned_fasttree_ladderized.newick")

# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
# Now we load can load the summary table generated from this and use to 
# generate the iTOL table for tree visualisation
top_genus_silva_summary.df <- read.table("Additional_results/SILVA_extract/genera_1_rep_summary.tsv", sep = "\t", header = T)

# top_10_genera.df$Genus[!top_10_genera.df$Genus %in% top_genus_silva_summary.df$taxonomy_genus]

# We need to determine the presence absence of each genera in each Commodity, Sample treatment and Sample type
head(top_genus_silva_summary.df)

genus_data_filtered.df <- genus_data.df[genus_data.df$taxonomy_genus %in% top_genus_silva_summary.df$taxonomy_genus,]

genus_data_filtered.df <- genus_data_filtered.df[c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
                                                   "taxonomy_phylum", "taxonomy_class", "taxonomy_order", "taxonomy_family", "taxonomy_genus",
                                                   "Commodity", "Sample_type", "Sample_treatment")]
genus_data_filtered.df <- unique(genus_data_filtered.df)

# process commodity
genus_commodity.df <- df2matrix(dcast(data = genus_data_filtered.df, taxonomy_genus~Commodity,fill = 0))
for (name in colnames(genus_commodity.df)){
  assigned_colour <- as.character(subset(unique(metadata.df[c("Commodity", "Commodity_colour")]), Commodity == name)$Commodity_colour)
  genus_commodity.df[,name][genus_commodity.df[,name] > 0] <- assigned_colour
  genus_commodity.df[,name][genus_commodity.df[,name] == 0] <- "#ffffff"
}
genus_commodity.df <- m2df(genus_commodity.df, "taxonomy_genus")

# process sample type
genus_sample_type.df <- df2matrix(dcast(data = genus_data_filtered.df, taxonomy_genus~Sample_type,fill = 0))
for (name in colnames(genus_sample_type.df)){
  assigned_colour <- as.character(subset(unique(metadata.df[c("Sample_type", "Sample_type_colour")]), Sample_type == name)$Sample_type_colour)
  genus_sample_type.df[,name][genus_sample_type.df[,name] > 0] <- assigned_colour
  genus_sample_type.df[,name][genus_sample_type.df[,name] == 0] <- "#ffffff"
}
genus_sample_type.df <- m2df(genus_sample_type.df, "taxonomy_genus")

# process sample treatment
genus_sample_treatment.df <- df2matrix(dcast(data = genus_data_filtered.df, taxonomy_genus~Sample_treatment,fill = 0))
for (name in colnames(genus_sample_treatment.df)){
  assigned_colour <- as.character(subset(unique(metadata.df[c("Sample_treatment", "Sample_treatment_colour")]), Sample_treatment == name)$Sample_treatment_colour)
  genus_sample_treatment.df[,name][genus_sample_treatment.df[,name] > 0] <- assigned_colour
  genus_sample_treatment.df[,name][genus_sample_treatment.df[,name] == 0] <- "#ffffff"
}
genus_sample_treatment.df <- m2df(genus_sample_treatment.df, "taxonomy_genus")


# Merge all together
itol_data.df <- left_join(left_join(genus_commodity.df, genus_sample_type.df, by = "taxonomy_genus"), genus_sample_treatment.df, by = "taxonomy_genus")

itol_data.df <- left_join(top_genus_silva_summary.df, itol_data.df, by = "taxonomy_genus")

# Add taxonomy data (partially done in bash)
itol_data.df$taxonomy_phylum <- with(itol_data.df, paste0(Domain, Phylum))
itol_data.df$taxonomy_class <- with(itol_data.df, paste0(Domain, Phylum, Class))
itol_data.df$taxonomy_order <- with(itol_data.df, paste0(Domain, Phylum, Class, Order))
itol_data.df$taxonomy_family <- with(itol_data.df, paste0(Domain, Phylum, Class, Order, Family))

# Assign colours for each taxa level
my_colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
domain_palette <- setNames(colorRampPalette(my_colour_palette_15)(length(unique(itol_data.df$Domain))), unique(itol_data.df$Domain))
phylum_palette <- setNames(colorRampPalette(my_colour_palette_15)(length(unique(itol_data.df$taxonomy_phylum))), unique(itol_data.df$taxonomy_phylum))
class_palette <- setNames(colorRampPalette(my_colour_palette_15)(length(unique(itol_data.df$taxonomy_class))), unique(itol_data.df$taxonomy_class))
order_palette <- setNames(colorRampPalette(my_colour_palette_15)(length(unique(itol_data.df$taxonomy_order))), unique(itol_data.df$taxonomy_order))
family_palette <- setNames(colorRampPalette(my_colour_palette_15)(length(unique(itol_data.df$taxonomy_family))), unique(itol_data.df$taxonomy_family))
genus_palette <- setNames(colorRampPalette(my_colour_palette_15)(length(unique(itol_data.df$taxonomy_genus))), unique(itol_data.df$taxonomy_genus))

itol_data.df$Domain_colour <- as.character(lapply(as.character(itol_data.df$Domain), function(x) as.character(domain_palette[x])))
itol_data.df$Phylum_colour <- as.character(lapply(as.character(itol_data.df$taxonomy_phylum), function(x) as.character(phylum_palette[x])))
itol_data.df$Class_colour <- as.character(lapply(as.character(itol_data.df$taxonomy_class), function(x) as.character(class_palette[x])))
itol_data.df$Order_colour <- as.character(lapply(as.character(itol_data.df$taxonomy_order), function(x) as.character(order_palette[x])))
itol_data.df$Family_colour <- as.character(lapply(as.character(itol_data.df$taxonomy_family), function(x) as.character(family_palette[x])))
itol_data.df$Genus_colour <- as.character(lapply(as.character(itol_data.df$taxonomy_genus), function(x) as.character(genus_palette[x])))

# Reorder for ease of use
# itol_data.dfgrep("taxonomy_", names(itol_data.df), value = T)

write.csv(itol_data.df, "Result_tables/combined/other/itol_silva_metadata.csv", quote = F, row.names = F)


# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------