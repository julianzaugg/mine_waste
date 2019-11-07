# Construct tree for abundant features and collapse

# See https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
#  or paper "Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses"


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

# devtools::install_github("GuangchuangYu/ggtree")
library(ggtree)


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

# Since it takes a long time to calculate, and since it was already calculated, load the sequences for most abundant features per sample across all projects
# This should be the unique set of sequences for the top 10 features by relative abundance for each sample across all projects
# In the following steps, we will filter these features further as we don't want to build a tree on all of them as
# many are going to be very low abundance, from Unknown commoditity and, if we decide to filter by region, from projects targetting different regions
seqs <- getSequences("Result_tables/combined/other/combined_most_abundant_assigned_features.fasta")
# names(seqs) <- seqs

# Load all the OTU metadata + abundance data
otu_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_OTU_counts_abundances_and_metadata.csv",header = T)

# And load the genus data. We load the genus data as we may want to filter to those features that are only in the most abundant genera.
# This may correspond to other results we have generated that are limited to the most abundant genera
genus_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_Genus_counts_abundances_and_metadata.csv",header = T)

# Remove unknown commodities
otu_data.df <- subset(otu_data.df, Commodity != "Unknown")
genus_data.df <- subset(genus_data.df, Commodity != "Unknown")

# Summarise the genus data for each study_accession
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df, taxa_column = "taxonomy_genus", group_by_columns = c("Commodity", "study_accession"))

# Get the top 10 genera for each study_accession
genus_taxa_summary_filtered.df <- filter_summary_to_top_n(taxa_summary = genus_taxa_summary.df, grouping_variables = c("Commodity", "study_accession"),
                                                          abundance_column = "Mean_relative_abundance", my_top_n = 10)

dim(otu_data.df)
# Filter the feature data to those features that are in the most abundant genera for each study_accession
otu_data_top_features.df <- subset(otu_data.df, taxonomy_genus %in% unique(genus_taxa_summary_filtered.df$taxonomy_genus))
dim(otu_data_top_features.df)

# Filter the feature data to those features that are at least 0.05% abundance
otu_data_top_features.df <- subset(otu_data_top_features.df, Relative_abundance >= 0.0005)
dim(otu_data_top_features.df)

# Filter the feature data to those most abundant features that remain
otu_data_top_features.df <- subset(otu_data_top_features.df, OTU.ID %in% names(seqs))
dim(otu_data_top_features.df)

# Summarise the remaining features in each study_accession
otu_data_top_features_taxa_summary.df <- generate_taxa_summary(mydata = otu_data_top_features.df, taxa_column = "OTU.ID", group_by_columns = c("Commodity", "study_accession"))

# Filter to those features that are in at least 10% of samples for a study_accession
otu_data_top_features_taxa_summary_filtered.df <- otu_data_top_features_taxa_summary.df %>% filter(Percent_group_samples >= 10)
otu_data_top_features.df <- subset(otu_data_top_features.df, OTU.ID %in% otu_data_top_features_taxa_summary_filtered.df$OTU.ID)

# Filter the sequences to the final list of features
seqs_filtered <- seqs[names(seqs) %in% unique(otu_data_top_features_taxa_summary_filtered.df$OTU.ID)]

# Ten most relatively abundant features for each sample were filtered to those at least 0.05% abundant and that were in at least 10% of samples 
# for a bioproject. Features that were not in any of the top ten relatively abundant genera were also removed.

# dim(otu_data_top_features_taxa_summary_filtered.df)
# length(unique(otu_data_top_features_taxa_summary_filtered.df$OTU.ID))
# length(seqs)
# length(seqs_filtered)
# summary(unique(subset(otu_data.df, OTU.ID %in% names(seqs_filtered))[c("study_accession", "Final_16S_region")]))
# temp <- unique(subset(otu_data.df, OTU.ID %in% names(seqs_filtered))[c("OTU.ID", "study_accession", "Final_16S_region")])
# length(unique(temp$OTU.ID))
# length(unique(subset(temp, Final_16S_region == "V4")$OTU.ID))
# summary(temp %>% group_by(OTU.ID) %>% summarise(Count = factor(n_distinct(Final_16S_region))) %>% arrange(desc(Count)))
# temp2 <- temp %>% group_by(OTU.ID) %>% summarise(Count = factor(n_distinct(Final_16S_region))) %>% arrange(desc(Count))
# temp3 <- temp %>% filter(Final_16S_region == "V4")
# length(unique(temp2$OTU.ID))
# length(unique(temp3$OTU.ID))

# Align the feature sequences using the DECIPHER aligner
alignment <- AlignSeqs(DNAStringSet(seqs_filtered), anchor=NA)

# Write the aligned sequences to file
writeXStringSet(alignment, file="Result_tables/combined/other/combined_most_abundant_assigned_features_aligned.fasta")


phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR_test <- update(fit, k=4, inv=0.2)
# This is will take a long time to run, e.g. ~3-4 hours for ~3500 features !
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
# fitGTR_test <- optim.pml(fitGTR, model="GTR", optInv=F, optGamma=F,
                         # rearrangement = "NNI", control = pml.control(trace = 0))

write.tree(fitGTR$tree, file = "fitted_GTR.newick")
detach("package:phangorn", unload=TRUE)

# Now that we have the tree, we want to collapse the tips to the genus level
# phyloseq has a function, tax_glom, that can do this.
# First we need to create a phyloseq object, which requires the following:
# otu_table - matrix, features can be either rows or columns (needs to be specified)
# sample_data - data.frame, rownames are the sample names in the otu_table
# tax_table - matrix, rownames must match the OTU names
# tree 

# Create the OTU table
my_otu_data.m <- subset(otu_data.df[c("OTU.ID","Sample", "Relative_abundance")], OTU.ID %in% names(seqs_filtered))
my_otu_data.m <- my_otu_data.m %>% spread(Sample,Relative_abundance,fill = 0)
my_otu_data.m <- df2matrix(my_otu_data.m)

# Sample data table
my_sample_data.df <- metadata.df[colnames(my_otu_data.m),]

# Taxonomy table
my_tax_data.m <- subset(otu_taxonomy_map.df, OTU.ID %in% names(seqs_filtered))[c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "taxonomy_phylum", "taxonomy_class", "taxonomy_order", "taxonomy_family", "taxonomy_genus", "taxonomy_species", "OTU.ID")]
rownames(my_tax_data.m) <- my_tax_data.m$OTU.ID
my_tax_data.m$OTU.ID <- NULL
my_tax_data.m <- as.matrix(my_tax_data.m)

# Create the phyloseq object  
ps <- phyloseq(otu_table(my_otu_data.m, taxa_are_rows = T),
               sample_data(my_sample_data.df),
               tax_table(my_tax_data.m),
               phy_tree(fitGTR_test$tree))
colnames(tax_table(ps))

rank_names(ps)
table(tax_table(ps)[, "taxonomy_phylum"], exclude = NULL)
table(tax_table(ps)[, "taxonomy_class"], exclude = NULL)
# table(tax_table(ps)[, "taxonomy_genus"], exclude = NULL)

# # Compute prevalence of each feature, store as data.frame
# prevdf = apply(X = otu_table(ps),
#                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
#                FUN = function(x){sum(x > 0)})
# # Add taxonomy and total read counts to this data.frame
# prevdf = data.frame(Prevalence = prevdf,
#                     TotalAbundance = taxa_sums(ps),
#                     tax_table(ps))

# How many genera would be present after filtering?
# length(get_taxa_unique(ps, taxonomic.rank = "taxonomy_genus"))

# Collapse the tree down to the genus level
ps2 <- tax_glom(ps, taxrank =  "taxonomy_genus", NArm = TRUE)


genus_tree <- phy_tree(ps2)
# length(genus_tree$tip.label) == length(unique(subset(otu_taxonomy_map.df, OTU.ID %in% genus_tree$tip.label)$taxonomy_genus))
# genus_tree$tip.label <- unlist(lapply(genus_tree$tip.label, function(x) subset(otu_taxonomy_map.df, OTU.ID == x)[,"Genus"]))

colnames(tax_table(ps2))
ggtree(genus_tree)

plot_tree(ps2,shape = NULL, ladderize = "left",label.tips = "Genus") #+ coord_polar(theta="y")
plot_tree(ps, ladderize = "left",label.tips = "Genus") #+ coord_polar(theta="y")

# Write the collapse genus tree to file
write.tree(genus_tree, file = "genus_tree.newick")
