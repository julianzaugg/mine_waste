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



# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")
source("Code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/other/combined_processed_metadata.csv", sep =",", header = T)

# Remove unknown commodity samples
metadata.df <- subset(metadata.df, Commodity != "Unknown")

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Since it takes a long time to calculate now, and since it was already generate, load the most abundant features per sample across all projects
# This should be the unique set of sequences for the top 10 features by relative abundance for each sample across all projects
seqs <- getSequences("Result_tables/combined/other/combined_most_abundant_assigned_features.fasta")
names(seqs) <- seqs

# Building a tree from all these features is unnecessary. Many are going to be very low abundance.
# Load all the OTU metadata + abundance data
otu_data.df <- read.csv("Result_tables/combined/combined_counts_abundances_and_metadata_tables/combined_OTU_counts_abundances_and_metadata.csv",header = T)

# Remove unknown commodities
otu_data.df <- subset(otu_data.df, Commodity != "Unknown")

# Filter the data to those features that are at least 0.05% abundance
otu_data_top_features.df <- subset(otu_data.df, Relative_abundance >= 0.0005)

# Filter the data to those most abundant features that remain
otu_data_top_features.df <- subset(otu_data_top_features.df, OTU.ID %in% names(seqs))

# Summarise the features in each study_accession
otu_data_top_features_taxa_summary.df <- generate_taxa_summary(mydata = otu_data_top_features.df, taxa_column = "OTU.ID", group_by_columns = c("Commodity", "study_accession"))

# Filter to those features that are in at least 10% of samples for a study_accession
otu_data_top_features_taxa_summary_filtered.df <- otu_data_top_features_taxa_summary.df %>% filter(Percent_group_samples >= 10)

# Filter the sequences to the final list of features
seqs_filtered <- seqs[names(seqs) %in% unique(otu_data_top_features_taxa_summary_filtered.df$OTU.ID)]

# Ten most relatively abundant features for each sample were filtered to those at least 0.05% abundant and that were in at least 10% of samples 
# for a bioproject

dim(otu_data_top_features_taxa_summary_filtered.df)
length(unique(otu_data_top_features_taxa_summary_filtered.df$OTU.ID))
length(seqs)
length(seqs_filtered)
# At this point we have 


alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# Write the aligned sequences to file
writeXStringSet(alignment, file="Result_tables/combined/other/combined_most_abundant_assigned_features_aligned.fasta")


phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
# This is will take a long time to run, e.g. ~3-4 hours for ~3500 features !
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
write.tree(fitGTR$tree, file = "test.tree")
detach("package:phangorn", unload=TRUE)