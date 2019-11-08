# Construct tree for abundant features and collapse

# See https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
#  or paper "Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses"

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

# Since it takes a long time to calculate, and since it was already calculated, load the sequences for most abundant features per sample across all projects
# This should be the unique set of sequences for the top 10 features by relative abundance for each sample across all projects
# In the following steps, we will filter these features further as we don't want to build a tree on all of them as
# many are going to be very low abundance, from Unknown commoditity and, if we decide to filter by region, from projects targetting different regions
seqs <- getSequences("Result_other/combined/sequences/combined_most_abundant_assigned_features.fasta")
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

# Filter the feature data to those features that are at least 0.1% abundance
otu_data_top_features.df <- subset(otu_data_top_features.df, Relative_abundance >= 0.001)
dim(otu_data_top_features.df)

summary(unique(genus_taxa_summary_filtered.df$taxonomy_genus) %in% unique(otu_data_top_features.df$taxonomy_genus))
# Filter the feature data to those most abundant features that remain
# otu_data_top_features.df <- subset(otu_data_top_features.df, OTU.ID %in% names(seqs))
# dim(otu_data_top_features.df)

# (optional) Filter to features that are only in projects that targetted just the V4 region
# otu_data_top_features.df <- subset(otu_data_top_features.df, Final_16S_region == "V4")

# Calculate the prevalence of the remaining features in the full data (prior to filtering)
N_samples_per_project <- otu_data.df %>% group_by(study_accession) %>% summarise(N_samples = n_distinct(Sample))
N_samples_per_feature <- otu_data_top_features.df %>% group_by(study_accession, OTU.ID) %>% summarise(In_N_samples = n_distinct(Sample))
# head(N_samples_per_feature)
# head(N_samples_per_project)
prevelances.df <- left_join(N_samples_per_feature, N_samples_per_project, by = "study_accession")
prevelances.df$Prevalence <- with(prevelances.df, In_N_samples/N_samples)
# subset(otu_data.df, OTU.ID == "9016f374255e870578c2fbb416ac42e6")
# unique(subset(otu_data.df, study_accession == "PRJNA339895")$Sample)
# subset(prevelances.df, OTU.ID == "9016f374255e870578c2fbb416ac42e6")

# Filter to those features that are in at least 20% of samples for a study_accession
otu_data_top_features.df <- otu_data_top_features.df %>% filter(OTU.ID %in% unique(prevelances.df[prevelances.df$Prevalence >= 0.2,]$OTU.ID))

# A number of the most abundant genera will likely no longer be represented by the remaining features at this point
# This is primarily due to the filtering by sequenced region, to those features in the top 10 per sample,
summary(unique(genus_taxa_summary_filtered.df$taxonomy_genus) %in% unique(otu_data_top_features.df$taxonomy_genus))

write.csv(otu_data_top_features.df, "Result_tables/combined/other/combined_otu_metadata_for_tree.csv",row.names = F, quote = F)

# Filter the sequences to the final list of features
seqs_filtered <- seqs[names(seqs) %in% unique(otu_data_top_features.df$OTU.ID)]
print(paste0("Number of feature sequences remaining from top-per-sample set: ", length(seqs_filtered), "/", length(seqs)))

# The ten most relatively abundant features for each sample were filtered to those that were at least 0.1% abundant in at least
# one sample and that were in at least 20% of samples for any bioproject. Features that were not in any of the top ten relatively abundant genera were also removed.

# Write the filtered feature sequences to file
writeXStringSet(DNAStringSet(seqs_filtered), file="Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered.fasta")

# Align the feature sequences using the DECIPHER aligner
alignment <- AlignSeqs(DNAStringSet(seqs_filtered), anchor=NA)

# Write the aligned feature sequences to file
writeXStringSet(alignment, file="Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered_aligned.fasta")

alignment <- DNAStringSet(readDNAMultipleAlignment("Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered_aligned.fasta"))

# Align with RAxML, also slow
# Requires input in different format, hence read.dna
rax_alignment <- read.dna("Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered_aligned.fasta",format="fasta",as.matrix=TRUE)
# alignment.rax.gtr <- raxml(rax_alignment,
#                            m="GTRGAMMAIX", # model
#                            f="a", # best tree and bootstrap
#                            p=1234, # random number seed
#                            x=2345, # random seed for rapid bootstrapping
#                            N=100, # number of bootstrap replicates
#                            file="alignment", # name of output files
#                            #exec="raxmlHPC-PTHREADS-SSE3", # name of executable
#                            exec = "/Applications/miniconda3/envs/raxml_8.2.12/bin/raxmlHPC-PTHREADS-SSE3",
#                            threads=2
# )

# Align with phangorn (can be slower) if optimising with optim.pml, e.g. ~3-4 hours for ~3500 features !
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
write.tree(fitGTR$tree, file = "Result_other/combined/trees/fitted_GTR.newick")

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

write.tree(fitGTR$tree, file = "Result_other/combined/trees/fitted_GTR_optim_pml.newick")

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
               phy_tree(fitGTR$tree))

rank_names(ps)
# table(tax_table(ps)[, "taxonomy_phylum"], exclude = NULL)
phy_tree(fitGTR$tree)

# Collapse the tree down to the genus level
ps2 <- tax_glom(ps, taxrank =  "taxonomy_genus", NArm = TRUE)

genus_tree <- phy_tree(ps2)

genus_tree$tip.label <- as.character(unlist(lapply(genus_tree$tip.label, function(x) subset(otu_taxonomy_map.df, OTU.ID == x)[,"Genus"])))
ggtree_data.df <- otu_data_top_features.df
names(ggtree_data.df)
ggtree_data.df <- ggtree_data.df[c("OTU.ID", "Commodity")]
ggtree_data.df <- as.data.frame(+(table(ggtree_data.df)!=0)) # binarise (presence / absence)

p <- ggtree(genus_tree, layout = "rectangular")  + geom_tiplab(size=3, align=F, linesize=.5)
p
gheatmap(p, data = ggtree_data.df)

ggtree_data.df <- dcast(ggtree_data.df, formula =OTU.ID~ Commodity)
ggtree_data.df <- df2matrix(ggtree_data.df)
ggtree_data.df[ggtree_data.df > 0] <- 1



ggtree_data.df$run_accession <- NULL
rownames(ggtree_data.df) <- ggtree_data.df$OTU.ID

# layout one of 'rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle' or 'daylight'
p <- ggtree(genus_tree, layout = "rectangular", branch.length='rate') 
gheatmap(p, data = )
?gheatmap

# seqs["27975adba200137bab8ad346917aee84"]
# length(genus_tree$tip.label) == length(unique(subset(otu_taxonomy_map.df, OTU.ID %in% genus_tree$tip.label)$taxonomy_genus))
# genus_tree$tip.label <- unlist(lapply(genus_tree$tip.label, function(x) subset(otu_taxonomy_map.df, OTU.ID == x)[,"Genus"]))

# Create metadata table for tree. Need columns (tracks) for each variable value. need 

x <- data.frame(label = genus_tree$tip.label,as.data.frame(my_tax_data.m)[genus_tree$tip.label,])

## convert the phylo object to a treeio::treedata object
genus_tree <- treedata(phylo = genus_tree)

## add the annotation
genus_tree <- full_join(genus_tree, x, by="label")

ggtree(ps2) +geom_text(aes(x=branch, label=Class, color = Phylum)) + coord_polar(theta="y") #geom_tiplab()

plot_tree(ps2,shape = NULL, ladderize = "left",label.tips = "Genus") #+ coord_polar(theta="y")
plot_tree(ps, ladderize = "left",label.tips = "Genus") #+ coord_polar(theta="y")

# Write the collapse genus tree to file
write.tree(phy_tree(ps2), file = "Result_other/combined/trees/genus_tree.newick")
