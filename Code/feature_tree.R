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
top_10_genera.df <- melt(unique(genus_taxa_summary_filtered.df$taxonomy_genus))
names(top_10_genera.df) <- "Genus"
top_10_genera.df$Genus_silva_format <- top_10_genera.df$Genus
top_10_genera.df$Genus_silva_format <- gsub("d__", "D_0__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("p__", "D_1__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("c__", "D_2__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("o__", "D_3__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("f__", "D_4__", top_10_genera.df$Genus_silva_format)
top_10_genera.df$Genus_silva_format <- gsub("g__", "D_5__", top_10_genera.df$Genus_silva_format)

write.csv(top_10_genera.df, 
          file = "Result_tables/combined/other/combined_study_accession_top_10_genera.csv",
          row.names = F,
          quote = F)

# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
## Generate iTOL table for tree visualisation for genera

# # TODO - make a separate script to do this for the silva entries rather than here. 
# # This is because leaf nodes will be the silva IDs rather than genera (ideally unique)
# 
# # Genus colour_for_each_commodity
# tree_summary_table.df <- subset(genus_data.df, taxonomy_genus %in% genus_taxa_summary_filtered.df$taxonomy_genus)
# 
# # Since projects are processed separately, colours are not in the abundance + metadata table. They need to be added back.
# # To do this, we need to collect the samples that each OTU.ID is found in and the corresponding Commodities etc.
# 
# tree_summary_table.df <- unique(tree_summary_table.df[c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
#                                                         "taxonomy_phylum", "taxonomy_class", "taxonomy_order", "taxonomy_family", "taxonomy_genus",
#                                                         "Commodity", "Sample_type", "Sample_treatment")])
# 
# # process commodity
# otu_commodity.df <- df2matrix(dcast(data = tree_summary_table.df, OTU.ID~Commodity,fill = 0))
# for (name in colnames(otu_commodity.df)){
#   assigned_colour <- as.character(subset(unique(metadata.df[c("Commodity", "Commodity_colour")]), Commodity == name)$Commodity_colour)
#   otu_commodity.df[,name][otu_commodity.df[,name] > 0] <- assigned_colour
#   otu_commodity.df[,name][otu_commodity.df[,name] == 0] <- "#ffffff"
# }
# 

# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------


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

summary(unique(genus_taxa_summary_filtered.df$taxonomy_genus) %in% unique(otu_data_top_features.df$taxonomy_genus))

# Filter to those features that are in at least 20% of samples for a study_accession
otu_data_top_features.df <- otu_data_top_features.df %>% filter(OTU.ID %in% unique(prevelances.df[prevelances.df$Prevalence >= 0.2,]$OTU.ID))

# A number of the most abundant genera will likely no longer be represented by the remaining features at this point
# This is primarily due to the filtering by : sequenced region, features in the top 10 per sample
summary(unique(genus_taxa_summary_filtered.df$taxonomy_genus) %in% unique(otu_data_top_features.df$taxonomy_genus))
length(unique(otu_data_top_features.df[unique(otu_data_top_features.df$taxonomy_genus) %in% unique(genus_taxa_summary_filtered.df$taxonomy_genus),]$taxonomy_family))
length(unique(otu_data_top_features.df[unique(otu_data_top_features.df$taxonomy_genus) %in% unique(genus_taxa_summary_filtered.df$taxonomy_genus),]$taxonomy_order))
length(unique(otu_data_top_features.df[unique(otu_data_top_features.df$taxonomy_genus) %in% unique(genus_taxa_summary_filtered.df$taxonomy_genus),]$taxonomy_class))
length(unique(otu_data_top_features.df[unique(otu_data_top_features.df$taxonomy_genus) %in% unique(genus_taxa_summary_filtered.df$taxonomy_genus),]$taxonomy_phylum))

write.csv(otu_data_top_features.df, "Result_tables/combined/other/combined_otu_metadata_for_tree.csv",row.names = F, quote = F)

# Filter the sequences to the final list of features
seqs_filtered <- seqs[names(seqs) %in% unique(otu_data_top_features.df$OTU.ID)]
print(paste0("Number of feature sequences remaining from top-per-sample set: ", length(seqs_filtered), "/", length(seqs)))

# Write the filtered feature sequences to file
writeXStringSet(DNAStringSet(seqs_filtered), file="Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered.fasta")

# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
## Generate iTOL table for tree visualisation for features

tree_summary_table.df <- subset(otu_data_top_features.df, OTU.ID %in% names(seqs_filtered))

# Since projects are processed separately, colours are not in the abundance + metadata table. They need to be added back.
# To do this, we need to collect the samples that each OTU.ID is found in and the corresponding Commodities etc.

tree_summary_table.df <- unique(tree_summary_table.df[c("OTU.ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", 
                                                        "taxonomy_phylum", "taxonomy_class", "taxonomy_order", "taxonomy_family", "taxonomy_genus", "taxonomy_species",
                                                        "Commodity", "Sample_type", "Sample_treatment")])

# process commodity
otu_commodity.df <- df2matrix(dcast(data = tree_summary_table.df, OTU.ID~Commodity,fill = 0))
for (name in colnames(otu_commodity.df)){
  assigned_colour <- as.character(subset(unique(metadata.df[c("Commodity", "Commodity_colour")]), Commodity == name)$Commodity_colour)
  otu_commodity.df[,name][otu_commodity.df[,name] > 0] <- assigned_colour
  otu_commodity.df[,name][otu_commodity.df[,name] == 0] <- "#ffffff"
}
otu_commodity.df <- m2df(otu_commodity.df, "OTU.ID")

# process sample type
otu_sample_type.df <- df2matrix(dcast(data = tree_summary_table.df, OTU.ID~Sample_type,fill = 0))
for (name in colnames(otu_sample_type.df)){
  assigned_colour <- as.character(subset(unique(metadata.df[c("Sample_type", "Sample_type_colour")]), Sample_type == name)$Sample_type_colour)
  otu_sample_type.df[,name][otu_sample_type.df[,name] > 0] <- assigned_colour
  otu_sample_type.df[,name][otu_sample_type.df[,name] == 0] <- "#ffffff"
}
otu_sample_type.df <- m2df(otu_sample_type.df, "OTU.ID")

# process sample treatment
otu_sample_treatment.df <- df2matrix(dcast(data = tree_summary_table.df, OTU.ID~Sample_treatment,fill = 0))
for (name in colnames(otu_sample_treatment.df)){
  assigned_colour <- as.character(subset(unique(metadata.df[c("Sample_treatment", "Sample_treatment_colour")]), Sample_treatment == name)$Sample_treatment_colour)
  otu_sample_treatment.df[,name][otu_sample_treatment.df[,name] > 0] <- assigned_colour
  otu_sample_treatment.df[,name][otu_sample_treatment.df[,name] == 0] <- "#ffffff"
}
otu_sample_treatment.df <- m2df(otu_sample_treatment.df, "OTU.ID")

# Merge all together
itol_data.df <- left_join(left_join(otu_commodity.df, otu_sample_type.df, by = "OTU.ID"), otu_sample_treatment.df, by = "OTU.ID")

# Add taxonomy data
temp <- unique(tree_summary_table.df[,!names(tree_summary_table.df) %in% c("Commodity", "Sample_treatment", "Sample_type")])
itol_data.df <- left_join(itol_data.df, temp, by = "OTU.ID")

# Create additional labels
# itol_data.df$Label <- 

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

write.csv(itol_data.df, "Result_tables/combined/other/itol_metadata.csv", quote = F, row.names = F)

# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------



# Align the feature sequences using the DECIPHER aligner
# Note - I have tried this and the alignment with default parameters was terrible (tree had long branches for some clades)
# alignment <- AlignSeqs(DNAStringSet(seqs_filtered), anchor=NA)

# Write the aligned feature sequences to file
# writeXStringSet(alignment, file="Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered_aligned.fasta")

# alignment <- DNAStringSet(readDNAMultipleAlignment("Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered_aligned.fasta"))
# alignment <- DNAStringSet(readDNAMultipleAlignment("Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered_aligned_MAFFT.fasta"))

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# Build tree with RAxML, also slow
# Requires input in different format, hence read.dna
# rax_alignment <- read.dna("Result_other/combined/sequences/combined_most_abundant_assigned_features_filtered_aligned.fasta",format="fasta",as.matrix=TRUE)
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
# phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
# dm <- dist.ml(phangAlign)
# treeNJ <- NJ(dm) # Note, tip order != sequence order
# fit = pml(treeNJ, data=phangAlign)
# fitGTR <- update(fit, k=4, inv=0.2)
# write.tree(fitGTR$tree, file = "Result_other/combined/trees/fitted_GTR.newick")
# write.tree(fitGTR$tree, file = "Result_other/combined/trees/fitted_GTR_MAFFT.newick")

# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
# write.tree(fitGTR$tree, file = "Result_other/combined/trees/fitted_GTR_optim_pml.newick")

detach("package:phangorn", unload=TRUE)

# Load alignment built externally
alignment <- DNAStringSet(readDNAStringSet("Additional_results/sina_aligner/sina_alignment_cleaned_man_filtered.fasta",format = "fasta"))

# Load tree built externally
# mytree <- read_tree("Additional_results/raxml/RAxML_bestTree.alignment")
mytree <- read_tree("Additional_results/sina_aligner/sina_tree_10col_10seq.newick")
write.tree(ladderize(mytree_relabeled),file = "Additional_results/sina_aligner/sina_tree_10col_10seq_ladderized.newick")

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
               phy_tree(mytree))

mytree_relabeled <- mytree
mytree_relabeled$tip.label <- as.character(unlist(lapply(mytree_relabeled$tip.label, function(x) subset(otu_taxonomy_map.df, OTU.ID == x)[,"Genus"])))
write.tree(ladderize(mytree_relabeled),file = "test.newick")

rank_names(ps)
# table(tax_table(ps)[, "taxonomy_phylum"], exclude = NULL)

# Collapse the tree down to the genus level
phyloseq_genus_tree <- tax_glom(ps, taxrank =  "taxonomy_genus", NArm = TRUE)

genus_tree <- phy_tree(phyloseq_genus_tree)

# genus_tree$tip.label <- as.character(unlist(lapply(genus_tree$tip.label, function(x) subset(otu_taxonomy_map.df, OTU.ID == x)[,"Genus"])))
ggtree_data.df <- otu_data_top_features.df
names(ggtree_data.df)
ggtree_data.df <- ggtree_data.df[c("OTU.ID", "Commodity")]
ggtree_data.df <- as.data.frame(+(table(ggtree_data.df)!=0)) # binarise (presence / absence)

p <- ggtree(genus_tree, layout = "circular")  + geom_tiplab(size=3, align=F, linesize=.5)
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


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
genus_tree <- phy_tree(phyloseq_genus_tree)
genus_tree$tip.label <- as.character(unlist(lapply(genus_tree$tip.label, function(x) subset(otu_taxonomy_map.df, OTU.ID == x)[,"taxonomy_genus"])))

otu_taxonomy_map.df[with(otu_taxonomy_map.df, grepl("36e33cb85", OTU.ID)),]
otu_taxonomy_map.df[with(otu_taxonomy_map.df, grepl("48c2df70d", OTU.ID)),]

otu_taxonomy_map.df[with(otu_taxonomy_map.df, grepl("4b1f67b", OTU.ID)),]



heatmap.m <- genus_taxa_summary.df[c("study_accession", "taxonomy_genus","Mean_relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_genus %in% genus_tree$tip.label ,]
heatmap.m <- heatmap.m %>% spread(study_accession, Mean_relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)
length(genus_tree$tip.label)
dim(heatmap.m)
genus_tree$tip.label %in% rownames(heatmap.m)

heatmap_metadata.df <- unique(metadata.df[,c("Commodity", "study_accession","Sample_type","Sample_treatment","Final_16S_region", 
                                             "Primers_for_16S_samples_from_manually_checking_database_or_publication",
                                             "Top_region_from_BLAST_raw_combined",
                                             grep("colour", names(metadata.df), value =T)), drop = F])
names(heatmap_metadata.df)[names(heatmap_metadata.df) == "Primers_for_16S_samples_from_manually_checking_database_or_publication"] <- "Published_16S_region"
names(heatmap_metadata.df)[names(heatmap_metadata.df) == "Top_region_from_BLAST_raw_combined"] <- "Inferred_16S_region"
heatmap_metadata.df <- subset(heatmap_metadata.df, Commodity != "Unknown")
rownames(heatmap_metadata.df) <- heatmap_metadata.df$study_accession

make_heatmap(heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("test_heatmap.pdf"),
             variables = c("Commodity","Sample_type","Sample_treatment", "Published_16S_region", "Inferred_16S_region", "Final_16S_region"),
             column_title = "Study accession",
             row_title = "Genus",
             plot_height = 30,
             plot_width = 15,
             cluster_columns = T,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             legend_title = "Mean relative abundance %",
             discrete_legend = T,
             palette_choice = 'purple',
             show_row_dend = F,
             row_dend_width = unit(25, "cm")
)


