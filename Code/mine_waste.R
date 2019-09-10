#************************************
# This script is the main data preparation script. Data is formatted for use elsewhere
# and relative abundances calculated at different taxonomy levels.
#
# NOTE - Although the below script refers to each representative sequence as an OTU, in reality
# these sequences are what you would call an amplicon sequence variant (ASV).
# For more information, see : https://www.nature.com/articles/ismej2017119
#************************************

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

# Uncomment and run to install the libraries that might be needed 
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# install.packages("ggplot2")
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("RColorBrewer")
# install.packages("vegan")
# install.packages("reshape2")
# install.packages("gplots")

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(gplots)

library(seqinr) # For writing fasta files


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


# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")


###############################################################
# Create result directories if they are missing
dir.create(file.path(".", "Result_figures"), showWarnings = FALSE)
dir.create(file.path(".", "Result_tables"), showWarnings = FALSE)
dir.create(file.path(".", "Result_objects"), showWarnings = FALSE)

dir.create(file.path("./Result_tables","combined"), showWarnings = FALSE)
dir.create(file.path("./Result_figures","combined"), showWarnings = FALSE)


create_project_result_dirs <- function(project_name){
  
  dir.create(file.path("./Result_figures",project_name, "abundance_analysis"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_figures",project_name, "DESeq"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_figures",project_name, "diversity_analysis"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_figures",project_name, "exploratory_analysis"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_figures",project_name, "heatmaps"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_figures",project_name, "ordination"), showWarnings = FALSE,recursive = T)
  
  dir.create(file.path("./Result_tables",project_name, "abundance_analysis_tables"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_tables",project_name, "count_tables"), showWarnings = FALSE,recursive = T)
  
  dir.create(file.path("./Result_tables",project_name, "DESeq_results"), showWarnings = FALSE,recursive = T)
  
  dir.create(file.path("./Result_tables",project_name, "diversity_analysis"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_tables",project_name, "other"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_tables",project_name, "relative_abundance_tables"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_tables",project_name, "stats_various"), showWarnings = FALSE,recursive = T)
  dir.create(file.path("./Result_tables",project_name, "contaminant_analysis"), showWarnings = FALSE,recursive = T)
}



# ------------------------------------------------------------
# ------------------------------------------------------------

# Make row names of a dataframe the value of a defined column (default 1st) and then remove the column
clean_dataframe <- function(mydf, rowname_col = 1){
  my_clean.df <- mydf
  rownames(my_clean.df) <- my_clean.df[,rowname_col]
  my_clean.df[,1] <- NULL
  return(my_clean.df)
}

create_combined_dataframe <- function(counts.df, counts_rare.df, abundances.df, abundances_rare.df, mymetadata, mylevel = "OTU", otu_map.df = NULL){
  counts <- clean_dataframe(counts.df)
  rel_abundances <- clean_dataframe(abundances.df)
  counts_rare <- clean_dataframe(counts_rare.df)
  rel_abundances_rare <- clean_dataframe(abundances_rare.df)
  
  # Ensure ordering is the same
  rel_abundances <- rel_abundances[rownames(counts),,drop = F]
  counts_rare <- counts_rare[rownames(counts),,drop = F]
  rel_abundances_rare <- rel_abundances_rare[rownames(counts),,drop = F]
  
  # Combine the datasets. Passing as.matrix(counts) captures the rownames as a column. This can be renamed after
  combined_data <- cbind(melt(as.matrix(counts), variable.name = "sample", value.name = "Read_count"),
                         melt(rel_abundances, value.name = "Relative_abundance")[,2, drop = F],
                         melt(counts_rare, value.name = "Read_count_rarefied")[,2, drop = F],
                         melt(rel_abundances_rare, value.name = "Relative_abundance_rarefied")[,2, drop = F])
  
  # Remove samples with a read count of zero
  combined_data <- combined_data[combined_data$Read_count > 0,]
  
  # Calculate logged read counts
  combined_data$Read_count_logged <- log(combined_data$Read_count, 10)
  combined_data$Read_count_rarefied_logged <- log(combined_data$Read_count_rarefied, 10)
  
  # Fix the Var2 column
  names(combined_data)[2] <- "Sample"
  
  # Merge with metadata. Assumes an Index column matching Sample
  combined_data <- merge(combined_data, mymetadata, by.x = "Sample", by.y = "Index")
  
  if (mylevel == "OTU.ID"){
    names(combined_data)[names(combined_data) == "Var1"] <- "OTU.ID"
    combined_data <- merge(combined_data, otu_map.df, by.x = "OTU.ID", by.y = "OTU.ID")
  }
  else if (mylevel == "Species"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_species"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species", "taxonomy_species")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_species", by.y = "taxonomy_species")
  }
  else if (mylevel == "Genus"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_genus"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "Genus", "taxonomy_genus")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_genus", by.y = "taxonomy_genus")
  }
  else if (mylevel == "Family"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_family"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "taxonomy_family")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_family", by.y = "taxonomy_family")
  }
  else if (mylevel == "Order"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_order"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "taxonomy_order")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_order", by.y = "taxonomy_order")
  }
  else if (mylevel == "Class"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_class"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "taxonomy_class")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_class", by.y = "taxonomy_class")
  }
  else if (mylevel == "Phylum"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_phylum"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "taxonomy_phylum")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_phylum", by.y = "taxonomy_phylum")
  }
  return(combined_data)
}


# Taxonomy-sample matrix to dataframe convertor
# Just converts a matrix where the row name is a taxonomy label (really can be anything)
m2df <- function(mymatrix, name_of_taxonomy_col = "taxonomy"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, name_of_taxonomy_col] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(name_of_taxonomy_col,cur_names)]
  return(mydf)
}
# ------------------------------------------------------------
# ------------------------------------------------------------

# metaA <- read.table("data/metadata_A.tsv", header = T, sep = "\t")
# metaB <- read.table("data/metadata_B.tsv", header = T, sep = "\t")

# dim(metaA)
# dim(metaB)
# metadata.df <- left_join(metaA, metaB,by = c("study_accession" = "BioProject_Accession"))

# temp <- read.table(pipe("pbpaste"), sep = "\t", header = T)
# metadata.df <- left_join(metadata.df, temp,by = c("study_accession" = "Bioproject"))
# metadata.df <- left_join(metadata.df, temp,by = "run_accession")

# dim(metadata.df)
# write.csv(x = metadata.df, file = "data/metadata.csv", quote = F, row.names = F)
metadata.df <- read.table("data/mine_waste_metadata.tsv", header = T, sep = "\t")

# --------------------------------------------------------------------
# Note, the following projects are the same, either can be removed
# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA493908/
# https://www.ncbi.nlm.nih.gov/bioproject/PRJEB28611/
metadata.df <- metadata.df[metadata.df$study_accession != "PRJNA493908",]

# We are only interested in a samples that were targetted at the V4 region
metadata.df <- subset(metadata.df, Top_region_from_BLAST == "V4")

# Remove samples where it was explicit that the ITS region was sequenced
metadata.df <- subset(metadata.df, region_ITS_targeted == "no")

# Remove samples where it was explicity that the 18S region was targeted
# metadata.df[metadata.df$region_18S_targeted == "",]$region_18S_targeted <- "no"
metadata.df <- subset(metadata.df, region_18S_targeted == "no")

dim(metadata.df)
length(unique(metadata.df$study_accession))

# --------------------------------------------------------------------
# Make empty cells NA
metadata.df[metadata.df == ''] <- NA


# Make the index the rowname
rownames(metadata.df) <- metadata.df$Index



# Remove samples that are not in the list of projects
# mydir <- "data/feature_stats_AP19_separate_downsampled"
# myfiles <- list.files(mydir,pattern = "*.csv")
# projects <- gsub("_features.*csv","",myfiles)
# metadata.df[metadata.df$study_accession %in% projects,]

# --------------------------------------------------------------------------------------------------------------------------------

feature_tables_dir <- "data/feature_stats_AP19_separate_downsampled"
my_feature_table_files <- list.files(feature_tables_dir)
# my_feature_table_files<- my_feature_table_files[grepl("PRJNA255485",my_feature_table_files)]
# my_feature_table_files<- my_feature_table_files[grepl("PRJNA450848",my_feature_table_files)]
# my_feature_table_files <- my_feature_table_files[grepl("PRJEB30328",my_feature_table_files)]
# my_feature_table_files <- my_feature_table_files[grepl("PRJNA470773",my_feature_table_files)]
# raw_project_otu_table.df <- data.frame()
cnt = 0
# Load and process the OTU table
for (feature_table_file in my_feature_table_files){
  
  # Load and process the OTU table
  full_path <- paste0(feature_tables_dir, "/",feature_table_file)
  project_otu_table.df <- read.csv(full_path)
  # raw_project_otu_table.df <- project_otu_table.df
  project_name <- gsub("_.*","", feature_table_file)
  
  
  if (! project_name %in% metadata.df$study_accession){
    print(paste0("Not processing project : ", project_name))
    next
  }
  # else{
    print(paste0("Processing project : ", project_name))
    cnt = cnt + 1
    # next
  # }
  # Create output dirs for the project
  create_project_result_dirs(project_name)
  
  # Get the just metadata for the project
  project_metadata.df <- subset(metadata.df, study_accession == project_name)

    # Fix name of first column
  names(project_otu_table.df)[1] <- "OTU.ID"
  
  
  # Get the sample ids from the OTU table
  sample_ids_original <- names(project_otu_table.df)[!names(project_otu_table.df) %in% c("OTU.ID",
                                                            "Frequency",
                                                            "Domain","Phylum","Class","Order","Family","Genus","Species",
                                                            "Confidence","RepSeq",
                                                            "Taxon",
                                                            "taxonomy_species",
                                                            "taxonomy_genus",
                                                            "taxonomy_family",
                                                            "taxonomy_order",
                                                            "taxonomy_class",
                                                            "taxonomy_phylum")]
  sample_ids <- sample_ids_original
  # sample_ids_original <- grep("R[0-9].*|S[AB][0-9].*|S[0-9].*", names(project_otu_table.df), value = T)
  # sample_ids <- grep("R[0-9].*|S[AB][0-9].*|S[0-9].*", names(project_otu_table.df), value = T)
  
  # print(paste0("There are ", length(sample_ids_original), " samples in the data"))
  
  # ------------------------------------------------------------------------------------------
  # Results from the ACE amplicon pipeline `should' contain at least one observation/count in every row, however just to be sure
  # remove any rows containing all zeros. To do this, simply keep any row where there is any value not equal to zero.
  # project_otu_table.df[sample_ids] will return all columns with names matching the sample ids
  # The command below will take each row (MARGIN = 1) for the sample columns and check if any value is not zero.
  project_otu_table.df <- project_otu_table.df[apply(project_otu_table.df[sample_ids], MARGIN = 1, function(z) any(z!=0)),]
  
  # Split the Taxon column into Domain, Phylum...Species
  project_otu_table.df <- separate(project_otu_table.df, "Taxon", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"), remove =F, sep = ";")
  
  # Splitting taxa strings that are not specified at certain taxa levels will produce NA entries at those levels. 
  # NA entries should be changed to "Unassigned"
  project_otu_table.df[is.na(project_otu_table.df)] <- "Unassigned"
  
  # Replace D_# at beginning of taxon rank string to corresponding taxa label, e.g. D_0 = d (domain) D_1 = p (phylum), etc.
  # This doesn't work (well) for non-bacterial entries, which have more than the standard 7 levels
  project_otu_table.df$Domain <- as.character(lapply(project_otu_table.df$Domain, FUN = function(x) gsub("D_0", "d", x)))
  project_otu_table.df$Phylum <- as.character(lapply(project_otu_table.df$Phylum, FUN = function(x) gsub("D_1", "p", x)))
  project_otu_table.df$Class <- as.character(lapply(project_otu_table.df$Class, FUN = function(x) gsub("D_2", "c", x)))
  project_otu_table.df$Order <- as.character(lapply(project_otu_table.df$Order, FUN = function(x) gsub("D_3", "o", x)))
  project_otu_table.df$Family <- as.character(lapply(project_otu_table.df$Family, FUN = function(x) gsub("D_4", "f", x)))
  project_otu_table.df$Genus <- as.character(lapply(project_otu_table.df$Genus, FUN = function(x) gsub("D_5", "g", x)))
  project_otu_table.df$Species <- as.character(lapply(project_otu_table.df$Species, FUN = function(x) gsub("D_6", "s", x)))
  
  # Recreate the full taxonomy string with the 'prettier' taxa labels
  project_otu_table.df$taxonomy_species <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep =";"))
  
  # Also create a taxonomy string up to the genus level, since species are very rarely characterised at the Specie level in amplicon data
  project_otu_table.df$taxonomy_genus <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, Family, Genus, sep =";"))
  
  # And just for easier plotting later, create taxonomy strings for phylum, class, order and family levels
  project_otu_table.df$taxonomy_family <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, Family, sep =";"))
  project_otu_table.df$taxonomy_order <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, sep =";"))
  project_otu_table.df$taxonomy_class <- with(project_otu_table.df, paste(Domain, Phylum, Class, sep =";"))
  project_otu_table.df$taxonomy_phylum <- with(project_otu_table.df, paste(Domain, Phylum, sep =";"))
  
  # Store a version of the unfiltered project table
  project_otu_table_unfiltered.df <- project_otu_table.df
  
  # ------------------------------------------------------------------------------------------
  ## Now ensure the metadata and the OTU table matches
  
  # Sanity check whether all the samples are in the metadata, this should return nothing. 
  # Otherwise fix the sample names in the metadata or simply remove the offending samples.
  # Possibly have filtered out these samples at this point based on their sampletype
  missing_samples.v <- sample_ids[!sample_ids %in% project_metadata.df$Index]
  # if ( length(missing_samples.v) == 0) {
  #   print("No samples missing")
  # } else {
  #   print("Samples missing from OTU table")
  # }
  
  # Also the reverse, are there samples in the metadata missing from the data (possibly filtered earlier)
  missing_samples_metadata.v <- as.character(project_metadata.df$Index[!project_metadata.df$Index %in% sample_ids])

  # if ( length(missing_samples_metadata.v) == 0) {
  #   print("No samples missing")
  # } else {
  #   print("Samples missing from metadata")
  # }
  # Remove samples from the project table that are not in the metadata
  # print(dim(project_otu_table.df))
  project_otu_table.df <- project_otu_table.df[, !colnames(project_otu_table.df) %in% missing_samples.v]
  # print(dim(project_otu_table.df))
  
  # Reassign the sample ids 
  sample_ids <- names(project_otu_table.df)[!names(project_otu_table.df) %in% c("OTU.ID",
                                                                                "Frequency",
                                                                                "Domain","Phylum","Class","Order","Family","Genus","Species",
                                                                                "Confidence","RepSeq",
                                                                                "Taxon",
                                                                                "taxonomy_species",
                                                                                "taxonomy_genus",
                                                                                "taxonomy_family",
                                                                                "taxonomy_order",
                                                                                "taxonomy_class",
                                                                                "taxonomy_phylum")]
  # sample_ids <- grep("R[0-9].*|S[AB][0-9].*|S[0-9].*", names(project_otu_table.df), value = T)
  
  
  # ------------------------------------------------------------------------------------------
  # ----------- Remove unwanted lineages -----------
  
  # Remove OTUs that are Unassigned
  # project_otu_table.df <- project_otu_table.df[project_otu_table.df$Taxon != "Unassigned",]
  
  # Discard anything not Bacterial or Fungal
  project_otu_table.df <- project_otu_table.df[grepl("D_0__Bacteria|D_3__Fungi", project_otu_table.df$Taxon),]
  
  # Discard Chloroplast entries
  # dim(project_otu_table.df)
  # D_3__Chloroplast distribution? Arachis?
  # project_otu_table.df <- project_otu_table.df[!grepl("D_3__Chloroplast", project_otu_table.df$Taxon),]
  
  # dim(project_otu_table.df)
  # unique(grep("D_3__Chloroplast", project_otu_table.df$Taxon, value = T))
  # "c__Oxyphotobacteria;o__Chloroplast"
  
  # Discard anything not Bacterial or Fungal or Unassigned
  # project_otu_table.df <- project_otu_table.df[grepl("D_0__Bacteria|D_3__Fungi|^Unassigned$", project_otu_table.df$Taxon),]
  
  # Discard anything not Bacterial
  # project_otu_table.df <- project_otu_table.df[grepl("D_0__Bacteria", project_otu_table.df$Taxon),]
  # ------------------------------------------------
  
  # Remove old Taxon column
  project_otu_table.df$Taxon <- NULL
  project_otu_table_unfiltered.df$Taxon <- NULL
  
  # Store the OTUs and corresponding taxonomy information in a separate dataframe
  otu_taxonomy_map.df <- project_otu_table.df[c("OTU.ID",
                                                "Domain", 
                                                "Phylum", 
                                                "Class", 
                                                "Order", 
                                                "Family",
                                                "Genus",
                                                "Species",
                                                "taxonomy_species", 
                                                "taxonomy_genus",
                                                "taxonomy_family",
                                                "taxonomy_order",
                                                "taxonomy_class",
                                                "taxonomy_phylum",
                                                "RepSeq")]
  
  # Can save this table for later use if required
  write.table(otu_taxonomy_map.df, file = paste0("Result_tables/",project_name,"/other/", project_name,"_otu_taxonomy_map.csv"), sep = ",", quote = F, row.names = F)
  
  # Also save the unfiltered table, to avoid processing the original data table again 
  write.table(project_otu_table_unfiltered.df, file = paste0("Result_tables/", project_name,"/other/", project_name,"_otu_table_unfiltered.csv"), sep = ",", quote = F, row.names = F)
  
  
  # ---------------------------------------------------------------------------------------------------------------
  # ---------------------------------------------------------------------------------------------------------------
  # Now we can generate the tables that we will need for different analyses at both the OTU and various taxa levels
  
  # -----------------------------
  # -----------OTU LEVEL---------
  # Dataframe containing the counts for each OTU, for each sample
  otu.df <- project_otu_table.df[c("OTU.ID", sample_ids)]
  otu_unfiltered.df <- project_otu_table_unfiltered.df[c("OTU.ID", sample_ids)]
  
  ## Create a matrix version for ease of processing
  # For filtered OTU matrix
  otu.m <- otu.df
  rownames(otu.m) <- otu.m$OTU.ID
  otu.m$OTU.ID <- NULL
  otu.m <- as.matrix(otu.m)
  
  # And unfiltered OTU matrix
  otu_unfiltered.m <- otu_unfiltered.df
  rownames(otu_unfiltered.m) <- otu_unfiltered.m$OTU.ID
  otu_unfiltered.m$OTU.ID <- NULL
  otu_unfiltered.m <- as.matrix(otu_unfiltered.m)
  
  # Create relative abundance matrix from counts matrix
  otu_rel.m <- t(t(otu.m)/ colSums(otu.m))
  otu_unfiltered_rel.m <- t(t(otu_unfiltered.m) / colSums(otu_unfiltered.m))
  
  # Change nans to 0. Occurs when a sample has no hits at this point.
  otu_rel.m[is.nan(otu_rel.m)] <- 0
  otu_unfiltered_rel.m[is.nan(otu_unfiltered_rel.m)] <- 0

  # print(colSums(otu_unfiltered.m[,"ERR3182717_J001", drop =F]))
  # print(colSums(otu_unfiltered_rel.m[,"ERR3182717_J001", drop =F]))
  # print(colSums(otu.m[,"ERR3182717_J001", drop =F]))
  # print(colSums(otu_rel.m[,"ERR3182717_J001", drop =F]))
  
  # --------------------------------------------------------------------------------
  # Filter those OTUs that are low abundance in all samples
  # Check how many OTUs would be removed if we filtered any whose abundance is less than 0.05% (0.0005) in all samples
  filter_percentage = 0.0005
  otu_rel_low_abundance_otus.m <- otu_rel.m[apply(otu_rel.m[,sample_ids,drop =F],1,function(z) all(z < filter_percentage)),]
  
  
  # TODO - Collect low abundance OTUs that will be filtered out at this stage
  # FIXME write.table(as.data.frame(otu_rel_low_abundance_otus.m), file = "Result_tables/count_tables/low_abundance_OTUs.csv", sep = ",", quote = F, col.names = T, row.names = F)
  percent_removed_before_taxa_filtered <- round(dim(otu_rel_low_abundance_otus.m)[1]/length(rownames(otu_unfiltered.m)) * 100,2)
  percent_removed_after_taxa_filtered <- round(dim(otu_rel_low_abundance_otus.m)[1]/length(rownames(otu.m)) * 100,2)
  
  print(paste("There are a total of", dim(otu_rel.m)[1], "OTUs before filtering"))
  print(paste("A total of", dim(otu_rel_low_abundance_otus.m)[1], 
              "OTUs will be filtered at a threshold of", 
              filter_percentage * 100, "percent"))
  
  # print(colSums(otu_rel.m[,"ERR3182717_J001", drop =F]))
  # print(colSums(otu.m[,"ERR3182717_J001", drop =F]))
  # If you are happy with the filtering threshold, apply it.
  otu_rel.m <- otu_rel.m[apply(otu_rel.m[,sample_ids, drop =F],1,function(z) any(z>=filter_percentage)),,drop = F]
  
  # Re-normalise the matrix after filtering
  otu_rel.m <- t(t(otu_rel.m) / colSums(otu_rel.m))
  
  # Change nans to 0. Occurs when a sample has no hits at this point.
  otu_rel.m[is.nan(otu_rel.m)] <- 0
  
  # Also remove low abundance OTUs from the original OTU count matrix
  otu.m  <- otu.m[rownames(otu_rel.m),,drop = F]
  
  # Remove those samples with less than # reads
  # print(colSums(otu_prior_to_removing_low_read_count_samples.m[,"ERR3182717_J001", drop =F]))
  otu_prior_to_removing_low_read_count_samples.m <- otu.m
  otu_prior_to_removing_low_read_count_samples_rel.m <- otu_rel.m
  otu.m <- otu.m[,colSums(otu.m) >= 5000, drop = F]
  # The might be many rows whos maximum is 0 at this point. Remove them.
  otu.m <- otu.m[apply(otu.m, 1, max) != 0,,drop = F]
  
  
  # -------------------------------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------------------------------
  # Get the most abundant unassigned features
  # Project table with just unassigned
  unassigned_project_otu_table_unfiltered.df <- project_otu_table_unfiltered.df[project_otu_table_unfiltered.df$Domain == "Unassigned",]
  # Convert to dataframe
  unassigned_otu_unfiltered_rel.df <- m2df(otu_unfiltered_rel.m[as.character(unassigned_project_otu_table_unfiltered.df$OTU.ID),,drop =F],
                                           name_of_taxonomy_col = "OTU.ID")
  # Melt
  unassigned_otu_unfiltered_rel.df <- melt(unassigned_otu_unfiltered_rel.df, variable.name = "Sample", value.name = "Relative_abundance")
  
  # Get the top most abundant unassigned features per sample
  most_abundant_unassigned.df <- unassigned_otu_unfiltered_rel.df %>% 
    group_by(Sample) %>%
    filter(Relative_abundance > 0) %>%
    top_n(n = 10, wt = Relative_abundance) %>% 
    mutate(Relative_abundance = round(Relative_abundance*100,3)) %>%
    as.data.frame() 
  
  # Get corresponding representative sequence
  most_abundant_unassigned.df$RepSeq <- unlist(lapply(most_abundant_unassigned.df$OTU.ID, function(x) as.character(unassigned_project_otu_table_unfiltered.df[unassigned_project_otu_table_unfiltered.df$OTU.ID == x,]$RepSeq)))
  write.csv(x = most_abundant_unassigned.df, file = paste0("Result_tables/",project_name,"/other/",project_name,"_most_abundant_unassigned.csv"), row.names = F)
  
  # Get unique set of features
  unique_most_abundant_unassigned.df <- unique(most_abundant_unassigned.df[c("OTU.ID", "RepSeq")])
  
  # Write fasta file
  write.fasta(sequences = as.list(unique_most_abundant_unassigned.df$RepSeq),open = "w", 
              names = as.character(unique_most_abundant_unassigned.df$OTU.ID),
              file.out = paste0("Result_tables/",project_name,"/other/",project_name,"_most_abundant_unassigned_features.fasta"))
  
  
  # -------------------------------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------------------------------
  #         Rarefying
  
  # Note - For samples with a read count lower then the sample=# parameter, 
  # the rrarefy function in vegan will return the this sample with its existing read count.
  # Sample with counts higher than the sample parameter will be rarefied as normal and their counts will be capped at the sample parameter value
  
  # Normally samples with OTU counts lower than our desired threshold need to be removed if we want to rarefy.
  # Here we are simply using rrarefy to cap the read depth so the fold difference between the highest and lowest is not extreme.
  
  # Looking at rarefaction curve and the counts for each sample can be a good
  # way to determine a good minimum/maximum count threshold for samples
  
  # Counts for each sample
  column_sums <- colSums(otu.m)
  column_sums.df <- melt(column_sums[order(column_sums)])
  column_sums.df$sample <- rownames(column_sums.df)
  rownames(column_sums.df) <- NULL
  column_sums.df <- column_sums.df[c("sample", "value")]
  column_sums.df$sample <- factor(column_sums.df$sample, levels = column_sums.df$sample)
  
  myplot <- ggplot(column_sums.df, aes(x = sample, y = value)) + 
    geom_histogram(stat = "identity") +
    geom_hline(yintercept = 30000, color = 'red')+
    # geom_hline(yintercept = 20000, color = 'red')+
    geom_hline(yintercept = mean(column_sums.df$value), color = 'blue')+
    geom_hline(yintercept = median(column_sums.df$value), color = 'purple')+
    scale_y_continuous(breaks = seq(0,max(column_sums.df$value), 4000)) +
    xlab("Sample") +
    ylab("Read count") +
    common_theme +
    theme(axis.text.x = element_text(angle = 90,vjust = .5, size = 3))
  ggsave(plot = myplot, filename = paste0("./Result_figures/",project_name,"/exploratory_analysis/sample_read_depth_distribution.pdf"), width=100, height=15, units = "cm")
  
  myplot <- ggplot(column_sums.df, aes(x = value)) + 
    xlab("Read count") +
    ylab("Number of samples") +
    scale_x_continuous(breaks = seq(4000,max(column_sums.df$value), 4000)) +
    scale_y_continuous(breaks = seq(0,length(column_sums.df$value), 5)) +
    geom_histogram(stat = "bin", bins = 60, colour = "black",fill = "grey") +
    geom_vline(xintercept = 30000, color = 'red')+
    common_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
  ggsave(plot = myplot, filename = paste0("./Result_figures/",project_name,"/exploratory_analysis/sample_read_depth_distribution_2.pdf"), width=20, height=20, units = "cm")

  otu_rare_count.m <- t(rrarefy(x = t(otu.m), sample=50000))
  
  # -------------------------------------------------------------------------------------------------------------------------------------
  
  
  # And re-calculate the abundances after filtering
  otu_rel.m <- t(t(otu.m)/ colSums(otu.m))
  otu_rel.m[is.nan(otu_rel.m)] <- 0
  otu_rel_rare.m <- t(t(otu_rare_count.m) /colSums(otu_rare_count.m))
  otu_rel_rare.m[is.nan(otu_rel_rare.m)] <- 0
  
  # Reassign sample IDs
  sample_ids <- colnames(otu_rel.m)
  
  ## (Re)build dataframes
  
  # Un-rarefied
  otu_rel.df <- data.frame("OTU.ID" = rownames(otu_rel.m))
  otu_rel.df <- cbind(otu_rel.df, otu_rel.m[,colnames(otu_rel.m),drop =F])
  rownames(otu_rel.df) <- c()
  
  otu.df <- data.frame("OTU.ID" = rownames(otu.m))
  otu.df <- cbind(otu.df, otu.m[,colnames(otu.m),drop =F])
  rownames(otu.df) <- c()
  
  # Rarified
  otu_rel_rare.df <- data.frame("OTU.ID" = rownames(otu_rel_rare.m))
  otu_rel_rare.df <- cbind(otu_rel_rare.df, otu_rel_rare.m[,colnames(otu_rel_rare.m),drop =F])
  rownames(otu_rel_rare.df) <- c()
  
  otu_rare.df <- data.frame("OTU.ID" = rownames(otu_rare_count.m))
  otu_rare.df <- cbind(otu_rare.df, otu_rare_count.m[,colnames(otu_rare_count.m),drop =F])
  rownames(otu_rare.df) <- c()
  
  # Write the final otu counts and abundances to file
  write.table(otu.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_OTU_counts.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_rel.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_OTU_relative_abundances.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_rare.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_OTU_counts_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_rel_rare.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_OTU_relative_abundances_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  # NOTE - otu.df, otu.m, otu_rel.df and otu_rel.m are the final, filtered OTU count / relative abundance dataframes and matrices. These can be used
  # elsewhere for a variety of analyses at the OTU level, or, as is shown below, used to calculate abundances at different taxa levels
  
  # Sample removed
  samples_retained <- colnames(otu.df)[2:length(colnames(otu.df))]
  samples_lost <- project_metadata.df$Index[!project_metadata.df$Index %in% samples_retained]
  print(paste(length(samples_retained), "samples retained"))
  print(paste(length(samples_lost), "samples lost"))
  
  # Label those samples from the metadata.df that are not in the OTU table. A sample that has been filtered out
  # will generally not be of interest for the study, though we may wish to have the metadata at hand
  project_metadata.df$Sample_retained <- "no"
  project_metadata.df[project_metadata.df$Index %in% samples_retained,]$Sample_retained <- "yes"
  write.table(project_metadata.df[project_metadata.df$Sample_retained == "yes",], file = paste0("Result_tables/",project_name,"/other/",project_name,"_processed_metadata.csv"), sep = ",", quote = F, row.names = F)
  write.table(project_metadata.df[project_metadata.df$Sample_retained == "no",], file = paste0("Result_tables/",project_name,"/other/",project_name,"_metadata_samples_removed.csv"), sep = ",", quote = F, row.names = F)
  
  # ------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------
  #           Reads stats
  
  # Rarefied counts                                  = otu_rare_count.m 
  # Rarefied relative abundances                     = otu_rare_rel.m
  # Relative abundances                              = otu_rel.m
  # Counts (filtered)                                = otu.m
  # Unfiltered counts                                = otu_unfiltered.m
  # Unfiltered relative abundances                   = otu_unfiltered_rel.m
  # Counts (with low abundance samples)              = otu_prior_to_removing_low_read_count_samples.m
  # Relative abundances (with low abundance samples) = otu_prior_to_removing_low_read_count_samples_rel.m
  
  
  # Only focus on those samples that have passed QC.
  samples_passing_QC <- colnames(otu.m)
  samples_in_unfiltered <- colnames(otu_unfiltered.m)
  # stats.df <- data.frame(Sample = samples_passing_QC)
  stats.df <- data.frame(Sample = samples_in_unfiltered)
  rownames(stats.df) <- stats.df$Sample
  
  stats.df$Sample_retained <- "no"
  stats.df[samples_passing_QC,]$Sample_retained <- "yes"
  samples_not_retained <- rownames(stats.df[stats.df$Sample_retained == "no",])
  
  stats.df$study_accession <- project_metadata.df[samples_in_unfiltered,"study_accession"]
  stats.df$run_accession <- project_metadata.df[samples_in_unfiltered,"run_accession"]
  stats.df$scientific_name <- project_metadata.df[samples_in_unfiltered,"scientific_name"]
  stats.df$Commodity <- project_metadata.df[samples_in_unfiltered,"Commodity"]
  stats.df$Sample_type <- project_metadata.df[samples_in_unfiltered,"Sample_type"]
  stats.df$Sample_treatment <- project_metadata.df[samples_in_unfiltered,"Sample_treatment"]
  stats.df$Top_region_from_BLAST <- project_metadata.df[samples_in_unfiltered,"Top_region_from_BLAST"]
  
  stats.df[,"Original_read_counts"] <- colSums(otu_unfiltered.m[,samples_in_unfiltered, drop = F]) # Read counts prior to filtering out taxonomy / contaminants / low abundance features
  stats.df[samples_passing_QC,"Filtered_read_counts"] <- colSums(otu.m[,samples_passing_QC, drop = F])  # Read counts after filtering (contaminants, low read depth samples and low abundance features removed)
  stats.df[samples_not_retained,"Filtered_read_counts"] <- colSums(otu_prior_to_removing_low_read_count_samples.m[,rownames(stats.df[samples_not_retained,, drop = F]), drop = F])
  stats.df[samples_passing_QC,"Filtered_rarefied_read_counts"] <- colSums(otu_rare_count.m[,samples_passing_QC, drop = F]) # "" and rarefied
  
  # Reads removed from filtering
  stats.df[,"Reads_removed_from_filtering"] <- stats.df[,"Original_read_counts"] - stats.df[,"Filtered_read_counts"]
  stats.df[,"Proportion_reads_removed_from_filtering"] <- stats.df[,"Reads_removed_from_filtering"] / stats.df[,"Original_read_counts"]
  
  # Reads removed from filtering and rarefaction
  stats.df[,"Reads_removed_from_filtering_rarefied"] <- stats.df[,"Original_read_counts"] - stats.df[,"Filtered_rarefied_read_counts"]
  stats.df[,"Proportion_reads_removed_from_filtering_rarefied"] <- stats.df[,"Reads_removed_from_filtering_rarefied"] / stats.df[,"Original_read_counts"]
  
  # ---------------------------------------------
  # Read counts and proportions original, unprocessed data. Summing each domain + unassigned should equal one
  
  # Proportion of reads that are mammal (generally human)
  stats.df[,"Mammal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("Mammal", project_otu_table_unfiltered.df$taxonomy_species),samples_in_unfiltered, drop = F])
  stats.df[,"Mammal_proportion_original"] <- stats.df[,"Mammal_read_count_original"] / stats.df[,"Original_read_counts"]
  
  # Proportion of reads that are fungal
  stats.df[,"Fungal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("Fungi", project_otu_table_unfiltered.df$taxonomy_species),samples_in_unfiltered, drop = F])
  stats.df[,"Fungal_proportion_original"] <- stats.df[,"Fungal_read_count_original"] / stats.df[,"Original_read_counts"]
  
  # Proportion of reads that are fungal
  stats.df[,"Bacterial_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("d__Bacteria", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered, drop = F])
  stats.df[,"Bacterial_proportion_original"] <- stats.df[,"Bacterial_read_count_original"] / stats.df[,"Original_read_counts"]
  
  # Proportion of reads that are Archaea
  stats.df[,"Archaeal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("d__Archaea", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered, drop = F])
  stats.df[,"Archaeal_proportion_original"] <- stats.df[,"Archaeal_read_count_original"] / stats.df[,"Original_read_counts"]
  
  # Proportion of reads that are Eukaryota
  stats.df[,"Eukaryal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("d__Eukaryota", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered, drop = F])
  stats.df[,"Eukaryal_proportion_original"] <- stats.df[,"Eukaryal_read_count_original"] / stats.df[,"Original_read_counts"]
  
  # Proportion of reads that are Unassigned a taxonomy
  stats.df[,"Unassigned_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("Unassigned", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered, drop = F])
  stats.df[,"Unassigned_proportion_original"] <- stats.df[,"Unassigned_read_count_original"] / stats.df[,"Original_read_counts"]
  
  # -------------------------------------
  # Read counts and proportions filtered data
  # temp <- otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID %in% rownames(otu.m),]
  # otu_prior_to_removing_low_read_count_samples.m should be the same as otu.m except the low abundance samples have not been removed
  temp <- otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID %in% rownames(otu_prior_to_removing_low_read_count_samples.m),]
  bacterial_otus_filtered <- temp[temp$Domain == "d__Bacteria",]$OTU.ID
  fungal_otus_filtered <- temp[grepl("Fungi", temp$taxonomy_species),]$OTU.ID
  # names(colSums(otu_prior_to_removing_low_read_count_samples.m[!rownames(otu_prior_to_removing_low_read_count_samples.m) %in% rownames(otu.m),])[colSums(otu_prior_to_removing_low_read_count_samples.m[!rownames(otu_prior_to_removing_low_read_count_samples.m) %in% rownames(otu.m),]) > 0]) %in% colnames(otu.m)
  
  # Proportion of reads in the filtered data that are bacterial
  # temp <- stats.df
  stats.df[,"Bacterial_read_count_after_filtering"] <- colSums(otu_prior_to_removing_low_read_count_samples.m[which(rownames(otu_prior_to_removing_low_read_count_samples.m) %in% bacterial_otus_filtered),, drop = F])
  # summary(stats.df[samples_passing_QC,"Bacterial_read_count_after_filtering"] == temp[samples_passing_QC,"Bacterial_read_count_after_filtering"])
  stats.df[,"Bacterial_proportion_after_filtering"] <- stats.df[,"Bacterial_read_count_after_filtering"] / stats.df[,"Filtered_read_counts"]
  
  # Proportion of reads in the filtered data that are Fungal
  stats.df[,"Fungal_read_count_after_filtering"] <- colSums(otu_prior_to_removing_low_read_count_samples.m[which(rownames(otu_prior_to_removing_low_read_count_samples.m) %in% fungal_otus_filtered),, drop = F])
  stats.df[,"Fungal_proportion_after_filtering"] <- stats.df[,"Fungal_read_count_after_filtering"] / stats.df[,"Filtered_read_counts"]
  
  # -------------------------------------
  # Read counts and proportions filtered + rarefied data
  temp <- otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID %in% rownames(otu_rare_count.m),]
  bacterial_otus_filtered_rarefied <- temp[temp$Domain == "d__Bacteria",]$OTU.ID
  fungal_otus_filtered_rarefied <- temp[grepl("Fungi", temp$taxonomy_species),]$OTU.ID
  
  # Proportion of reads in the filtered data that are bacterial
  stats.df[samples_passing_QC,"Bacterial_read_count_after_filtering_rarefied"] <- colSums(otu_rare_count.m[which(rownames(otu_rare_count.m) %in% bacterial_otus_filtered_rarefied),samples_passing_QC, drop = F])
  stats.df[samples_passing_QC,"Bacterial_proportion_after_filtering_rarefied"] <- stats.df[samples_passing_QC,"Bacterial_read_count_after_filtering_rarefied"] / stats.df[samples_passing_QC,"Filtered_rarefied_read_counts"]
  
  # Proportion of reads in the filtered data that are Fungal
  stats.df[samples_passing_QC,"Fungal_read_count_after_filtering_rarefied"] <- colSums(otu_rare_count.m[which(rownames(otu_rare_count.m) %in% fungal_otus_filtered_rarefied),samples_passing_QC, drop = F])
  stats.df[samples_passing_QC,"Fungal_proportion_after_filtering_rarefied"] <- stats.df[samples_passing_QC,"Fungal_read_count_after_filtering_rarefied"] / stats.df[samples_passing_QC,"Filtered_rarefied_read_counts"]
  
  ## ------------------------------------------------------------------------------------
  ## This will calculate the total number of features across all samples and the number of non-zero features in each sample
  # Can either calculate the feature numbers on just samples passing QC
  
  # Or calculate on all samples
  stats.df[,"Features_total"] <- length(which(rowSums(otu_unfiltered.m) > 0 ))
  stats.df[,"Features_original"] <- apply(otu_unfiltered.m, 2, function(x) { length(which(x > 0)) } )
  ## ------------------------------------------------------------------------------------
  stats.df[samples_passing_QC,"Features_filtered"] <- apply(otu.m[,samples_passing_QC,drop=F], 2, function(x) { length(which(x > 0)) } )
  stats.df[samples_not_retained,"Features_filtered"] <- apply(otu_prior_to_removing_low_read_count_samples.m[,samples_not_retained,drop=F], 2, function(x) { length(which(x > 0)) } )
  stats.df[samples_passing_QC,"Features_filtered_rarefied"] <- apply(otu_rare_count.m[,samples_passing_QC,drop=F], 2, function(x) { length(which(x > 0)) } )
  
  stats.df[,"Features_removed_from_filtering"] <- stats.df[,"Features_original"] - stats.df[,"Features_filtered"] 
  stats.df[,"Features_removed_from_filtering_rarefied"] <- stats.df[,"Features_original"] - stats.df[,"Features_filtered_rarefied"] 
  
  stats.df[,"Proportion_features_removed_from_filtering"] <- stats.df[,"Features_removed_from_filtering"] / stats.df[,"Features_original"]
  stats.df[,"Proportion_features_removed_from_filtering_rarefied"] <- stats.df[,"Features_removed_from_filtering_rarefied"] / stats.df[,"Features_original"]
  
  write.csv(stats.df, paste0("Result_tables/",project_name,"/other/",project_name,"_QC_summary.csv"), row.names = F, quote = F)
  
  
  # ------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------

  # Above we processed the frequencies for each OTU to calculate the relative abundances.
  # However, we often want the abundances at not just the individual OTU level, but also different taxonomy levels.
  # For example, we may want to know the abundance of a particular Family.
  # Now we will generate the abundance tables at each taxonomy level from Phylum, Class, Order, Family and Genus
  
  # Merge the otu counts back with the taxonomy data
  otu_metadata_merged.df <- merge(otu.df, otu_taxonomy_map.df, by.x = "OTU.ID", by.y = "OTU.ID")
  # temp <- left_join(otu.df, otu_taxonomy_map.df, by = "OTU.ID")
  
  # Do the same for the rarefied data
  otu_metadata_merged_rare.df <- merge(otu_rare.df, otu_taxonomy_map.df, by.x = "OTU.ID", by.y = "OTU.ID")
  
  # We are then going to create the dataframe for each taxonomy level containing the relative abundances.
  # e.g. otu_phylum_rel.df,otu_class_rel.df etc.
  # and the counts
  # e.g. otu_phylum.df,otu_class.df etc.
  
  # We use the 'full' taxonomy strings, e.g. taxonomy_genus, so that "Unassigned" or "Uncultured" from different lineages are kept separate!
  for (tax_string_level in c("taxonomy_species", "taxonomy_genus", "taxonomy_family", "taxonomy_class", "taxonomy_order", "taxonomy_phylum")){
    # Collapse the dataframe by summing the counts for each unique taxonomy string within each sample
    otu_taxa_level.df <- otu_metadata_merged.df[c(tax_string_level, sample_ids)] %>% 
      group_by_(tax_string_level) %>% # Group the dataframe by the taxonomy string 
      # Summarise each group by applying the the 'sum' function to the counts of each member of the group, 
      # i.e. duplicate entries for each taxa level are collapsed into a single entry and their counts summed together
      dplyr::summarise_all(funs(sum)) %>%  
      as.data.frame() # convert back to dataframe
    
    otu_taxa_level_rare.df <- otu_metadata_merged_rare.df[c(tax_string_level, sample_ids)] %>%
      group_by_(tax_string_level) %>%
      dplyr::summarise_all(funs(sum)) %>%
      as.data.frame()
    
    # Now create the relative abundance matrix at the current taxa level
    otu_taxa_level_rel.m <- otu_taxa_level.df
    rownames(otu_taxa_level_rel.m) <- otu_taxa_level_rel.m[[tax_string_level]]
    otu_taxa_level_rel.m[tax_string_level] <- NULL
    otu_taxa_level_rel.m <- as.matrix(otu_taxa_level_rel.m)
    otu_taxa_level_rel.m <- t(t(otu_taxa_level_rel.m) / colSums(otu_taxa_level_rel.m))
    
    otu_taxa_level_rel_rare.m <- otu_taxa_level_rare.df
    rownames(otu_taxa_level_rel_rare.m) <- otu_taxa_level_rel_rare.m[[tax_string_level]]
    otu_taxa_level_rel_rare.m[tax_string_level] <- NULL
    otu_taxa_level_rel_rare.m <- as.matrix(otu_taxa_level_rel_rare.m)
    otu_taxa_level_rel_rare.m <- t(t(otu_taxa_level_rel_rare.m) / colSums(otu_taxa_level_rel_rare.m))
    
    if (grepl("phylum", tax_string_level)){
      # otu_phylum_rel.m <- otu_taxa_level_rel.m
      otu_phylum_rel.df <- m2df(otu_taxa_level_rel.m, tax_string_level)
      otu_phylum.df <- otu_taxa_level.df
      otu_phylum_rel_rare.df <- m2df(otu_taxa_level_rel_rare.m, tax_string_level)
      otu_phylum_rare.df <- otu_taxa_level_rare.df
    } 
    else if (grepl("class", tax_string_level)){
      # otu_class_rel.m <- otu_taxa_level_rel.m
      otu_class_rel.df <- m2df(otu_taxa_level_rel.m, tax_string_level)
      otu_class.df <- otu_taxa_level.df
      otu_class_rel_rare.df <- m2df(otu_taxa_level_rel_rare.m, tax_string_level)
      otu_class_rare.df <- otu_taxa_level_rare.df
    }
    else if (grepl("order", tax_string_level)){
      # otu_order_rel.m <- otu_taxa_level_rel.m
      otu_order_rel.df <- m2df(otu_taxa_level_rel.m, tax_string_level)
      otu_order.df <- otu_taxa_level.df
      otu_order_rel_rare.df <- m2df(otu_taxa_level_rel_rare.m, tax_string_level)
      otu_order_rare.df <- otu_taxa_level_rare.df
    }
    else if (grepl("family", tax_string_level)){
      # otu_family_rel.m <- otu_taxa_level_rel.m
      otu_family_rel.df <- m2df(otu_taxa_level_rel.m, tax_string_level)
      otu_family.df <- otu_taxa_level.df
      otu_family_rel_rare.df <- m2df(otu_taxa_level_rel_rare.m, tax_string_level)
      otu_family_rare.df <- otu_taxa_level_rare.df
    }
    else if (grepl("genus", tax_string_level)){
      # otu_genus_rel.m <- otu_taxa_level_rel.m
      otu_genus_rel.df <- m2df(otu_taxa_level_rel.m, tax_string_level)
      otu_genus.df <- otu_taxa_level.df
      otu_genus_rel_rare.df <- m2df(otu_taxa_level_rel_rare.m, tax_string_level)
      otu_genus_rare.df <- otu_taxa_level_rare.df
    }
    else if (grepl("species", tax_string_level)){
      # otu_species_rel.m <- otu_taxa_level_rel.m
      otu_species_rel.df <- m2df(otu_taxa_level_rel.m, tax_string_level)
      otu_species.df <- otu_taxa_level.df
      otu_species_rel_rare.df <- m2df(otu_taxa_level_rel_rare.m, tax_string_level)
      otu_species_rare.df <- otu_taxa_level_rare.df
    }
  }
  
  ### Write the final counts and abundances for each taxonomy level to file
  # Not rarefied
  write.table(otu_species.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Specie_counts.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_species_rel.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Specie_relative_abundances.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_genus.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Genus_counts.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_genus_rel.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Genus_relative_abundances.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_family.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Family_counts.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_family_rel.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Family_relative_abundances.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_order.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Order_counts.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_order_rel.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Order_relative_abundances.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_class.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Class_counts.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_class_rel.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Class_relative_abundances.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_phylum.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Phylum_counts.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_phylum_rel.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Phylum_relative_abundances.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  # Rarefied
  write.table(otu_species_rare.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Specie_counts_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_species_rel_rare.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Specie_relative_abundances_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_genus_rare.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Genus_counts_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_genus_rel_rare.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Genus_relative_abundances_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_family_rare.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Family_counts_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_family_rel_rare.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Family_relative_abundances_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_order_rare.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Order_counts_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_order_rel_rare.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Order_relative_abundances_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_class_rare.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Class_counts_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_class_rel_rare.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Class_relative_abundances_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  write.table(otu_phylum_rare.df, file = paste0("Result_tables/",project_name,"/count_tables/",project_name,"_Phylum_counts_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(otu_phylum_rel_rare.df, file = paste0("Result_tables/",project_name,"/relative_abundance_tables/",project_name,"_Phylum_relative_abundances_rarefied.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  
  # ------------------------------------------------------------------------------------------------------------------------------
  # Finally create and save a dataframe, separately for each Phylum, Class, Order, Family, Genus, Species and OTU,
  # containing the abundances/counts/log(counts, 10) in each sample, metadata and taxonomy information.
  
  # Include both the count, log(count), abundance, rarefied count, log(rarefied count) and rarefied abundance.
  
  # There should four tables at this point each level with the same entries, e.g.
  # Counts = otu.df
  # Relative abundance = otu_rel.df
  # Rarified counts = otu_rare.df
  # Rarified relative abundance = otu_rel_rare.df
  
  # df2matrix <- function(mydf, rowname_col = 1){
  #   temp <- mydf[,c(1:length(names(mydf)))[c(1:length(names(mydf))) != rowname_col]]
  #   mymatrix <- as.matrix(temp)
  #   rownames(mymatrix) <- mydf[,rowname_col]
  #   # mymatrix[,rowname_col] <- NULL
  #   return(mymatrix)
  # }
  
  reduced_tax_map <- otu_taxonomy_map.df
  reduced_tax_map$RepSeq <- NULL
  
  otu_combined <- create_combined_dataframe(counts.df = otu.df, 
                                            counts_rare.df = otu_rare.df,
                                            abundances.df = otu_rel.df,
                                            abundances_rare.df = otu_rel_rare.df,
                                            mylevel = "OTU.ID",
                                            mymetadata = project_metadata.df,
                                            otu_map.df = reduced_tax_map)
  
  species_combined <- create_combined_dataframe(counts.df = otu_species.df, 
                                                counts_rare.df = otu_species_rare.df,
                                                abundances.df = otu_species_rel.df,
                                                abundances_rare.df = otu_species_rel_rare.df,
                                                mylevel = "Species",
                                                mymetadata = project_metadata.df,
                                                otu_map.df = reduced_tax_map)
  
  genus_combined <- create_combined_dataframe(counts.df = otu_genus.df, 
                                              counts_rare.df = otu_genus_rare.df,
                                              abundances.df = otu_genus_rel.df,
                                              abundances_rare.df = otu_genus_rel_rare.df,
                                              mylevel = "Genus",
                                              mymetadata = project_metadata.df,
                                              otu_map.df = reduced_tax_map)
  
  family_combined <- create_combined_dataframe(counts.df = otu_family.df, 
                                               counts_rare.df = otu_family_rare.df,
                                               abundances.df = otu_family_rel.df,
                                               abundances_rare.df = otu_family_rel_rare.df,
                                               mylevel = "Family",
                                               mymetadata = project_metadata.df,
                                               otu_map.df = reduced_tax_map)
  
  order_combined <- create_combined_dataframe(counts.df = otu_order.df, 
                                              counts_rare.df = otu_order_rare.df,
                                              abundances.df = otu_order_rel.df,
                                              abundances_rare.df = otu_order_rel_rare.df,
                                              mylevel = "Order",
                                              mymetadata = project_metadata.df,
                                              otu_map.df = reduced_tax_map)
  
  class_combined <- create_combined_dataframe(counts.df = otu_class.df, 
                                              counts_rare.df = otu_class_rare.df,
                                              abundances.df = otu_class_rel.df,
                                              abundances_rare.df = otu_class_rel_rare.df,
                                              mylevel = "Class",
                                              mymetadata = project_metadata.df,
                                              otu_map.df = reduced_tax_map)
  
  phylum_combined <- create_combined_dataframe(counts.df = otu_phylum.df, 
                                               counts_rare.df = otu_phylum_rare.df,
                                               abundances.df = otu_phylum_rel.df,
                                               abundances_rare.df = otu_phylum_rel_rare.df,
                                               mylevel = "Phylum",
                                               mymetadata = project_metadata.df,
                                               otu_map.df = reduced_tax_map)
  
  write.table(otu_combined, file = paste0("Result_tables/",project_name,"/other/",project_name,"_OTU_counts_abundances_and_metadata.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(species_combined, file = paste0("Result_tables/",project_name,"/other/",project_name,"_Specie_counts_abundances_and_metadata.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(genus_combined, file = paste0("Result_tables/",project_name,"/other/",project_name,"_Genus_counts_abundances_and_metadata.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(family_combined, file = paste0("Result_tables/",project_name,"/other/",project_name,"_Family_counts_abundances_and_metadata.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(order_combined, file = paste0("Result_tables/",project_name,"/other/",project_name,"_Order_counts_abundances_and_metadata.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(class_combined, file = paste0("Result_tables/",project_name,"/other/",project_name,"_Class_counts_abundances_and_metadata.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  write.table(phylum_combined, file = paste0("Result_tables/",project_name,"/other/",project_name,"_Phylum_counts_abundances_and_metadata.csv"), sep = ",", quote = F, col.names = T, row.names = F)
}


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

# Now each individual project has been processed, generate combined files

# First dataframes that can be simply combined without melting
combine_dataframes <- function(base_location, mypattern){
  filenames <- list.files(base_location,pattern = mypattern,recursive = T,full.names = T)
  All <- lapply(filenames,function(i){read.csv(i, header=T)})
  df <- do.call(rbind.data.frame, All)
}

write.csv(x = combine_dataframes("Result_tables","P.*_processed_metadata.csv"), 
          file = "Result_tables/combined/combined_processed_metadata.csv", row.names = F, quote = F)

temp <- read.csv("Result_tables/combined/combined_processed_metadata.csv")
length(unique(temp$study_accession))
length(unique(temp$run_accession))

write.csv(x = combine_dataframes("Result_tables","P.*_OTU_counts_abundances_and_metadata.csv"),
          file = "Result_tables/combined/combined_OTU_counts_abundances_and_metadata.csv", row.names = F, quote = F)

write.csv(x = combine_dataframes("Result_tables","P.*_Genus_counts_abundances_and_metadata.csv"), 
          file = "Result_tables/combined/combined_Genus_counts_abundances_and_metadata.csv", row.names = F, quote = F)

write.csv(x = combine_dataframes("Result_tables","P.*_Family_counts_abundances_and_metadata.csv"), 
          file = "Result_tables/combined/combined_Family_counts_abundances_and_metadata.csv", row.names = F, quote = F)

write.csv(x = combine_dataframes("Result_tables","P.*_QC_summary.csv"), 
          file = "Result_tables/combined/combined_QC_summary.csv", row.names = F, quote = F)

write.csv(x = unique(combine_dataframes("Result_tables","P.*_otu_taxonomy_map.csv")), 
          file = "Result_tables/combined/combined_otu_taxonomy_map.csv", row.names = F, quote = F)


# write.csv(x = combine_dataframes("Result_tables","P.*_most_abundant_unassigned.csv"), 
          # file = "Result_tables/combined_most_abundant_unassigned.csv", row.names = F, quote = F)
temp <- combine_dataframes("Result_tables","P.*_most_abundant_unassigned.csv")
temp$study_accession <- unlist(lapply(as.character(temp$Sample), function(x) as.character(metadata.df[x,]$study_accession)))
temp <- temp[c("OTU.ID", "Sample", "study_accession", "Relative_abundance", "RepSeq")]
write.csv(x = temp, file = "Result_tables/combined/combined_most_abundant_unassigned.csv", row.names = F, quote = F)

# Combine matrices. These need to be melted together and re-spread

combine_matrices <- function(base_location, mypattern){
  myfiles <- list.files(base_location, pattern = mypattern, recursive = T, full.names = T)
  my_data_frame <- NULL
  for (myfile in myfiles){
    if (is.null(my_data_frame)){
      temp <- read.csv(myfile)
      my_data_frame <- melt(temp)
    } else{
      my_data_frame <- rbind(my_data_frame, melt(read.csv(myfile)))
    }
  }
  my_data_frame <- my_data_frame %>% spread(variable,value,fill = 0)
  my_data_frame
}
# temp <- combine_matrices("Result_tables", "P.*_Phylum_counts_rarefied.csv")


write.csv(x = combine_matrices("Result_tables","P.*_Genus_counts.csv"), 
          file = "Result_tables/combined/combined_Genus_counts.csv", row.names = F, quote = F)

write.csv(x = combine_matrices("Result_tables","P.*_OTU_counts.csv"), 
          file = "Result_tables/combined/combined_OTU_counts.csv", row.names = F, quote = F)

write.csv(x = combine_matrices("Result_tables","P.*_Genus_relative_abundances.csv"), 
          file = "Result_tables/combined/combined_Genus_relative_abundances.csv", row.names = F, quote = F)

write.csv(x = combine_matrices("Result_tables","P.*_Family_relative_abundances.csv"), 
          file = "Result_tables/combined/combined_Family_relative_abundances.csv", row.names = F, quote = F)

write.csv(x = combine_matrices("Result_tables","P.*_Class_relative_abundances.csv"), 
          file = "Result_tables/combined/combined_Class_relative_abundances.csv", row.names = F, quote = F)

write.csv(x = combine_matrices("Result_tables","P.*_Phylum_relative_abundances.csv"), 
          file = "Result_tables/combined/combined_Phylum_relative_abundances.csv", row.names = F, quote = F)

# write.csv(x = combine_matrices("Result_tables","P.*_Genus_counts_rarefied.csv"), 
#           file = "Result_tables/combined_Genus_counts_rarefied.csv", row.names = F, quote = F)
# 
# write.csv(x = combine_matrices("Result_tables","P.*_OTU_counts_rarefied.csv"), 
#           file = "Result_tables/combined_OTU_counts_rarefied.csv", row.names = F, quote = F)
# 
# write.csv(x = combine_matrices("Result_tables","P.*_Genus_relative_abundances_rarefied.csv"), 
#           file = "Result_tables/combined_Genus_relative_abundances_rarefied.csv", row.names = F, quote = F)
# 
# write.csv(x = combine_matrices("Result_tables","P.*_Family_relative_abundances_rarefied.csv"), 
#           file = "Result_tables/combined_Family_relative_abundances_rarefied.csv", row.names = F, quote = F)
# 
# write.csv(x = combine_matrices("Result_tables","P.*_Class_relative_abundances_rarefied.csv"), 
#           file = "Result_tables/combined_Class_relative_abundances_rarefied.csv", row.names = F, quote = F)
# 
# write.csv(x = combine_matrices("Result_tables","P.*_Phylum_relative_abundances_rarefied.csv"), 
#           file = "Result_tables/combined_Phylum_relative_abundances_rarefied.csv", row.names = F, quote = F)

# ------------------------------------------------------------
# Combine fasta files
my_files <- list.files("Result_tables",pattern = "P.*most_abundant_unassigned_features.fasta", recursive = T,include.dirs = T,full.names = T)
# my_files <- my_files[grepl("PRJEB30328", my_files)]
combined_fasta_info.df <- data.frame("OTU.ID" = character(), "RepSeq" = character())
for (fastafile in my_files){
  # print(fastafile)
  temp <- read.fasta(fastafile,as.string = T,forceDNAtolower = F)
  for (seq_name in names(temp)){
    combined_fasta_info.df <- unique(rbind(combined_fasta_info.df, data.frame("OTU.ID" = seq_name,
                                                                       "RepSeq" = as.character(temp[seq_name]))))
  }
}
# combined_fasta_info.df[combined_fasta_info.df$OTU.ID == "a790bd7b8fea7a2ae52dc0c6aba93e51",]
# Write fasta file
write.fasta(sequences = as.list(combined_fasta_info.df$RepSeq),open = "w", 
            names = as.character(combined_fasta_info.df$OTU.ID),
            file.out = paste0("Result_tables/combined/combined_most_abundant_unassigned_features.fasta"))
# ------------------------------------------------------------


