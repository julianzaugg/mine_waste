# Diversity calculations (Shannon, Chao1, Simpson) for each sample.
# Significance tests of each discrete group comparing diversity indices.
# Generate boxplots for discrete 
# Comparison of continuous variables vs diversity indices


library(vegan)
library(reshape2)
library(dplyr)
library(ggplot2)
library(FSA)
library(phyloseq)
library(nlme)

source("Code/helper_functions.R")


common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
  panel.background = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white", size = 1),
  strip.text = element_text(size = 6),
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

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/combined/other/combined_processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

otu.m <- as.matrix(read.csv("Result_tables/combined/count_tables/combined_OTU_counts.csv", row.names =  1))

# Create the rarefied matrix
otu_rare.m <- t(rrarefy(t(otu.m[,colSums(otu.m) >= 5000]), 5000))

# create phyloseq object
otu_rare_phyloseq <- otu_table(otu_rare.m, taxa_are_rows=TRUE)


# Calculate number of unique features/OTUs per sample
temp <- otu_rare.m
temp[temp > 0] <- 1
temp <- melt(colSums(temp))
metadata.df$Number_of_features <- temp[rownames(metadata.df), ,]

# Estimate alpha diversities
otu_rare_alpha.df <- estimate_richness(otu_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_rare_alpha.df <- otu_rare_alpha.df[rownames(metadata.df),]

# Combine the metadata and the diversity metrics into a single dataframe
full=cbind(metadata.df, otu_rare_alpha.df)


generate_diversity_boxplot <- function(mydata, variable, metric, variable_colours_available = T){
  internal_data.df <- mydata[!is.na(mydata[variable]),]
  variable_values <- factor(as.character(unique(internal_data.df[[variable]])))
  if (variable_colours_available == T){
    color_col_name <- paste0(variable, "_colour")
    variable_colours <- setNames(as.character(unique(internal_data.df[[color_col_name]])), as.character(unique(internal_data.df[[variable]])))
  } else{
    variable_colours <- setNames(my_colour_palette_206_distinct[1:length(variable_values)], variable_values)  
  }
  myplot <- ggplot(internal_data.df, aes(x = get(variable), y = get(metric))) +
    geom_boxplot(outlier.shape = NA, aes(fill = get(variable))) +
    scale_fill_manual(values = variable_colours, name = variable) +
    geom_jitter(size=0.5, width = 0.10, height=0) +
    guides(fill=FALSE) +
    # scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5,.5)) +
    xlab("") +
    ylab(metric)  +
    common_theme +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = 0.5),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "white", colour = "white", size = 1),
          strip.text = element_text(size = 6),
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
          plot.title = element_text(size = 6))
  myplot
}

generate_diversity_boxplot(full, "Commodity", "Chao1",variable_colours_available = T)
