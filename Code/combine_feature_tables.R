library(dplyr)
library(tidyr)
library(reshape2)


# ------------------------------------------------------------------------------------------------------
# Load and process the OTU tables into a single one

# Set the working directory
#setwd("/Users/julianzaugg/Desktop/ACE/major_projects/mine_waste/analysis/")
#mydir <- "data/feature_stats_AP19_separate_downsampled"
mydir <- "/srv/projects1/mine_waste/8_acepipe/feature_statistics"
myfiles <- list.files(mydir)
my_data_frame <- NULL
for (myfile in myfiles){
  if (is.null(my_data_frame)){
    mydata.df <- read.csv(file = paste0(mydir,"/",myfile), header =T)
    temp <- melt(mydata.df[names(mydata.df)[!names(mydata.df) %in% c("Frequency", "Confidence","RepSeq")]])
    temp <- temp[temp$value != 0,]
    my_data_frame <- temp
  }else{
    mydata.df <- read.csv(file = paste0(mydir,"/",myfile), header =T)
    temp <- melt(mydata.df[names(mydata.df)[!names(mydata.df) %in% c("Frequency", "Confidence","RepSeq")]])
    temp <- temp[temp$value != 0,]
    # my_data_frame <- full_join(my_data_frame, new.df, by = "X.OTU.ID")
    my_data_frame <- rbind(my_data_frame, temp)
  }
}

# Requires a lot of memory. May need to be run on server.
my_data_frame_spread <- my_data_frame %>% spread(variable,value,fill = 0)

#write.csv(x = my_data_frame_spread, file = "Data/feature_statistics.csv", quote = F, row.names = F)
write.csv(x = my_data_frame_spread, file = "/srv/projects1/mine_waste/8_acepipe/feature_statistics.csv", quote = F, row.names = F)