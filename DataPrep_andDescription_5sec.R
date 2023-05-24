
####--------------------------------------------------------------------#
# 
# Loading accelerometer data and prep for features
#  Predictive features are generated in accordance with: 
#     Khataeipour SJ et al. BMJ Open Sport & Exercise Medicine 2022;8:e001242.doi: 10.1136/bmjsem-2021-001242
#  Another set of predictive featured,not used in this particular study was generated in this script, following the paper 
#     Barua A et al. Biosensors (Basel). 2022 Jul 21;12(7):549. doi: 10.3390/bios12070549.
#
#
# The original data consist of frequency of 30-50Hz. THey are strandardized to 30Hz, and aggregated to 5 second time window. 
# The resulting file used to fit and test BART and random forest is called 
#   bmjFeatures_reSampled_5sec.csv 
# and the original data is named: 
#   accel_data_no_id.csv
#
# Both data are found in the website:     
#   https://doi.org/10.7910/DVN/LXVZRC
#
# Several other data are generated in this script but can be ignored
# THis script was originally run in multicore unix machine. If running in windows, parallel code (mclapply) can be switched to lappy
#--------------------------------------------------------------------#

rm(list=ls())
setwd("~/scratch/BART_test/")
set.seed(1)

# Glimpse file - works in unix  
# System("zcat ~/data/accel_data.zip | wc -l") # peek
# System("zcat  ~/data/accel_data.zip | head -n 1")

library(data.table)
library(R.utils)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(activityCounts)
library(magrittr)
library(zoo)
library(e1071)
library(entropy)
library(progress)
library(parallelly)
library(imputeTS)
library(parallel)
library(kableExtra)
library(tictoc)
library(facetscales)

# Some run configuration
bool_plot = F # plot figures or not (some images are large)
bool_saveplot = F # save image to file 
bool_load = T # Create data or load already processed data, the former takes long due due to data size 
bool_BMJFeatureGenerate = T 

# Functions 
source("functions/generateFeatures.R")# Function to generate features from the BMJ paper
source("functions/generateBiosensorFeatures.R")  # Function to generate features from the second (Biosensor) paper
source("functions/funcPlot.R") # Function to make descriptive plots 

# Varaibles to subset by location, people, and class     
devLocation = "all" # For now, all 3 device locations - pocket, arm, and backpack are loaded  


# people with raw adxis values far greater than plus or minus 50 i.e. outliers in acclerometer magnitude 
removeVecLargeGravityID <- c(111, 129, 130, 138, 139) 

# People with serious autocorrelation in data, malfunction of sensor?
removeAutocorID <- c(146, 114)

# People where class label were inaccurate due to some error in lab equipment to measure METs
removeVecID2 <- c(112, 121, 122, 132) # From the BMJ study # Participants 112 , 121, 122 , and 132 are not properly classified so shoudl not be used 

# Participant IDs with unusual patterns of accelerometer signals, included in the analysis 
removeVecID1 <- c(112, 127, 129, 130, 131, 138) 

# Temporal resolution for raw data (number of data per second), to be resampeld at this frequency.  
Hz = 30
windowSec = 5 # window to generate features from the BMJ study 

cat("system.name")
cat(Sys.info()[4])
cat("\n")

# Load data - these data are stored in unit machine (hence using unzip) - rewrite if using windows 
if(grep("gra", Sys.info()[4])){ #if loading from Graham cluster 
  dAll <- fread(cmd = 'unzip -cq ~/scratch/BART_test/data/accel_data.zip', 
                select=c("record_time", 
                         "x_axis", 
                         "y_axis", 
                         "z_axis", 
                         "participant_id", 
                         "wear_location",
                         "trimmed_activity"
                )
  )# requires unzip on path
}else{ # if loading from Cedar cluster
  #dAll <- fread(cmd = 'unzip -cq ~/data/accel_data.zip', 
  dAll <- fread(cmd = 'unzip -cq ~/projects/def-dfuller/interact/exchanges/Accel_Labelled/accel_data.zip', 
                select=c("record_time", 
                         "x_axis", 
                         "y_axis", 
                         "z_axis", 
                         "participant_id", 
                         "wear_location",
                         "trimmed_activity"
                )
  )# requires unzip on path 
}



#### processing data 
# Rename variables            
dAll <- dAll[, .(record_time,
                 x = x_axis, 
                 y = y_axis, 
                 z = z_axis, 
                 id = participant_id, 
                 loc = wear_location, 
                 class=trimmed_activity)]

# remove some people defined above
#dAll <- dAll[dAll[, !(id %in% removeVecID1)]]
dAll <- dAll[dAll[, !(id %in% removeVecID2)]]
dAll <- dAll[dAll[, !(id %in% removeVecLargeGravityID)]]
dAll <- dAll[dAll[, !(id %in% removeAutocorID)]]
rm(list = ls()[grep("remove", ls())])

# Remove transit Class  
dAll <- dAll[class != "transit"] # rm transit

# Add a variable coding time ordering within each second
dAll[, timeVar := .(seq_len(.N)), by = c("record_time", "id", "loc")]





### Lead variables, first order and second order lag for each axis to differenciate in later plots -----
ori <- c("x", "y", "z") # Vector of variable names for the axis of accelerometer daa 

# Creat lagged (lead) values
diffVar <- paste( ori, "d1", sep="_")
dAll[,  (diffVar) := data.table::shift(.SD, n = 1, fill = NA, type = "lead") - .SD, .SDcols=ori, by = c("id", "loc")]

diffVar2 <- paste(ori, "d2", sep="_")
dAll[,  (diffVar2) := data.table::shift(.SD, n= 1, fill = NA, type =  "lead") - .SD, .SDcols=diffVar, by = c("id", "loc")]

# Creating a vector of column names for X, Y, Z axis of acclerometer data 
rollSDVar30_d1 <- paste(ori, "mvSD_d1_30", sep="_") # MOving SD on first order lagged data 
rollSDVar30_d2 <- paste(ori, "mvSD_d2_30", sep="_") # Moving SD on second order lagged ata 
rollSDVar30 <- paste(ori, "mvSD_30", sep="_") #Moving SD
rollMeanVar30 <- paste(ori, "mean30", sep="_") # Moving mean 


# Generate window statistics mean and SD. Takes long time to compute, so load the pre-calcualted ata
# This computation takes time, so the processed data are saved if bool_load==FALSE. If bool_load == TRUE, load the pre-processed data 
if(!bool_load){
  cat("Creating SD and mean of time-series. Lenghty process")
  ### SD (standard deviation) of first and second order difference 
  dAll[, (rollSDVar30_d2) := frollapply(.SD, n=30, FUN=sd, fill=NA, na.rm= T, align=c("right")), 
       .SDcols=diffVar2, by = c("id", "loc")]
  
  dAll[, (rollSDVar30_d1) := frollapply(.SD, n=30, FUN=sd, fill=NA, na.rm= T, align=c("right")), 
       .SDcols=diffVar, by = c("id", "loc")]
  
  ### Moving SD
  dAll[, (rollSDVar30) := frollapply(.SD, n=30, FUN=sd, fill=NA, na.rm= T, align=c("right")),
       .SDcols=ori, by = c("id", "loc")]
  
  ### Moving M=mean 
  dAll[, (rollMeanVar30) := frollmean(.SD, n=30, fill=NA, na.rm= T, align=c("right"), algo= "fast"), .SDcols=ori, by = c("id", "loc")]
  
  # save for next time 
  saveRDS(dAll, "~/scratch/BART_test/data/dAll_windowStats_raw.rds")
}else{ # load pre-processed data   
  dAll <- readRDS("~/scratch/BART_test/data/dAll_windowStats_raw.rds")
}





# Skip this codes unless generating descriptive plots 
if(1==2){ # not run currently 
  #### Descriptive time-series plots - plots are saved as too large in filesize and dimension
  ### TS plots preperation 
  # select person ID to gauge 
  plotPersonID <- 108
  plotDevLocation <- "hand"
  timePoint <- c(1, 300) # to zeeom 
  
  # Create DF of a person 
  dTs <- dAll[id == plotPersonID & loc == plotDevLocation]
  dTs$seq <- 1:nrow(dTs)
  #cat("Plot activity types to check patterns of signals within seconds")
  sampDataList <- dTs[,.SD[timePoint[1]:timePoint[2]], by = class] %>% split(by =c("class"))
  
  ### Plot zoomed time-series for each class and each axis  (3 axis and 6 classes)------
  # function to plot xyz by class

  
  p <- lapply(sampDataList, funcPlotAxis)
  imagePath = "~/scratch/BART_test/image/descriptive/axisUnsampled/"
  g <-arrangeGrob(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], ncol = 1)
  #g <-arrangeGrob(grobs = p, ncol = 1)
  if(bool_plot){
    ggsave(file = paste(imagePath, "Plot of each Class__PersonID", plotPersonID, plotDevLocation, "timeRange_", timePoint[1], ":", timePoint[2], ".pdf", sep = "_" ), g, height = 20, width = 15)
  }
  
  
  ### Same plots of xyz without separrating by class, entire time series, so it is large  ----- 
  d <- melt(dTs, id.vars = c("seq", "class"), measure.vars = c("x", "y", "z"))
  p <- ggplot(d, aes(y = value, x = seq, color= factor(class))) +
    #geom_line(size=0.3) +
    geom_path(aes(group = 1), size = 0.1) +
    facet_wrap(~variable, scales = "free", ncol = 1) +  
    theme_minimal()
  if(bool_plot){
    ggsave(file = paste(imagePath, "Plot of AllRawPoints_ColorCodedClass", plotPersonID, plotDevLocation, ".pdf", sep = "_" ), p, height = 20, width = 90, limitsize = FALSE)
  }
  
  # Different person, to see lead of class definition occur or not 
  plotPersonID2 <- 113; plotDevLocation <- "hand"
  dTs2 <- dAll[id == plotPersonID2 & loc == plotDevLocation]
  dTs2$seq <- 1:nrow(dTs2)
  d2 <- melt(dTs2, id.vars = c("seq", "class"), measure.vars = c("x", "y", "z"))
  p <- ggplot(d2, aes(y = value, x = seq, color= factor(class))) +
    geom_path(aes(group = 1), size = 0.1) +
    facet_wrap(~variable, scales = "free", ncol = 1) +  
    theme_minimal()
  
  if(bool_plot){
    ggsave(file = paste(imagePath, "Plot of AllRawPoints_ColorCodedClass", plotPersonID2, plotDevLocation, ".pdf", sep = "_" ), p, height = 20, width = 90, limitsize = FALSE)
  }
  
  plotPersonID2 <- 156; plotDevLocation <- "hand"
  dTs2 <- dAll[id == plotPersonID2 & loc == plotDevLocation]
  dTs2$seq <- 1:nrow(dTs2)
  d2 <- melt(dTs2, id.vars = c("seq", "class"), measure.vars = c("x", "y", "z"))
  p <- ggplot(d2, aes(y = value, x = seq, color= factor(class))) +
    geom_path(aes(group = 1), size = 0.1) +
    facet_wrap(~variable, scales = "free", ncol = 1) +  
    theme_minimal()
  if(bool_plot){
    ggsave(file = paste(imagePath, "Plot of AllRawPoints_ColorCodedClass", plotPersonID2, plotDevLocation, ".pdf", sep = "_" ), p, height = 20, width = 90, limitsize = FALSE)
  }
  
  
  ### Same as above, but display and save each axis separately, with 1st and 2nd order delta of axis ------------------  
  d <- melt(dTs, id.vars = c("seq", "class"), measure.vars = c(ori, diffVar, diffVar2))
  d$variable <- factor(d$variable, levels=c(c(ori, diffVar, diffVar2)[c(1, 4, 7)+0], 
                                            c(ori, diffVar, diffVar2)[c(1, 4, 7)+1], 
                                            c(ori, diffVar, diffVar2)[c(1, 4, 7)+2]))
  
  
  p <- funcPlotClassByColor(d[variable %in% c(ori, diffVar, diffVar2)[c(1, 4, 7)+0]])
  p <- p+ggtitle(paste("Raw (not resampled) axis data, with first and second order difference: PersonID", plotPersonID, "Device Location", plotDevLocation, ".pdf", sep = "_" ))
  
  if(bool_plot){
    ggsave(file = paste(imagePath, "X_RawPoints_Deltas_Person_", plotPersonID, plotDevLocation, ".pdf", sep = "_" ), p, height = 30, width = 90, limitsize = FALSE)
  }
  
  
  p <- funcPlotClassByColor(d[variable %in% c(ori, diffVar, diffVar2)[c(1, 4, 7)+1]])
  p <- p+ggtitle(paste("Raw (not resampled) axis data, with first and second order difference: PersonID", plotPersonID, "Device Location", plotDevLocation, ".pdf", sep = "_" ))
  if(bool_plot){
    ggsave(file = paste(imagePath, "Y_RawPoints_Deltas_Person_", plotPersonID, plotDevLocation, ".pdf", sep = "_" ), p, height = 30, width = 90, limitsize = FALSE)
  }
  
  p <- funcPlotClassByColor(d[variable %in% c(ori, diffVar, diffVar2)[c(1, 4, 7)+2]])
  p <- p+ggtitle(paste("Raw (not resampled) axis data, with first and second order difference: PersonID", plotPersonID, "Device Location", plotDevLocation, ".pdf", sep = "_" ))
  if(bool_plot){
    ggsave(file = paste(imagePath, "Z_RawPoints_Deltas_Person_", plotPersonID, plotDevLocation, ".pdf", sep = "_" ), p, height = 30, width = 90, limitsize = FALSE)
  }
  
  
  
  ### show spike in y axis on hands, person 108, spike occur about 30 data points width so not very short up and down ------  
  varSet <- c(ori, diffVar, diffVar2, rollSDVar30, rollSDVar30_d1, rollSDVar30_d2)
  varSetIndex <- c(1,4,7,10,13,16)
  
  d <- melt(dTs, id.vars = c("seq", "class"), measure.vars = varSet)
  d$variable <- factor(d$variable, 
                       levels=c(varSet[varSetIndex+0], 
                                varSet[varSetIndex+1], 
                                varSet[varSetIndex+2])
  )
  
  
  ### Plot all Features 
  p <- funcPlotClassByColor(d[variable %in% varSet[varSetIndex+1]])
  p <- p+ggtitle(paste("Raw (not resampled) axis data, with first and second order difference and the SD of raw data, first and second order differenciated data: PersonID", plotPersonID, "Device Location", plotDevLocation, ".pdf", sep = "_" ))
  
  if(bool_plot){
    imagePath = "~/scratch/BART_test/image/descriptive/axisUnsampled/"
    ggsave(file = paste(imagePath, "Y_RawPoints_DeltasAndWindowStats_Person_", plotPersonID, plotDevLocation, ".pdf", sep = "_" ), p, height = 80, width = 90, limitsize = FALSE)
  }
  
  
  
  ### magnify 
  ### More difficult points 
  #### Magnify 
  timePoint = c(10800:16300)
  p <- funcPlotClassByColor(d[seq %in% timePoint][variable %in% varSet[varSetIndex+1]])
  p <- p+ggtitle(paste("Raw (not resampled) axis data, with first and second order difference: PersonID", plotPersonID, "Device Location", plotDevLocation, "TimePoint", min(timePoint), max(timePoint),  ".pdf", sep = "_" ))
  p <- p + scale_color_manual(values=c("black"))
  
  if(bool_plot){
    imagePath = "~/scratch/BART_test/image/descriptive/axisUnsampled/"
    ggsave(file = paste(imagePath, "Y_TimePoint", min(timePoint), max(timePoint), "RawPoints_DeltasAndWindowStats_Person_", plotPersonID, plotDevLocation,  ".pdf", sep = "_" ), p, height = 80, width = 90, limitsize = FALSE)
  }
  
  ### magnify another point 
  timePoint = c(3600:4800)
  p <- funcPlotClassByColor(d[seq %in% timePoint][variable %in% varSet[varSetIndex+1]])
  p <- p+ggtitle(paste("Raw (not resampled) axis data, with first and second order difference: PersonID", plotPersonID, "Device Location", plotDevLocation, "TimePoint", min(timePoint), max(timePoint),  ".pdf", sep = "_" ))
  p <- p + scale_color_manual(values=c("black"))
  
  if(bool_plot){
    imagePath = "~/scratch/BART_test/image/descriptive/axisUnsampled/"
    ggsave(file = paste(imagePath, "Y_TimePoint", min(timePoint), max(timePoint), "RawPoints_DeltasAndWindowStats_Person_", plotPersonID, plotDevLocation,  ".pdf", sep = "_" ), p, height = 80, width = 90, limitsize = FALSE)
  }
}





#### Removal of anomaliy signal (perfectly linear increase or decrease of signals) that should be discarded 
# Window sum of D-2 
ori <- c("x", "y", "z")
plotDevLocation <- "hand"
plotPersonID <- 108

# variable names to calcualte window statistics 
diffVar2 <- paste(ori, "d2", sep="_")
rollSum_window30_d2 <- paste(ori, "Sum_d2_30", sep="_")

dAll[, (rollSum_window30_d2) := 
       frollsum(abs(.SD), n=30, algo=c("fast") ,fill=NA, na.rm= T,align=c("right")), 
     .SDcols=diffVar2, by = c("id", "loc")]


if(1==2){ # Vizualization of signals to be removed, not run currently 
  # Chenage range of delta D and double-delta D (first and second order difference) to remove, 
  # for a person 108, anomaly is around 15000
  dTs <- dAll[id == plotPersonID & loc == plotDevLocation]
  dTs$seq <- 1:nrow(dTs)
  rollSum_window30_d2_LOG <- paste(rollSum_window30_d2, "LOG", sep = "_")
  
  # Make it natural log for interpretation / viz 
  dTs[, (rollSum_window30_d2_LOG) := log(.SD), .SDcols = c(rollSum_window30_d2)]
  
  # Time range to visualize as an example 
  timePoint = c(15000:15800)
  
  # Vector for var names to display 
  varSet <- c(ori, diffVar, diffVar2, rollSDVar30, rollSDVar30_d1, rollSDVar30_d2)
  varSet <- c(varSet, rollSum_window30_d2, rollSum_window30_d2_LOG)
  varSetIndex <- which(varSet %in% c("y", "y_Sum_d2_30", "y_Sum_d2_30_LOG", "y_d2"))
  
  # distribution of log transformed 2nd-order Y
  cat("Distribution of second-order Y to decide cut-off value to remove record_time (seconds) that have anomaly ")
  hist(dTs$y_Sum_d2_30_LOG, main = "Distribution of standard deviatio of second order difference in raw data, in log")
  hist(dTs$y_Sum_d2_30, main = "Distribution of standard deviatio of second order difference in raw data, non-log, small values indicate linear trend in signal", xlim = c(0, 100))
  
  d <- melt(dTs, id.vars = c("seq", "class"), measure.vars = varSet[varSetIndex])
  d$variable <- factor(d$variable, levels=c(varSet[varSetIndex]))
  
 
  
  ### Plot to see cutoff value of SD and log SD to removing lienar raw axis signal , 
  ### Plot function, this one gives pre-set axis scale for Y
  funcCustomScale <- function(df){        
    scales_y <- list(
      `y`           = scale_y_continuous(limits = c(-20, 20)),
      `y_Sum_d2_30` = scale_y_continuous(limits = c(10^-15, 10^-1)), 
      `y_Sum_d2_30_LOG`  = scale_y_continuous(limits = c(-30, -0)),
      `y_d2`        = scale_y_continuous(limits = c(-30, 30))
    )
    
    p <- ggplot(df, aes(y = value, x = seq, color= factor(class))) +
      geom_path(aes(group = 1), size = 0.1) +
      #https://stackoverflow.com/questions/66658012/how-to-change-the-limits-and-divisions-of-the-axes-of-each-facet-using-ggplot2
      #facet_grid_sc(scales = list(y = scales_y), rows=vars(variable)) +  
      facet_grid_sc(scales ="free" , rows=vars(variable)) +  
      #facet_wrap(~ variable, ncol = 1) +  
      theme_minimal() + 
      theme(strip.text.x = element_text(size = 60)) + 
      theme(text = element_text(size = 60)) +  
      ylab("Raw value") + xlab("Raw (not resampled) time-series ") + 
      scale_x_continuous(breaks = round(seq(min(d$seq), max(d$seq), by = 5000),0)) + 
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+ # graph order lines 
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
    return(p)
  }
  
  ### title and save 
  p <- funcCustomScale(d)
  p <- p+ggtitle(paste("Raw (not resampled) axis data, second order difference and window sum of second order difference and its log: PersonID", plotPersonID, "Device Location", plotDevLocation, "TimePoint", min(timePoint), max(timePoint),  ".pdf", sep = "_" ))
  #p <- p + scale_color_manual(values=c("black"))
  
  if(bool_plot){
    imagePath = "~/scratch/BART_test/image/descriptive/axisUnsampled/"
    ggsave(file = paste(imagePath, "Y_WIndowD2", min(timePoint), max(timePoint), "RawPoints_ WindowDelta_Person_", plotPersonID, plotDevLocation,  ".pdf", sep = "_" ), p, height = 40, width = 90, limitsize = FALSE)
  }
  

  #Flag to remove these values based on the window sum of SD of second order difference, At the level of Second 
  dAll[, removeD2Low := ifelse(y_Sum_d2_30 < 10^-5 ,1, 0) ][
    , `:=` (sumRemoveInd = sum(removeD2Low)), by = c("id", "loc", "record_time")][
      , indTrimSecond := ifelse(sumRemoveInd != 0 ,1, 0)      
      ]
  
  # Newly added codes to make sure to remove only data that consecutive 30 points iwth linear increase (moreve by second)
  plotData <- dAll[1:800, ]
  plotData$timeID <- 1:nrow(plotData)
  p1 <- ggplot(data=plotData, aes(x=timeID, y=y_Sum_d2_30)) + geom_point()
  p2 <- ggplot(data=plotData, aes(x=timeID, y=removeD2Low)) + geom_point()
  p3 <- ggplot(data=plotData, aes(x=timeID, y=indTrimSecond)) + geom_point()
  p4 <- ggplot(data=plotData, aes(x=timeID, y=sumRemoveInd)) + geom_point()
  p5 <- ggplot(data=plotData, aes(x=timeID, y=y)) + geom_line()
  
  grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
  
  cat("Record-time (seconds) to remove from the ad-hoc filter above")
  table(dAll$indTrimSecond) #remove   1804699 out of 14558961  
  dAll$removeD2Low <- NULL
  dAll$sumRemoveInd <- NULL 
  
  
  # visually verify data to remove 
  dTs <- dAll[id == plotPersonID & loc == plotDevLocation]
  dTs$seq <- 1:nrow(dTs)
  
  d <- data.frame(y = dTs$y, remove = as.factor(dTs$indTrimSecond), x = 1:nrow(dTs))
  p <- ggplot(d, aes(y = y, x = x, color=remove)) + 
    geom_path(aes(group = 1), size = 0.1) + 
    xlab("Raw time points") + 
    scale_color_manual(values = c("grey", "black")) + 
    theme_minimal() + 
    theme(strip.text.x = element_text(size = 60)) + 
    theme(text = element_text(size = 60)) +  
    ylab("Raw value") + xlab("Raw (not resampled) time-series ") 
  
  p <- p+ggtitle(paste("Raw (not resampled) axis data, data to be removed (black lines)", plotPersonID, "Device Location", plotDevLocation, "TimePoint", min(timePoint), max(timePoint),  ".pdf", sep = "_" ))
  
  
  if(bool_plot){
    imagePath = "~/scratch/BART_test/image/descriptive/axisUnsampled/"
    ggsave(file = paste(imagePath, "Y_TimePoint", min(timePoint), max(timePoint), "RawPoints_IndicatorTobeRemoved", plotPersonID, plotDevLocation,  ".pdf", sep = "_" ), p, height = 30, width = 300, limitsize = FALSE)
  }
  
  
  # Remove and compare sample sizes by location and class 
  cat("count of raw data by device location before removing seconds with anomaly ") 
  knitr::kable(dAll[, .N, by=loc], 
               format.args = list(big.mark = ",",
                                  scientific = FALSE)) %>%
    kable_styling(full_width = F)
  
  # Finally, remove records (record_time, seconds) with wired signals
  dAll <- dAll[indTrimSecond ==0, ]
  
  cat("count by device location after removing seconds with anomaly ") 
  knitr::kable(dAll[, .N, by=loc], 
               format.args = list(big.mark = ",",
                                  scientific = FALSE)) %>%
    kable_styling(full_width = F)
  
  
  # remove window variables (mean and SD or xyz that were used to remove from original raw time-series due to abnormal time signals )
  dAll[, grep("mvSD", names(dAll)) := NULL] 
  dAll[, grep("mean", names(dAll)) := NULL] 
  dAll[, grep("Sum_d2", names(dAll)) := NULL] 
  dAll[, grep("Remove", names(dAll)) := NULL] 
  dAll[, grep("remove", names(dAll)) := NULL] 
  
  #### Resampling to 30Hz from data with 30-40Hz, a couple different approach 
  #Drawing original data poitns from each seconds
  #Calculting mean from 30tile per sec (so aggregation, rather than resample)
  ### Function that will select N values and preserve the order
  mySample2 <- function(v, n){v[sample(c(rep(TRUE, n), rep(FALSE, length(v) - n)))]} #https://stackoverflow.com/questions/63565123/get-an-uniform-distributed-sample-but-keep-the-order
  mySample <- function(values,N){
    size <- length(values)
    values[sapply(1:size, function(i){
      select <- as.logical(rbinom(1,1,N/(size+1-i)))
      if(select) N <<- N - 1
      select
    })]
  }
}        













# Here data is standardized to 30Hz (30 points per sec) two appraoches 
# Resampling data to 30Hz make aggregated axis at 30th of sec ##################################   

dAll[, originalHz :=.N, by = c("record_time", "id", "loc", "class")]

if(!bool_load){
  funcQuant <- function(x)findInterval(x$timeVar, c(-Inf, quantile(x$timeVar, probs=c(seq(30)/31))), Inf)
  dAll[, quant := funcQuant(.SD), by = c("record_time", "id", "loc")] # 
  vars <- c("x", "y", "z", "timeVar")
  dSampAgg <- dAll[ , lapply(.SD, mean), by = c("record_time", "id", "loc", "quant", "class", "originalHz"), .SDcols = vars]
  
  #dSampAgg <- dSampAgg[, head(.SD, 30), by = c("record_time", "id", "loc")]
  dSampAgg = dSampAgg[dSampAgg[, .I[1:30], by = c("record_time", "id", "loc")]$V1]
  dSampAgg <- dSampAgg[order(id, loc, record_time, timeVar)]
  
  #class <- dAll[, originalHz :=.N, by = c("record_time", "id", "loc", "class")]
  
  saveRDS(dSampAgg, "~/scratch/BART_test/data/dSampAgg.rds")
  #save.image("dAll_And_DAgg.RData")
}else{
  dSampAgg <- readRDS("~/scratch/BART_test/data/dSampAgg.rds")
  #load("dAll_And_DAgg.RData")  
}

# Check no lost records due to resampling???
#nrow(dAll[, .N, by = c("record_time", "id", "loc")])*30
#nrow(dSampAgg)



# 30Hz data (2) # Resample by drawing data 30 points at each sec ################### 
#takes long #dSamp =  dAll[,.SD[sample(.N, min(30,.N))], by = c("record_time", "id", "loc")]
# much faster
dSamp <- dAll[dAll[ , .I[sample(.N,30)] , by = c("record_time", "id", "loc", "originalHz")]$V1] 
#https://stackoverflow.com/questions/33887083/from-data-table-randomly-select-one-row-per-group
# order data by time var 
dSamp <- dSamp[order(id, loc, record_time, timeVar)]


# A few hundred duplicated records in the raw data - remove?   
dup <- dAll[duplicated(dAll[,c("id", "loc", "x", "y", "z")]) | duplicated(dAll[,c("id", "loc", "x", "y", "z")], fromLast = T)]

# make sure that resampling was not done on the same data point 
dSamp[duplicated(dSamp[, c("x", "y", "z")]) | duplicated(dSamp[, c("x", "y", "z")], fromLast = TRUE)]




# remove some vars 
dSamp[, c("timeVar", "quant") := NULL] 
dSampAgg[, c("timeVar", "quant") := NULL] 

# ...and new time variable   
dSamp <- dSamp[,timeVarNew := 1:.N, by = c("id", "loc", "record_time")] 
dSampAgg <- dSampAgg[, timeVarNew := 1:.N, by = c("id", "loc", "record_time")] 







### BMJ study features , 5 HZ so aggregated by 5seconds -------
### Features from the BMJ paper, for dSamp data ------------------------------------------------
# Subset by location, if "all", use all data 
if(devLocation != "all"){
  dSamp <- dSamp[loc == devLocation]
}


# Multicore function is used. Remove this line if using non-unix multicore machine. 
availCore <- parallelly::availableCores()
numCore <- 30

# ID for each time window
dSampAgg[, windowID := rep(seq_len(ceiling(.N/(Hz*windowSec))), each=(Hz*windowSec), length.out=.N),by=c("id", "loc")]
dSamp[, windowID := rep(seq_len(ceiling(.N/(Hz*windowSec))), each=(Hz*windowSec), length.out=.N),by=c("id", "loc")]



### Create features for each person and location - again, run in unix multicore machine. Rewrite to lapply from mclapply if necessary 
if(.Platform$OS.type == "unix"  & bool_BMJFeatureGenerate){
  dSampSplit <- split(dSampAgg, by = c("id", "loc"))
  # function to create featues, "GenerateFeatures" is borrowed from 
        # https://github.com/walkabillylab/Smartphone_accelerometers-Pocket_location/blob/88a3bdc23a3217e4772d62888bc38260af198403/Misc_codes/generate-features.R
  dSampSplitFeatures <- parallel::mclapply(dSampSplit, function(x) GenerateFeatures(data.frame(x), window_size_sec = windowSec), mc.cores = numCore)
  dSampFeatures <- bind_rows(dSampSplitFeatures, .id = "id.label")
  setDT(dSampFeatures)
  
  # num records in each list 
  unlist(lapply(dSampSplitFeatures, nrow))
  
  #save.image()
  # Create label for id and device location 
  #dSampFeatures$time <- dSampFeatures$...1 
  dSampFeatures <- dSampFeatures %>%  
    separate(id.label, c("id", "loc")) %>% 
    mutate(id = as.numeric(id))
  
  ####  Class label ---------------------------------------
  d <- dSampAgg[,c("class", "id", "loc", "record_time", "windowID")]
  
  class <-  dSampAgg  %>%  group_by(id, loc, windowID) %>%
    summarise_all(first) 
  
  # number of counts (frequency per sec) in each window 
  c <- dSampAgg  %>%  group_by(id, loc, windowID) %>%
    summarise(countInWindow = n()) 
  
  class$countInWindow <- c$countInWindow
  
  #sum(
  #  data.frame(class) %>% filter(countInWindow==150) %>% select(id)!= 
  #    setDF(dSampFeatures) %>% select(id)
  #)
  
  #sum(
  #  data.frame(class) %>% filter(countInWindow==150) %>% select(loc)!= 
  #    setDF(dSampFeatures) %>% select(loc)
  #)
  
  # Attach class to generated features 
  if(nrow(dSampFeatures) == nrow(class[class$countInWindow ==150, ])){
    dSampFeatures <- cbind(dSampFeatures, class[class$countInWindow ==150, c("windowID", "record_time", "class")])
  }
  
  dSampFeatures$class <- factor(dSampFeatures$class, levels=c("Lying", "Sitting", "Self Pace walk", "Running 3 METs", "Running 5 METs", "Running 7 METs"))
  
  # Save data 
  saveRDS(dSampFeatures, "~/scratch/BART_test/data/bmjFeatures_reSampledAgg_5sec.rds")
  rm(dSampSplit)
}




### For re-sampled data 1  
if(.Platform$OS.type == "unix" & bool_BMJFeatureGenerate ){
  dSampSplit <- split(dSamp, by = c("id", "loc"))
  dSampSplitFeatures <- parallel::mclapply(dSampSplit, function(x) GenerateFeatures(data.frame(x), window_size_sec = windowSec), mc.cores = numCore)
  dSampFeatures <- bind_rows(dSampSplitFeatures, .id = "id.label")
  setDT(dSampFeatures)
  
  # num records in each list 
  unlist(lapply(dSampSplitFeatures, nrow))
  
  #save.image()
  # Create label for id and device location 
  #dSampFeatures$time <- dSampFeatures$...1 
  dSampFeatures <- dSampFeatures %>%  
    separate(id.label, c("id", "loc")) %>% 
    mutate(id = as.numeric(id))
  
  ####  Class label ---------------------------------------
  d <- dSampAgg[,c("class", "id", "loc", "record_time", "windowID")]
  
  class <-  dSampAgg  %>%  group_by(id, loc, windowID) %>%
    summarise_all(first) 
  
  # number of counts (frequency per sec) in each window 
  c <- dSampAgg  %>%  group_by(id, loc, windowID) %>%
    summarise(countInWindow = n()) 
  
  class$countInWindow <- c$countInWindow
  
  #sum(
  #  data.frame(class) %>% filter(countInWindow==150) %>% select(id)!= 
  #    setDF(dSampFeatures) %>% select(id)
  #)
  
  #sum(
  #  data.frame(class) %>% filter(countInWindow==150) %>% select(loc)!= 
  #    setDF(dSampFeatures) %>% select(loc)
  #)
  
  # Attach class to generated features 
  if(nrow(dSampFeatures) == nrow(class[class$countInWindow ==150, ])){
    dSampFeatures <- cbind(dSampFeatures, class[class$countInWindow ==150, c("windowID", "record_time", "class")])
  }
  
  # Save data 
  saveRDS(dSampFeatures, "~/scratch/BART_test/data/bmjFeatures_reSampled_5sec.rds")
  rm(dSampSplit)
}



##### Plot all ID time-series plots of xyz and features, to be imported into another output
### Plot loop for raw xyz-------------

if(bool_plot){
  pXYZ <- funcPlotListTime(d = dSamp, 
                           title = "30Hz_RawXYZ_PersonID_",
                           target= c("x", "y", "z"), 
                           devLocation = "pock")
  
  # Save multiple images 
  imagePath <- "~/scratch/BART_test/image/descriptive/xyz/"
  # Plot loop for raw Biosensors first 4 features  
  imagePath <- "~/scratch/BART_test/image/descriptive/angleNorm/"
  pNormAngle <- funcPlotListTime(d = dSamp, 
                                 title = "30Hz_NormAngle_featuresFromBiosensors_PersonID_",
                                 target= c("en", "en_d1", "en_d2", "ang1", "ang2"), 
                                 devLocation = "pock")
}



if(bool_plot){
  availCore <- parallelly::availableCores()
  mcCore <- min(10, availCore)
  mclapply(pXYZ, funcPlotSave, imagePath, mc.cores = mcCore)
  mclapply(pNormAngle, funcPlotSave, imagePath, mc.cores = mcCore, h= 7.5, w = 9)
}




### Biosensor plot - not used in this study , no need to generate these data 
bool_biosensor = FALSE
if(bool_biosensor){
  dSampB <- funcVectorDotFeatures(dSamp, removeDiff = FALSE, deltaGenerate = FALSE)
  dSampAggB <- funcVectorDotFeatures(dSampAgg, removeDiff = FALSE, deltaGenerate = TRUE)
  

  #Missing records due to lead function used to generate norm and angles, and due to invalid denominator (zero values) to generate angle features"
  #sapply(dSampB, function(x) sum(is.na(x)))
  #sapply(dSampAggB, function(x) sum(is.na(x)))
  
  # remove two tailing data points that have NA
  dSampB <- dSampB[!dSampB[, last(.I), by = c("id", "loc")]$V1]
  dSampB <- dSampB[!dSampB[, last(.I), by = c("id", "loc")]$V1]
  dSampAggB <- dSampAggB[!dSampAggB[, last(.I), by = c("id", "loc")]$V1]
  dSampAggB <- dSampAggB[!dSampAggB[, last(.I), by = c("id", "loc")]$V1]
  
  # Remove the records with angle features that have invalid values (zero denominators) for now
  dSampB <- dSampB[complete.cases(dSampB)]
  dSampAggB <- dSampAggB[complete.cases(dSampAggB)]
  
  # Remove some variables 
  dSampB[, grep("lead", names(dSampB)) := NULL] 
  dSampAggB[, grep("lead", names(dSampAggB)) := NULL] 
  dSampB[, (diffVar) := NULL] 
  dSampB[, (diffVar2) := NULL] 
  dSampAggB[, (diffVar) := NULL] 
  dSampAggB[, (diffVar2) := NULL] 
  
  dSampAggB$ang2 <- NULL
  dSampB$ang2 <- NULL
  
  
  #saveRDS(dSamp, "~/scratch/BART_test/data/sensorFeaturesResampled.rds")
  #saveRDS(dSampAgg, "~/scratch/BART_test/data/sensorFeaturesResampledAgg.rds")
  
  # interporation of ang 1 and 2 varaibles - this part is not complete, as only ~2000 peopel are missing, but may need to be addressed later 
  #dSampFeatures$ang1 <- na_interpolation(dSampFeatures$ang1)
  #dSampFeatures$ang2 <- na_interpolation(dSampFeatures$ang2)
  #dSamp$ang1 <- na_interpolation(dSamp$ang1)
  #dSamp$ang2 <- na_interpolation(dSamp$ang2)

}



####Correlation plots, sample of data points from pocket
source("functions/corrPlot.R")
dSamp <- readRDS("~/scratch/BART_test/data/bmjFeatures_reSampled.rds")

featureCols = colnames(dSamp)[!(colnames(dSamp) %in% c("class", "time", "record_time", "id", "loc", "...1"))]
discreteCols <- featureCols[featureCols %in% c("ntile") ]
discreteCols <- c(discreteCols, featureCols[grep("pin_", featureCols)])
scaleCols <- featureCols[!featureCols %in% discreteCols]

# scale by location
dSamp[, (scaleCols) := lapply(.SD, function(x) as.vector(scale(x))), .SDcols = scaleCols, by=loc]

# param for corr
sampSize <- 10000
threshold <- 0.7

### for Pocket 
loc <- "pock"
sampleIndex <- sample(1:nrow(dSamp[loc==loc]), size = sampSize)
dCorrPock <- dSamp[loc==loc][sampleIndex, ..featureCols] %>% corr_simple( threshold, title = loc)

### Plot of histogram 
p <- dSamp[loc==loc][sampleIndex, ..featureCols] %>% 
  gather() %>%  
  ggplot(aes(value)) + 
  geom_histogram() +
  facet_wrap(~key, scales = 'free')
ggsave(paste("~/scratch/BART_test/image/descriptive/histogram", loc, ".png"))

### Biariate plot 
d <- setDF(dSamp[loc==loc][sampleIndex, ..featureCols]) 
for(i in 1:nrow(dCorrPock)){
  X <- d[, which(colnames(d) %in% dCorrPock[i, "Var1"])]
  Y <- d[, which(colnames(d) %in% dCorrPock[i, "Var2"])]
  png(paste("~/scratch/BART_test/image/descriptive/corrPlots/",loc,i,   dCorrPock$Var1[i], dCorrPock$Var2[i] ,".png"), height = 4, width =4, units = "in", res = 100)
  plot(x=X, y=Y, pch =".", main=paste(dCorrPock$Var1[i], "   :  ",  dCorrPock$Var2[i]))
  dev.off()
}

### for hand
loc <- "hand"
sampleIndex <- sample(1:nrow(dSamp[loc==loc]), size = sampSize)
dCorrHand <- dSamp[loc==loc][sampleIndex, ..featureCols] %>% corr_simple( threshold, title = loc)

### for back
loc <- "back"
sampleIndex <- sample(1:nrow(dSamp[loc==loc]), size = sampSize)
dCorrBack <- dSamp[loc==loc][sampleIndex, ..featureCols] %>% corr_simple( threshold, title = loc)


coorBMJFeatureDf <- rbind(
  data.frame(loc = "pock", dCorrPock),
  data.frame(loc = "hand", dCorrHand),
  data.frame(loc = "back", dCorrBack)
)
saveRDS(coorBMJFeatureDf, "~/scratch/BART_test/data/corrData_070.rds")


# no need to plot this, as all 4 features will be used for larning anyways
if(bool_plot){
  cat("Disaggregated 30 Hz data, sample of 50000 data points")
  cat("en: euclidian norm, en_d1: norm of first order difference, en_d2: norm of second_order difference, and angle 1 and 2")
  sampleIndex <- sample(1:nrow(dSamp[loc=="pock"]), size = 50000)
  dSamp[loc=="pock"][sampleIndex, c("x", "y", "z", "en", "en_d1", "en_d2", "ang1")] %>%  pairs(pch = ".")
  
  cat("Aaggregated 1 Hz data, sample of 50000 data points")
  sampleIndex <- sample(1:nrow(dSampFeatures[loc=="pock"]), size = 50000)
  dSampFeatures[loc=="pock"][sampleIndex, c("x", "y", "z", "en", "en_d1", "en_d2", "ang1", "ang2")] %>%  pairs(pch = ".")
}


