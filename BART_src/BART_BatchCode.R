#!/usr/bin/env Rscript
#--------------------------------------------------------------------#
#
#   BART codes for accelerator data, also running random forest 
# 
#   Predictive features are generated in accordance with: 
#     Khataeipour SJ,BMJ Open Sport & Exercise Medicine 2022;8:e001242.doi: 10.1136/bmjsem-2021-001242
#
#   This file can be run in batch a batch mode perfoming leave-one-out cross-validation for all participants using compute canada's SLURM
#      See scripts to call this file via SLURM,  "BART_job_test_argTest_looCV_array.sh"
#
#   Feature set are generated from raw accelerometer data calculated by a file: DataPrep_andDescription_5sec.R
#--------------------------------------------------------------------#


rm(list=ls())

# Libraries 
require(data.table)
require(R.utils)
require(BART)
require(ggplot2)
require(gridExtra)
require(dplyr)
require(tidyverse)
require(parallel)

set.seed(10000)

# Arguments for this script, read when this script is run as batch mode 
args=(commandArgs(TRUE))

# Run interactively or as a batch? 
# If the vector of arguments exists, run as batch mode in SLURM. Else run as interactively 
manualRun = ifelse(length(args)==0,  TRUE, FALSE) 

if(manualRun){cat("-------------------- Interactive ------------------")
  }else{
    cat("------------------------ Batch ---------------------")
    # Display arguments if batch mode 
    paste("Batch or not", as.integer(args[1]))
    paste("nInstance",  as.integer(args[2]))
    paste("ndPost", as.integer(args[3]))
    paste("keepEvery", as.integer(args[4]))
    paste("mcCore", as.integer(args[5]))
    paste("nSkip", as.integer(args[6]))
    paste("deviceLoc", as.character(args[7]))
    paste("dataSource", as.integer(args[8]))
    paste("Note", args[9])
    paste("personIDfile", as.character(args[10]))
    paste("bool_scale", as.integer(args[11]))
    paste("bool_dropCoorVar", as.integer(args[12]))
    paste("source", as.character(args[13]))
    paste("nTree",  as.integer(args[14]))
    paste("bool_loocv", as.integer(args[15]))
}





# Set params ------------------------------------------------------------------------------------------------
# Save training data for BART or not 
bool_saveData <- TRUE

# Save training data or not for random forest 
bool_saveTrainData <- FALSE

# Run BART and random forest (RT) model or not
bool_runBART<-TRUE
bool_runRF <- TRUE

# Apply BART object to test data? 
bool_BARTeval <- TRUE 

# If sample == TRUE, perform split sample evaluation (not recommended) rather than leave-one-out-cross-validation  (recommended)
bool_sample = FALSE 
bartType = 1 # mbart, or mbart2 - see BART manual 

# Save traceplots or not (should be set to TRUE, if running in batch model)
bool_tracePlots = FALSE

# Function to calculate the correlation of features 
source("~/scratch/BART_test/functions/corrPlot.R")


# Read arguments -------------------------------------------------------------------------
# Arguments, either provided in this script or externally defined when in batch mode 
if(manualRun){
    nInstance <- 10000      # NOT used currently. Number of data points to sample, if split sample validation is done (not recommended)
    ndPost <- 4000          # Number o sample - default  4000 
    keepEvery <- 100        # Thinning - default 100
    mcCore <- 8             # Number of cores, if using Unix like machine. 
    nSkip <- 3000           # Burn-in to exclude from inference 
    personIDfile <- "~/scratch/BART_test/BART_src/personIDVec.R" # Default File to pick person 
    deviceLoc <- "pock"     # Location of device, iether pocket, arm, or backpack
    dataSource <- 5         # 5 - indicating data with 5 sec window data 
    Note <- ""              # Any note to be attached to saved model object
    bool_scale <- 1         # Scale predictive feature (mean center and 1 SD scaling) or not 
    bool_dropCoorVar <- 1   # Drop strongly correlated feature or not? set 1/0
    nTree <- 50             # Number of BART trees, default 50 
    bool_loocv <- 0         # Perform leave-one-out-cross-validation if set to 1.
    xMonitorDiag =  c(1,50, 100, 150, 200, 400, 600)  # Data points to monitor for trace plots. Each person has about 650 time points 
  }else{ # Below same arguments as above, but provided by batch shell script that calls this file 
    nInstance <- as.integer(args[2])      
    ndPost <- as.integer(args[3])      
    keepEvery <- as.integer(args[4])      
    mcCore <- as.integer(args[5]) 
    nSkip <- as.integer(args[6])
    deviceLoc <- as.character(args[7])
    dataSource <- as.integer(args[8])
    Note <- args[9]
    personIDfile <- as.character(args[10])  # This is a person ID (range from 1-37). 
    bool_scale <- as.integer(args[11])
    bool_dropCoorVar <- as.integer(args[12])
    source(as.character(args[13]))
    nTree <- as.integer(args[14])
    bool_loocv <- as.integer(args[15])
}


# Print some args for logging if running batch 
cat("CHECK ARGS----------------------------------------------------\n")
cat(c("LOOCV", as.integer(args[15])))
cat(c("PERSON INDEX ARG ", as.integer(args[10])))
cat("-----------------------------------------------------------\n")



# Load data ----------------------------------------------------------------------------------
# Path where accelerator data features are saved 
datPath <- "~/scratch/BART_test/data/"
datName <- switch(dataSource, 
                  "bmjFeatures_reSampled.rds", 
                  "bmjFeatures_reSampledAgg.rds", 
                  "sensorFeaturesResampledAgg.rds",
                  "sensorFeaturesResampled.rds", 
                  "bmjFeatures_reSampled_5sec.rds", # Use only this option - other data are too large or used inappropriate resampling methods  
                  "bmjFeatures_reSampledAgg_5sec.rds"
)
datOriginal <- readRDS(paste(datPath, datName, sep = ""))

# Fix variable name 
datOriginal$vec_mag_g <- datOriginal$`vec_mag-g`
datOriginal$`vec_mag-g` <- NULL

# Subset Data by location pocket, hand, or arm 
datOriginal <- datOriginal[loc == deviceLoc]


# Split data by person ID - number of ID in test data is supplied instead of file name  -----
# Load file containing ID of people to include in train/test data 
cat(c("PERSON ID ----------------------------------------------------------\n"))
if(bool_loocv == 0){ # no leave-one-out-CV
  if(is.na(as.numeric(personIDfile))){
      source(personIDfile) #Load person ID(s) from this file, as training data - IDs not listed in this file will be used for testing 
      selectPerson <- personVec
      print(c("PersonID File", personIDfile))
      print(c("PersonIDs for training", selectPerson))
    }else{ # if the variable person ID file is numeric (rather than file name) with a single number, it then indicate the number of people to be sampled as training data, randomly  
      numTrainPerson <- as.numeric(personIDfile)
      selectPerson <- sample(unique(datOriginal$id), size = numTrainPerson)
      selectPerson <- sort(selectPerson)
      print(c("No LOOCV, PersonIDs automatically selected for training :    ", selectPerson))
    }
  }else if (bool_loocv ==1 & length(personIDfile)==1){ # Loaded file contain one value indicating person index (1 ... 45) to be used as TEST data, the rest of participants not listed indicate training. This index can be supplied as scalar or array in SLURM  
    if(as.numeric(personIDfile) > 40) stop("Hiroshi: Person Index needed ")
    print(c("LOOCV, TestPersonIndex:", personIDfile))
    print(c("PersonID:",  unique(datOriginal$id)[as.numeric(personIDfile)]))
    selectPerson <- unique(datOriginal$id)[-as.numeric(personIDfile)] # Assign everyone else as training 
    print(c("LOOCV, training ID:", selectPerson))
}
cat(c("PERSON ID ----------------------------------------------------------\n"))



### Print run parameters ----------------------------------------------------------------
print(c("nInstance", nInstance))
print(c("ndPost", ndPost))
print(c("mcCore", mcCore))
print(c("KeepEvery", keepEvery))
print(c("Note", Note))
print(c("SKIP", nSkip))
print(c("shortNote", Note))
print(c("selected person", selectPerson))
print(c("Device lcoation", deviceLoc))
print(c("dataSource", dataSource))
print(c("Scaling", bool_scale))
print(c("Drop Correlated features ? ", bool_dropCoorVar))
cat(c("----------------------------------------------------------"))


# Stop the script if the number of argument is incorrect, if runing in batch mode 
if(length(args)<15 & manualRun==F){
  stop("---15 arguments must be supplied ('name' (text) and 'numer' (integer) )", call.=FALSE)
}






# Scaling of variables to full data set  ----------------------------------------------
# Get name of features 
featureCols = colnames(datOriginal)[!(colnames(datOriginal) %in% c("class", "time", "record_time", "id", "loc", "...1", "windowID"))]
# Scale discrete variables? 
discreteCols <- featureCols[featureCols %in% c("ntile") ]
discreteCols <- c(discreteCols, featureCols[grep("pin_", featureCols)])
scaleCols <- featureCols[!featureCols %in% discreteCols]

# Data ID 1,2,5,6 are data with feature sets provided by Khataeipour, BMJ Sports Med (2022) paper listed above
if(dataSource %in% c(1,2, 5,6)){
  if(bool_scale){
    datOriginal[, (scaleCols) := lapply(.SD, function(x) as.vector(scale(x))), .SDcols = scaleCols, by = loc]
    #### Check standardization was done or not - histograms
    ##pD <- datOriginal[, ..scaleCols] %>% gather() 
    ##p <- ggplot(gather(pD), aes(value)) + 
    ##  geom_histogram(bins = 10) + 
    ##  facet_wrap(~key, scales = 'free_x')
    ##rm(pD)
    #datOriginal %>% summarise(across(where(is.numeric), mean))
    #datOriginal %>% summarise(across(where(is.numeric), sd))
  }
  cat("Scaling of features \n ", bool_scale)
}


# Give a unique index to original data
datOriginal[, indexOriginalData := 1:.N]


# Save all participant data, if performing model evaluation on test data in a separate script 
if(bool_saveData){
  fileNameData <- paste("BART_inputTrainTestData_allPerson","_",
                        "_dataSource_", dataSource, 
                        "_loc_", deviceLoc,
                        "_scale_", bool_scale, 
                        sep = "")
  
  datOriginal$...1 <- NULL
  saveRDS(datOriginal, paste("~/scratch/BART_test/data/", 
                             fileNameData, 
                             ".rds", 
                             sep = ""))
}




# Drop Correlated predictors  ------------------------------
cat("Drop correlated features? \n")
cat((dataSource %in% c(1,2,5,6) & bool_dropCoorVar==1))

# list correlation of variables and correlation matrix
if(dataSource %in% c(1,2, 5,6 ) & bool_dropCoorVar==1){
  corrSampleSize = 20000 #cannot exceed the sample size
  corrThreshold = 0.8 # threshold corrlation, to remove variables 
  sampleIndex <- sample(nrow(datOriginal), size = corrSampleSize)
  corDf <- datOriginal[sampleIndex, ..featureCols] %>% corr_simple(corrThreshold, title = "")
  
  datOriginal <- datOriginal %>% select(-any_of(corDf$Var2))
} 



# Remove a few records of 7 MTEs at the beginning of some participants' time-series. 
# Since, strangely, 4 people start with Running 7 METS in their activity label for first a few seconds - they need to be removed. 
datOriginal <- datOriginal %>% 
  group_by(id) %>% 
  filter(!(windowID  < 50 & class == "Running 7 METs")) %>% 
  ungroup()


# Create data table, before subsetting data into train and test 
setDT(datOriginal)









# Extract training data ------------------------------------
# Subset data by person, after scaling so that unsampled data (test data) will also be in scale as well 
dSubset <- datOriginal[id %in% c(selectPerson)]
print(c("Person ID------------", unique(dSubset$id)))
print(c("Device Location------------", unique(dSubset$loc)))


# randomly select data points, if variable sampling is set - not recommended but one should rather run leave-one-out CV
if(bool_sample){
  sampleID <- sample(1:nrow(dSubset), nInstance)
  dSubset <- dSubset[sampleID]
  print("Random sampling across person was performed")
}


# Class for training data 
y.train <- factor(dSubset$class, levels=c("Lying", "Sitting", "Self Pace walk", "Running 3 METs", "Running 5 METs", "Running 7 METs"))
y.train = as.numeric(y.train)


# Prepare training data made up of the  predictor variables
x.train = as.data.frame(dSubset[, intersect(names(dSubset), featureCols), with=FALSE])



# Timepoints to monitor for the mixing of MCMC 
dMonitor <- dSubset[xMonitorDiag]

# Person and time ID in training data, to be matched to the BART output later 
trainDataID <- dSubset %>% select(id, loc, record_time, windowID, indexOriginalData)

# FileName CV - if LOOCV is 0, it is left as zero. If LOOCV is 1, then Person ID should be returned \
if(bool_loocv==0){
  indicator_looCV <- 0
}else{
  indicator_looCV <- paste(unique(datOriginal$id)[as.numeric(personIDfile)], "_ArrayIndex_", personIDfile, sep = "")
}


print('Outcome frequency in training data ---------------------------------------------------------------------------')
table(y.train)
print('Exposure  ---------------------------------------------------------------------------')
names(x.train)
print('Data point to monitor')
print(xMonitorDiag)



### RUN BART model -------------------------------------------------------------------------
# Change the platform (e.g., win) if necessary 
if(.Platform$OS.type == "unix" & bool_runBART){
  start_time <- Sys.time()
  post=mc.mbart(x.train, 
                y.train,
                x.train[xMonitorDiag, ], 
                ndpost = ndPost, 
                keepevery = keepEvery, 
                mc.cores = mcCore, 
                printevery = 100, 
                nskip = nSkip, 
                ntree = nTree) 
  
  end_time <- Sys.time()
  mcmcTime <- end_time - start_time
  saveTime <- format(Sys.time(), "____%b_%d_%H:%M")
  
  # Save BART run object ---------------------------------
  listOut <- list(post=post, 
    size = nInstance, 
    mcmcTime=mcmcTime, 
    keepEvery=keepEvery, 
    ndPost =ndPost, 
    bartType = bartType, 
    indexMonitorDiag = xMonitorDiag, 
    #indexTrain = dat$index, 
    trainData = x.train,
    trainDataID = trainDataID, 
    trainMonitorLabel = y.train[xMonitorDiag], 
    trainY = y.train, 
    personID = selectPerson, 
    loc = deviceLoc, 
    dataSource = dataSource, 
    scaling = bool_scale, 
    corrFeatures = bool_dropCoorVar, 
    loocv = indicator_looCV
    )

  fileName <- paste("BART","_",
    "_ndPost_", ndPost, 
    "_Thin_", keepEvery, 
    "_Core_", mcCore,
    "_nSkip_", nSkip, 
    "_bartType_", bartType,
    "_note_", Note, 
    "_saveTime_", saveTime, 
    "_dataSource_", dataSource, 
    "_loc_", deviceLoc, 
    "_numTrain_", length(selectPerson), 
    "_Scale_", bool_scale, 
    "_dropCoorFeatures_", bool_dropCoorVar, 
    "_nTree_", nTree, 
    "loocv_", indicator_looCV,
    #"____MANUALRUN____", 
    sep = ""
    )

  saveRDS(listOut, paste("~/scratch/BART_test/data_output/", 
    fileName, 
    ".rds", 
    sep = "")
    )
  
  # Plot MCMC traces -------------------------------
  source("~/scratch/BART_test/functions/funcPlot.R")
    
  
  # Save plots into list 
    l <- list()
    l[[1]] <- funcGetSamplePlot(listOut, title = paste("LOC:", deviceLoc, "  :  ", "PersonID:", listOut$personID), 0)
    for(i in 2:length(xMonitorDiag)){
      l[[i]] <- funcGetSamplePlot(listOut, title = paste(xMonitorDiag[i]), i-1)
   } 
    
    p <- arrangeGrob(grobs = l, ncol = 1)
    
    if(bool_tracePlots){
      ggsave(paste("~/scratch/BART_test/image/mcmc/", fileName ,".png", sep = ""), p, width = 11, height = length(xMonitorDiag)*2)
    }
    
    print("TIME--------------------------------------------------------------------------------------------------")
    cat(mcmcTime)
    print("TIME--------------------------------------------------------------------------------------------------")
    
  }# End of runBART

  


# Save training data separately? 
if(bool_saveTrainData){
    fileName <- paste("BART_trainData","_",
                      "_dataSource_", dataSource, 
                      "_loc_", deviceLoc,
                      "_scale_", bool_scale, 
                      sep = "")
    
    datOriginal$...1 <- NULL
    saveRDS(x.train, paste("~/scratch/BART_test/data/", 
                        fileName, 
                        ".rds", 
                        sep = ""))
} # Save training data 









# Apply BART object to test data and select class with highest probability for each iteration ------------------------
if(bool_BARTeval){
  
  ### Functions ---- these are same as the lengthy codes below -----
  source("~/scratch/BART_test/functions/funcBARTPosterior.R")
  #names(post$treedraws$cutpoints) selected variables 
  trainPostSim<- funcSummrizePost(bartObj = listOut$post, 
                                  xData = listOut$trainData  %>%  select(names(post$treedraws$cutpoints)), 
                                  yData = listOut$trainY, 
                                  dataID = listOut$trainDataID, 
                                  cores = 8
                                  ) 
  cat(c("Posterior mean of accuracy: Training data:", sum(trainPostSim$matchMCMC)/nrow(trainPostSim)))
  
  # Posterior of out-of-sample test data 
  # Prepare data set for testing
  datOriginalTest <- datOriginal %>% filter(!id %in% listOut$personID) # not in the training 
  cat("Number of people in test data", length(unique(datOriginalTest$id)))
  
  # Create test data 
  # features in test data
  x.test <- datOriginalTest %>% 
    select(colnames(listOut$trainData))
  setDF(x.test)
  
  # Outcome in test data 
  y.test <- as.integer(factor(datOriginalTest$class, 
  levels=c("Lying", "Sitting", "Self Pace walk", "Running 3 METs", "Running 5 METs", "Running 7 METs")))
  
  # Person and time ID in test data
  dataID <- datOriginalTest %>%  
    select(colnames(listOut$trainDataID))

  rm(datOriginalTest)
  
  # Apply BART object to test data 
  test_PostSim<- funcSummrizePost(bartObj = listOut$post, 
                                  xData = x.test, 
                                  yData = y.test, 
                                  dataID = dataID, 
                                  cores = 8, #0.9347499
                                  #trainingPersonIDVec <- listOut$personID
                                  )
  cat(c("Posterior mean accuracy, Test data:", sum(test_PostSim$matchMCMC)/nrow(test_PostSim)))
  
  # Combine posterior distribution from train and test data 
  listOutPost <- list()
  listOutPost[[1]] <- trainPostSim
  listOutPost[[2]] <- test_PostSim
  listOutPost[[3]] <- c("loocv_", indicator_looCV)
  saveTime <- format(Sys.time(), "__%d_%H:%M")
  
  # File name 
  fileNamePost <- paste("BART_PostSimTestTrain_",
                    "_ndPost_", ndPost, 
                    "_Thin_", keepEvery, 
                    "_Core_", mcCore,
                    "_xTrain_","_numTrace_4",
                    "_nSkip_", nSkip, 
                    "_bartType_", bartType,
                    "_note_", Note, 
                    "_saveTime_", saveTime, 
                    "_dataSource_", dataSource, 
                    "_loc_", deviceLoc, 
                    "_ID_", length(selectPerson), 
                    "_Scale_", bool_scale, 
                    "_dropCoorFeatures_", bool_dropCoorVar, 
                    "_nTree_", nTree, 
                    "loocv_", indicator_looCV,
                    sep = ""
  )
  
  saveRDS(listOutPost, paste("~/scratch/BART_test/data_output/", fileNamePost, ".rds", sep = ""))
} # End of BART evaluation on test data 









### Random Forest (RF)------------------------
if(bool_runRF){
  source("~/scratch/BART_test/BART_src/RF.R")
  
  # Data frame for RF
  datTrainRF <- data.frame(x.train, y = y.train)
  datTrainRF$y <- factor(datTrainRF$y)
  levels(datTrainRF$y) <- c("L", "S", "W", "R3", "R5", "R7")
  
  
  # Create test data - Exclude data with training set
  dTest <- datOriginal %>% filter(!id %in% selectPerson) 
  
  # Select same 
  datTestRF <- dTest %>% select(colnames(x.train)) %>% data.table
  
  datTestRF$y<- factor(dTest$class, levels=c("Lying", "Sitting", "Self Pace walk", "Running 3 METs", "Running 5 METs", "Running 7 METs"))
  levels(datTestRF$y) <- c("L", "S", "W", "R3", "R5", "R7")
  
  print(c("Length and dimention of training data", dim(datTrainRF)))
  print(c("Length and dimention of test data", dim(datTestRF)))
  
  # Compare mean of features between the test and training data 
  cbind(
    datTrainRF  %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% t,
    datTestRF  %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% t
  )
  
  # Train RF
  resultTrainRF <- funcRF(trainData = datTrainRF, testData = datTrainRF, cv_folds = 5)
  
  #Accuracy of RF on training data 
  resultTrainRF$results$RF$cf_matrix
  resultTrainRF$plots
  cat(c("Overall TRAIN Performance", resultTrainRF$results$RF$cf_matrix$overall))
  
  # Get test results 
  resultTestRF <- funcRF(trainData = datTrainRF, testData = datTestRF, cv_folds = 5)
  
  # Accuracy of RF and test data 
  resultTestRF$results$RF$cf_matrix
  resultTestRF$plots
  cat("RF test data results --------------------------------------------------------------------")
  cat(c("Overall TEST Performance", resultTestRF$results$RF$cf_matrix$overall))
  

  # Save   
   if(!is.null(listOut) & !is.null(fileName)){ # Save in BART object above, if BART was run   
      listOut$resultTrainRF <- resultTrainRF 
      listOut$resultTestRF <- resultTestRF 
      saveRDS(listOut, paste("~/scratch/BART_test/data_output/", 
                             fileName, 
                             ".rds", 
                             sep = "")
      )
   }else{ # Create new object to save RF results, if BART was not run above 
     listOutRF <- list()
     listOutRF$resultTrainRF <- resultTrainRF 
     listOutRF$resultTestRF <- resultTestRF 
     rfileName <- paste("RF","_",
                     "_dataSource_", dataSource, 
                     "_loc_", deviceLoc, 
                     "_saveTime_", saveTime, 
                     "_ID_", length(selectPerson), 
                     "_Scale_", bool_scale, 
                     "_dropCoorFeatures_", bool_dropCoorVar, 
                     "loocv_", indicator_looCV,
                    # "__MANUALRUN__",
                     sep = "")
   saveRDS(listOutRF, paste("~/scratch/BART_test/data_output/", 
                          rfileName, 
                          ".rds", 
                          sep = "")
           )
  }
} # End of Random Forest train and test 



