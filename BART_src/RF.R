######################################################################################################
#
# Codes to run random forest on the same data as BART are fit 
# Codes to run ML are borowed from https://github.com/khataei/Beap-Archive/blob/master/R/ApplyModels.R
# 
######################################################################################################

# There are only two inputs to this file: Training and Test data frame, each containing target and predictor vairalbes 

#trainData <- datTrainRF
#testData <- datTestRF


# Import libraries
library(data.table)
library(tidyverse)
library(lubridate)
library(activityCounts)
library(magrittr)
library(scales)
library(caret)
library(stringr)
library(MLmetrics)
library(PRROC)
library(gridExtra)
library(kableExtra)
library(imputeTS)
library(tictoc)
library(doParallel)
library(ranger)


#### Run Params -----------------------------
#scale_center = FALSE
#cv_folds = 0
#shrink = 1
#save_results_on_disk = TRUE
#return_plots = TRUE
#RF_mtry = 2
#min_node_size = 50
#cores = 10

funcRF = function(trainData, 
                  testData, 
                  scale_center=FALSE, 
                  cv_folds=0, 
                  shrink = 1, 
                  save_results_on_disk = TRUE,
                  return_plots = TRUE, 
                  RF_mtry = 2,
                  min_node_size = 50, 
                  cores = 1){
  
  # Parallel and time to see if caret parallel works
  tic("Preprocessing")
  
  if (cores != 1) {
    message("Parallel computing is activated!")
    cl <- makePSOCKcluster(cores)
    registerDoParallel(cl)
  }else {
    message("parallel computing is not activated. For parallel computing set it to an integer greater than 1")
  }
  
  
  # Create train and test to train and evalute the model
  seed <- 2020
  set.seed(seed)
  
  if (scale_center) {
    preProcValues <- preProcess(training_df, method = c("center", "scale"))
    training_df <- predict(preProcValues, training_df)
    testing_df <- predict(preProcValues, testing_df)
    message("Data is scaled and centered")
  }
  
  
  # To store the model and all the performance metric
  results <- NULL
  # To store plots
  plts <- NULL
  # To store confusion matrix
  cf_mat <- NULL
  # To return main results for each model
  accuracies <- NULL
  # to return RF feature importance
  importance <- NULL
  # The end of preprocessing step
  toc()
  
  
  ### Cross_validation ---------------
  if (cv_folds <= 0) {
    fitControl <-
      trainControl(method = "none", classProbs = TRUE)
    message("Cross-validation is not being used, set cv_folds to a positive number to use cross-validation")
  } else if (cv_folds > 0){
    cv_folds  %<>%  ceiling()
    fitControl <-
      trainControl(method = "cv", number = cv_folds, classProbs = TRUE)
    message(paste0(cv_folds, " fold cross-validation is being used"))
  }
  
      
  ### Training -----------------------
  # Ranger is a fast implementation of random forests (Breiman 2001)
  # The method is none becuase we have test and train data
  tic("RF took")
  message("Starting RF")
  
  model_name <- "ranger"
  train_control_method <- "none"
  model_mtry <- RF_mtry
  model_splitrule <- "extratrees"
  model_min_node_size <- min_node_size
  
  model_A <- train(
    y ~ .,
    # data = training_df,
    #data = trainData,
    data = trainData, 
    method = model_name,
    trControl = fitControl,
    verbose = FALSE,
    importance = "impurity",
    tuneGrid = data.frame(
      mtry = model_mtry,
      splitrule = model_splitrule,
      min.node.size = model_min_node_size
    ),
    metric = "ROC"
  )
  
  
  pred <- stats::predict(model_A, newdata = testData)
  
  # To calculate area AUC we need probabilies and predicted classes in a single dataframe
  pred_prob <- data.frame(obs =  testData$y, pred = pred)
  
  pred <- stats::predict(model_A, newdata = testData, type = "prob")
  pred_prob <- bind_cols(pred_prob, pred)
  
  # Calculate different metrics
  metrics <- multiClassSummary(data = pred_prob,
    lev = levels(testData$y)) %>%
    as.data.frame()
  
  # Return the metric in a nicer format
  metric_names <- rownames(metrics)
  metrics  %<>% data.table::transpose()
  colnames(metrics) <- metric_names
  rownames(metrics) <- "RF"
  accuracies  %<>%  rbind(metrics)
  
  
  # CF need a different format of prediction results so recalcuate
  pred <- stats::predict(model_A, newdata = testData)
  
  # Calculate confusion matrix
  cf_matrix <-
    confusionMatrix(
      data = pred,
      reference = testData$y,
      mode = "prec_recall"
    )
  
  # Create a list of the model and the results to save
  results[["RF"]] <-
    list(
      split_seed = seed,
      model_name = model_name,
      model = model_A,
      train_control_method = train_control_method,
      tune_parameters = c(model_mtry, model_splitrule, model_min_node_size),
      cf_matrix = cf_matrix
    )
  
  
  if (return_plots) {
    plts[["RF"]] <- cf_matrix$table %>%
      data.frame() %>%
      ggplot2::ggplot(aes(Prediction, Reference)) +
      geom_tile(aes(fill = Freq), colour = "gray50") +
      scale_fill_gradient(low = "white", high = muted("gray80")) +
      geom_text(aes(label = Freq)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle("RF")
  }
  importance <- caret::varImp(model_A,scale = FALSE)
  
  return(list(results = results, plots = plts, cf_mat = cf_mat, acc = accuracies, ranking = importance))

  toc()
}
