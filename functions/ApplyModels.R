# ML on physical acitviity, from this. MOdified to incorporated BART train and test data 
# https://github.com/khataei/Beap-Archive/blob/master/R/ApplyModels.R

ApplyModels_Hiroshi <-
  function(#working_df,
    training_df,
    testing_df,
           model_names = c("RF"),
           #model_names = c("RF", "LDA", "NB", "SVM", "KNN" , "DT", "XGB"),
    bool_split = FALSE, 
           split_ratio = 0.66,
           scale_center = FALSE,
           cv_folds = 0,
           shrink = 1,
           save_results_on_disk = TRUE,
           return_plots = TRUE,
           RF_mtry = 2,
           min_node_size = 50,
           cores = 1){
    # # Parallel and time to see if caret parallel works
    tic("Preprocessing")
    
    if (cores != 1) {
      message("Parallel computing is activated!")
      cl <- makePSOCKcluster(cores)
      registerDoParallel(cl)
    }
    else {
      message("parallel computing is not activated. For parallel computing set it to an integer greater than 1")
    }
    
    
    # Create train and test to train and evalute the model
    seed <- 2020
    set.seed(seed)
    
    #modified by Hiroshi, as train and test data are already created in BART run codes 
    if(bool_split){
      working_df %<>% dplyr::sample_frac(shrink)
      message(paste0(shrink * 100, " % of the data will be used"))
      
      training_indices <-createDataPartition(working_df$trimmed_activity, p = split_ratio, list = FALSE)
      
      training_df <- working_df %>% dplyr::slice(training_indices)
      testing_df <- working_df %>% dplyr::slice(-training_indices)
    }
    
    if(scale_center) {
      preProcValues <- preProcess(training_df, method = c("center", "scale"))
      training_df <- predict(preProcValues, training_df)
      testing_df <- predict(preProcValues, testing_df)
      message("Data is scaled and centered")
    }
    
    message(paste0("Data is devided into training and test set "))
    message(paste0(
      "Training set has ",
      nrow(training_df),
      " rows and ",
      ncol(training_df),
      " columns"
    ))
    message(paste0(
      "Testing set has ",
      nrow(testing_df),
      " rows and ",
      ncol(testing_df),
      " columns"
    ))
    
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
    
    
    if (cv_folds <= 0){
      fitControl <- trainControl(method = "none", classProbs = TRUE)
      message("Cross-validation is not being used, set cv_folds to a positive number to use cross-validation")
    }else if (cv_folds > 0){
      cv_folds  %<>%  ceiling()
      fitControl <-
        trainControl(method = "cv", number = cv_folds, classProbs = TRUE)
      message(paste0(cv_folds, " fold cross-validation is being used"))
    }
    
    
    # ------------------------------------- LDA --------------------------------------
    if ("LDA" %in% model_names) {
      tic("LDA took")
      message("Starting LDA")
      model_name <- "lda"
      train_control_method <- "none"
      model_parameter <- 10
      
      model_A <- train(
        trimmed_activity ~ .,
        data = training_df,
        method = model_name,
        trControl = fitControl,
        verbose = FALSE,
        tuneGrid = data.frame(parameter = model_parameter),
        metric = "ROC"
      )
      
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # To calculate area AUC we need probabilies and predicted classes in a single dataframe
      pred_prob <-
        data.frame(obs =  testing_df$trimmed_activity,
                   pred = pred)
      pred <-
        stats::predict(model_A, newdata = testing_df, type = "prob")
      pred_prob <- bind_cols(pred_prob, pred)
      
      # Calculate different metrics
      metrics <-
        multiClassSummary(data = pred_prob,
                          lev = levels(testing_df$trimmed_activity)) %>%
        as.data.frame()
      # Return the metric in a nicer format
      metric_names <- rownames(metrics)
      metrics  %<>% data.table::transpose()
      colnames(metrics) <- metric_names
      rownames(metrics) <- "LDA"
      accuracies  %<>%  rbind(metrics)
      
      # CF need a different format of prediction results so recalcuate
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # Calculate confusion matrix
      cf_matrix <-
        confusionMatrix(
          data = pred,
          reference = testing_df$trimmed_activity,
          mode = "prec_recall"
        )
      
      
      # create a list of the model and the results to save
      results[["LDA"]] <-
        list(
          split_seed = seed,
          model_name = model_name,
          model = model_A,
          train_control_method = train_control_method,
          tune_parameters = c(model_parameter),
          cf_matrix = cf_matrix
        )
      
      if (return_plots) {
        plts[["LDA"]] <-  cf_matrix$table %>%
          data.frame() %>%
          ggplot2::ggplot(aes(Prediction, Reference)) +
          geom_tile(aes(fill = Freq), colour = "gray50") +
          scale_fill_gradient(low = "beige", high = muted("chocolate")) +
          geom_text(aes(label = Freq)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("LDA")
      }
      toc()
    }
    
    
    
    
    #------------------------------ Random Forest _ Ranger package----------------------
    if ("RF" %in% model_names) {
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
        trimmed_activity ~ .,
        data = training_df,
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
      
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # To calculate area AUC we need probabilies and predicted classes in a single dataframe
      pred_prob <-
        data.frame(obs =  testing_df$trimmed_activity,
                   pred = pred)
      pred <-
        stats::predict(model_A, newdata = testing_df, type = "prob")
      pred_prob <- bind_cols(pred_prob, pred)
      
      # Calculate different metrics
      metrics <-
        multiClassSummary(data = pred_prob,
                          lev = levels(testing_df$trimmed_activity)) %>%
        as.data.frame()
      # Return the metric in a nicer format
      metric_names <- rownames(metrics)
      metrics  %<>% data.table::transpose()
      colnames(metrics) <- metric_names
      rownames(metrics) <- "RF"
      accuracies  %<>%  rbind(metrics)
      
      
      # CF need a different format of prediction results so recalcuate
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # Calculate confusion matrix
      cf_matrix <-
        confusionMatrix(
          data = pred,
          reference = testing_df$trimmed_activity,
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
          scale_fill_gradient(low = "beige", high = muted("chocolate")) +
          geom_text(aes(label = Freq)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("RF")
      }
      importance <- caret::varImp(model_A,scale = FALSE)
      toc()
    }
    
    #------------------------------ Boosting Trees----------------------
    if ("XGB" %in% model_names) {
      tic("XGB took")
      message("XGB")
      model_name <- "xgbTree"
      train_control_method <- "none"
      model_max_depth <- 4
      model_nrounds<- 500
      model_eta <- 0.4
      model_min_child_weight <- 1
      model_subsample <- 1
      model_gamma <- 1
      model_colsample_bytree <- 1
      
      
      model_A <- train(
        trimmed_activity ~ .,
        data = training_df,
        method = model_name,
        trControl = fitControl,
        verbose = FALSE,
        importance = "impurity",
        tuneGrid = data.frame(
          max_depth = model_max_depth,
          nrounds = model_nrounds,
          eta = model_eta,
          min_child_weight = model_min_child_weight,
          subsample = model_subsample,
          gamma = model_gamma,
          colsample_bytree = model_colsample_bytree
          
        ),
        metric = "ROC"
      )
      
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # To calculate area AUC we need probabilies and predicted classes in a single dataframe
      pred_prob <-
        data.frame(obs =  testing_df$trimmed_activity,
                   pred = pred)
      pred <-
        stats::predict(model_A, newdata = testing_df, type = "prob")
      pred_prob <- bind_cols(pred_prob, pred)
      
      # Calculate different metrics
      metrics <-
        multiClassSummary(data = pred_prob,
                          lev = levels(testing_df$trimmed_activity)) %>%
        as.data.frame()
      # Return the metric in a nicer format
      metric_names <- rownames(metrics)
      metrics  %<>% data.table::transpose()
      colnames(metrics) <- metric_names
      rownames(metrics) <- "XGB"
      accuracies  %<>%  rbind(metrics)
      
      
      # CF need a different format of prediction results so recalcuate
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # Calculate confusion matrix
      cf_matrix <-
        confusionMatrix(
          data = pred,
          reference = testing_df$trimmed_activity,
          mode = "prec_recall"
        )
      
      
      # Create a list of the model and the results to save
      results[["XGB"]] <-
        list(
          split_seed = seed,
          model_name = model_name,
          model = model_A,
          train_control_method = train_control_method,
          tune_parameters = c(model_max_depth, model_nrounds, model_eta,
                              model_min_child_weight,model_subsample,
                              model_gamma,model_colsample_bytree),
          cf_matrix = cf_matrix
        )
      
      
      if (return_plots) {
        plts[["XGB"]] <- cf_matrix$table %>%
          data.frame() %>%
          ggplot2::ggplot(aes(Prediction, Reference)) +
          geom_tile(aes(fill = Freq), colour = "gray50") +
          scale_fill_gradient(low = "beige", high = muted("chocolate")) +
          geom_text(aes(label = Freq)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("XGB")
      }
      importance <- caret::varImp(model_A,scale = FALSE)
      toc()
    }
    
    
    
    
    #---------------------------------- Naive Bayes Classifier ----------------------------
    if ("NB" %in% model_names) {
      # The method is none becuase we have test and train data
      tic("NB took:")
      message("Starting NB")
      model_name <- "nb"
      train_control_method <- "none"
      model_fL <- 3
      model_usekernel <- TRUE
      model_adjust <- 1.5
      
      model_A <- train(
        trimmed_activity ~ .,
        data = training_df,
        method = model_name,
        trControl = fitControl,
        verbose = FALSE,
        tuneGrid = data.frame(
          fL = model_fL ,
          usekernel = model_usekernel ,
          adjust = model_adjust
        ),
        metric = "ROC"
      )
      
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # To calculate area AUC we need probabilies and predicted classes in a single dataframe
      pred_prob <-
        data.frame(obs =  testing_df$trimmed_activity,
                   pred = pred)
      pred <-
        stats::predict(model_A, newdata = testing_df, type = "prob")
      pred_prob <- bind_cols(pred_prob, pred)
      
      # Calculate different metrics
      metrics <-
        multiClassSummary(data = pred_prob,
                          lev = levels(testing_df$trimmed_activity)) %>%
        as.data.frame()
      # Return the metric in a nicer format
      metric_names <- rownames(metrics)
      metrics  %<>% data.table::transpose()
      colnames(metrics) <- metric_names
      rownames(metrics) <- "NB"
      accuracies  %<>%  rbind(metrics)
      
      # CF need a different format of prediction results so recalcuate
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # Calculate confusion matrix
      cf_matrix <-
        confusionMatrix(
          data = pred,
          reference = testing_df$trimmed_activity,
          mode = "prec_recall"
        )
      
      
      # Create a list of the model and the results to save
      results[["NB"]] <-
        list(
          split_seed = seed,
          model_name = model_name,
          model = model_A,
          train_control_method = train_control_method,
          tune_parameters = c(model_fL, model_usekernel, model_adjust),
          cf_matrix = cf_matrix
        )
      if (return_plots) {
        plts[["NB"]] <-  cf_matrix$table %>%
          data.frame() %>%
          ggplot2::ggplot(aes(Prediction, Reference)) +
          geom_tile(aes(fill = Freq), colour = "gray50") +
          scale_fill_gradient(low = "beige", high = muted("chocolate")) +
          geom_text(aes(label = Freq)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("NB")
      }
      toc()
    }
    
    
    #----------------------------------- k-Nearest Neighbors -----------------------------------
    if ("KNN" %in% model_names) {
      # using kknn package
      # The method is none becuase we have test and train data
      message("Starting KNN")
      tic("KNN took")
      
      
      model_name <- "kknn"
      train_control_method <- "none"
      model_kmax <- 3
      model_kernel <- "optimal" # Normal unweighted KNN
      model_distance <-
        1 # 1 for Manhatan , 2 for Euclidean distance
      
      
      model_A <- train(
        trimmed_activity ~ .,
        data = training_df,
        method = model_name,
        trControl = fitControl,
        verbose = FALSE,
        tuneGrid = data.frame(
          kmax = model_kmax,
          kernel = model_kernel,
          distance = model_distance
        ),
        metric = "ROC"
      )
      
      
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # To calculate area AUC we need probabilies and predicted classes in a single dataframe
      pred_prob <-
        data.frame(obs =  testing_df$trimmed_activity,
                   pred = pred)
      pred <-
        stats::predict(model_A, newdata = testing_df, type = "prob")
      pred_prob <- bind_cols(pred_prob, pred)
      
      # Calculate different metrics
      metrics <-
        multiClassSummary(data = pred_prob,
                          lev = levels(testing_df$trimmed_activity)) %>%
        as.data.frame()
      # Return the metric in a nicer format
      metric_names <- rownames(metrics)
      metrics  %<>% data.table::transpose()
      colnames(metrics) <- metric_names
      rownames(metrics) <- "KNN"
      accuracies  %<>%  rbind(metrics)
      
      # CF need a different format of prediction results so recalcuate
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # Calculate confusion matrix
      cf_matrix <-
        confusionMatrix(
          data = pred,
          reference = testing_df$trimmed_activity,
          mode = "prec_recall"
        )
      
      
      # Create a list of the model and the results to save
      results[["KNN"]] <-
        list(
          split_seed = seed,
          model_name = model_name,
          model = model_A,
          train_control_method = train_control_method,
          tune_parameters = c(model_kmax, model_kernel, model_distance),
          cf_matrix = cf_matrix
        )
      
      if (return_plots) {
        plts[["KNN"]] <- cf_matrix$table %>%
          data.frame() %>%
          ggplot2::ggplot(aes(Prediction, Reference)) +
          geom_tile(aes(fill = Freq), colour = "gray50") +
          scale_fill_gradient(low = "beige", high = muted("chocolate")) +
          geom_text(aes(label = Freq)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("KNN")
      }
      toc()
    }
    
    #------------------------------ Support Vector Machines with Polynomial Kernel ----------------------
    if ("SVM" %in% model_names) {
      # using  kernlab package
      # The method is none becuase we have test and train data
      tic("SVM took")
      message("Starting SVM")
      
      model_name <- "svmPoly"
      train_control_method <- "none"
      model_degree <- 3
      model_scale <- 1
      model_C <- 0.01
      
      model_A <- train(
        trimmed_activity ~ .,
        data = training_df,
        method = model_name,
        trControl = fitControl,
        verbose = FALSE,
        tuneGrid = data.frame(
          degree = model_degree,
          scale = model_scale,
          C = model_C
        ),
        metric = "ROC"
      )
      
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # To calculate area AUC we need probabilies and predicted classes in a single dataframe
      pred_prob <-
        data.frame(obs =  testing_df$trimmed_activity,
                   pred = pred)
      pred <-
        stats::predict(model_A, newdata = testing_df, type = "prob")
      pred_prob <- bind_cols(pred_prob, pred)
      
      # Calculate different metrics
      metrics <-
        multiClassSummary(data = pred_prob,
                          lev = levels(testing_df$trimmed_activity)) %>%
        as.data.frame()
      # Return the metric in a nicer format
      metric_names <- rownames(metrics)
      metrics  %<>% data.table::transpose()
      colnames(metrics) <- metric_names
      rownames(metrics) <- "SVM"
      accuracies  %<>%  rbind(metrics)
      
      # CF need a different format of prediction results so recalcuate
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # Calculate confusion matrix
      cf_matrix <-
        confusionMatrix(
          data = pred,
          reference = testing_df$trimmed_activity,
          mode = "prec_recall"
        )
      
      
      # Create a list of the model and the results to save
      results[["SVM"]] <-
        list(
          split_seed = seed,
          model_name = model_name,
          model = model_A,
          train_control_method = train_control_method,
          tune_parameters = c(model_degree, model_scale, model_C),
          cf_matrix = cf_matrix
        )
      
      if (return_plots) {
        plts[["SVM"]] <- cf_matrix$table %>%
          data.frame() %>%
          ggplot2::ggplot(aes(Prediction, Reference)) +
          geom_tile(aes(fill = Freq), colour = "gray50") +
          scale_fill_gradient(low = "beige", high = muted("chocolate")) +
          geom_text(aes(label = Freq)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("SVM")
      }
      toc()
    }
    
    # -------------------------- DT ------------------------
    if ("DT" %in% model_names) {
      tic("DT took")
      message("Starting DT")
      model_name <- "C5.0"
      train_control_method <- "none"
      model_trails <- 10
      model_model <- "C5.0"
      model_winnow <- FALSE
      
      model_A <- train(
        trimmed_activity ~ .,
        data = training_df,
        method = model_name,
        trControl = fitControl,
        verbose = FALSE,
        tuneGrid = data.frame(
          trials = model_trails,
          model = model_model,
          winnow = model_winnow
        ),
        metric = "ROC"
      )
      
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # To calculate area AUC we need probabilies and predicted classes in a single dataframe
      pred_prob <-
        data.frame(obs =  testing_df$trimmed_activity,
                   pred = pred)
      pred <-
        stats::predict(model_A, newdata = testing_df, type = "prob")
      pred_prob <- bind_cols(pred_prob, pred)
      
      # Calculate different metrics
      metrics <-
        multiClassSummary(data = pred_prob,
                          lev = levels(testing_df$trimmed_activity)) %>%
        as.data.frame()
      # Return the metric in a nicer format
      metric_names <- rownames(metrics)
      metrics  %<>% data.table::transpose()
      colnames(metrics) <- metric_names
      rownames(metrics) <- "DT"
      accuracies  %<>%  rbind(metrics)
      
      # CF need a different format of prediction results so recalcuate
      pred <- stats::predict(model_A, newdata = testing_df)
      
      # Calculate confusion matrix
      cf_matrix <-
        confusionMatrix(
          data = pred,
          reference = testing_df$trimmed_activity,
          mode = "prec_recall"
        )
      
      
      
      # Create a list of the model and the results to save
      results[["DT"]] <-
        list(
          split_seed = seed,
          model_name = model_name,
          model = model_A,
          train_control_method = train_control_method,
          tune_parameters = c(model_trails, model_model, model_winnow),
          cf_matrix = cf_matrix
        )
      
      if (return_plots) {
        plts[["DT"]] <- cf_matrix$table %>%
          data.frame() %>%
          ggplot2::ggplot(aes(Prediction, Reference)) +
          geom_tile(aes(fill = Freq), colour = "gray50") +
          scale_fill_gradient(low = "beige", high = muted("chocolate")) +
          geom_text(aes(label = Freq)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("DT")
      }
      toc()
    }
    
    
    # ---------------------------- save the results ----------------------
    
    if (save_results_on_disk) {
      fname <- paste0("Model_results_", as.numeric(now()), ".RData")
      save(results,
           file = fname)
      message(paste0("The models are stored in ", fname))
    }
    
    
    if (cores != 1) {
      # To stop parallel calculation
      stopCluster(cl)
    }
    
    
    
    output <-
      list(
        "Model-Accuracy" = accuracies,
        "Plots" = plts,
        "RF_Feature_Importance" = importance
      )
    return(output)
  }