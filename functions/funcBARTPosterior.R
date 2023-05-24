# Function to re-oirganize posterior sample, from wide to long , attach class label, and find class with max probability at each data point and each mcmc iteration 

funcSummrizePost <- function(bartObj, xData, yData, cores, dataID, classNum = 6){
  classN <- classNum # Number of classes
  if(nrow(xData) != length(yData))(stop("Data mismatch"))
  # generate posterior sample of prediction 
  postData <- predict(bartObj, xData, mc.cores = cores) 
  # number of iterations 
  mcmcN <- nrow(postData$prob.test) 
  # number of class*observations - probabilities of each class for each time point
  paramN <- ncol(postData$prob.test)
  # prob.test is a matrix, with row number of iterations and column for class probability and number of observation
  # Convert prob.test to list , where each list is time point and thus matrix of mcmcN*classN
  splitPostProb <- split.default(data.frame(postData$prob.test), rep(1:(paramN/classN), each = classN))
  # Attach column names again
  splitPostProb <- lapply(splitPostProb, function(x){colnames(x) <- paste("Class", c(1:classN), sep = ""); return(x)})
  # Combine all data vertically and convert to DT (data.table)  
  splitPostProbDT <- bind_rows(splitPostProb, .id = "obsID")
  setDT(splitPostProbDT)

  # Attach observed (true) Y, extract columns for class probabilities, get column name with max probability and and rename, get max probabilities as new variable, create variables for matching status of observed and estimated, and give mcmc iteration id 
  splitPostProbDT[, obsClass := rep(yData, each=mcmcN)]
  splitPostProbDT[, obsClassLabel := factor(obsClass)]
  levels(splitPostProbDT$obsClassLabel) <- c("Lying", "Sitting", "Walk", "Run3", "Run5", "Run7")
  colIndex <- grep("Class\\d", colnames(splitPostProbDT)) # Column index of class probabilities 
  # Class with max probability from each classes at each time point
  splitPostProbDT[, maxProbClass := names(.SD)[max.col(.SD)], .SDcols = colIndex]
  # Extract class ID (1, 6) from columsn name that have max probability for each time point and iteraction 
  splitPostProbDT[, maxProbClass := as.integer(gsub("Class","",maxProbClass))]
  # Make a column listing max prbability at each iteration 
  splitPostProbDT[, maxProbValue := do.call(pmax, .SD), .SDcols = colIndex]
  # Give an indicator, whether observed and predicted class (latter is class with max prbability) are matched or not 
  splitPostProbDT[, matchMCMC := ifelse(maxProbClass == obsClass, 1, 0)]
  # Give iteration ID 
  splitPostProbDT[, mcmcItr := rep( 1:mcmcN, (paramN/classN))]
  
  # Attach time ID 
  if(nrow(splitPostProbDT) != length(rep(seq_len(nrow(dataID)), each = mcmcN)))(stop("Data mismatch"))
  splitPostProbDT <- splitPostProbDT %>% bind_cols( 
    dataID[rep(seq_len(nrow(dataID)), each = mcmcN), ] 
  )
  return(splitPostProbDT)
}



funcPlotHistPost <- function(timeid, personid, df, bin, yMax = 300){
  ti <- df[obsID == timeid & id == personid, ]
  ti_long <- melt(ti, 
                  id.vars =  c("mcmcItr", "obsID", "obsClass", "obsClassLabel"),
                  measure.vars = className, 
                  variable.name = "class", 
                  value.name = "Probability")
  
  levels(ti_long$class) <- c("Lying", "Sitting", "Walking", "Running_3METs", "Running_5METs", "Running_7METs")
  
  ggplot(ti_long, aes(x=Probability, fill=class)) + 
    geom_histogram(alpha=0.5, position="identity", binwidth = bin) + 
    xlim(-0.1,1.1) + 
    coord_cartesian(ylim=c(0, yMax)) + 
    theme_classic() + 
    theme(
      text = element_text(size = 15), 
      plot.title = element_text(size = 15), 
      legend.text = element_text(size = 15), 
      legend.title = element_text(size = 15)
    ) +
    scale_fill_discrete(name = "Activitiy types & \n 3 intensities")
}

# Get Plot of each observation, for 6 classes from data frame of posterior samples 
funcGetSamplePlot2 <- function(postDf, title, obsID, classVec){
  trueClass <- classVec[obsID]
  obsID <- obsID-1
  df <- data.frame(postDf[,(obsID*classN
                            +1:classN
  )])
  className <- c("Lying", "Sit", "Walk", "Run3", "Run5", "Run7")
  colnames(df) <- className
  df$MCMCsample <- 1:nrow(df)
  df <- melt(setDT(df), id.vars =  "MCMCsample", measure.vars = className, variable.name = "class", value.name = "sample")
  ggplot(data = df, aes(x = MCMCsample, y = sample)) + 
    geom_line()+ 
    facet_wrap(~class, nrow=1) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    ) + 
    ggtitle(paste(title, "CLASS", trueClass))
}

# get post prob summary and 
funcGetPostMean <- function(postVec, title, obsID, classVec){
  trueClass <- classVec[obsID]
  obsID <- obsID-1
  v <- c(postVec[(obsID*6 +1:6)])
  a <- rep(FALSE, 6); a[trueClass] <- TRUE 
  c <- rep(FALSE, 6); 
  c[which.max(v)] <- max(v)
  rbind(v, c, a)
}

