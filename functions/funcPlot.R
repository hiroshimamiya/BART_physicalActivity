



### Plot zoomed time-series for each class and each axis  (3 axis and 6 classes)------
# function to plot xyz by class
funcPlotAxis <- function(d){
  d$seq <- 1:nrow(d)
  className <- d$class
  d <- melt(d, id.vars = c("seq", "class"), measure.vars = c("x", "y", "z"))
  p <- ggplot(d, aes(y = value, x = seq)) +
    geom_line(size = 0.1) + 
    facet_wrap(~variable, ncol = 1) + 
    ggtitle(paste("CLASS: ", className, ":    TimePoint(disaggregated, not resampled)= ", timePoint[1], "...", timePoint[2])) + 
    theme_minimal()
  p$class <- unique(d$class)
  p$title <- paste("RawPlot", p$class, plotPersonID, plotDevLocation, "time_", timePoint[1], ":", timePoint[2], sep = "_")
  return(p)
}

# function to save pltos above 
funcPlotSave <- function(p, path, w = 9, h = 4.5){
  f <- paste(path, p$title, ".jpg", sep = "")
  ggsave(filename = f, plot = p, width = w, height = h)
}

#p <- lapply(sampDataList, funcPlotAxis)
#imagePath = "~/scratch/BART_test/image/descriptive/axisUnsampled/"
#g <-arrangeGrob(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], ncol = 1)
#g <-arrangeGrob(grobs = p, ncol = 1)
#if(bool_plot){
#  ggsave(file = paste(imagePath, "Plot of each Class__PersonID", plotPersonID, plotDevLocation, "timeRange_", timePoint[1], ":", timePoint[2], ".pdf", sep = "_" ), g, height = 20, width = 15)
#}


funcPlotClassByColor <- function(df){                          
  p <- ggplot(df, aes(y = value, x = seq, color= factor(class))) +
    geom_path(aes(group = 1), size = 0.1) +
    facet_wrap(~variable, scales = "free", ncol = 1) +  
    theme_minimal() + 
    theme(strip.text.x = element_text(size = 60)) + 
    theme(text = element_text(size = 60)) +  
    ylab("Raw value") + xlab("Raw (not resampled) time-series ") + 
    scale_x_continuous(breaks = round(seq(min(df$seq), max(df$seq), by = 5000),0)) + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+ # graph order lines 
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  return(p)
}



funcGetSamplePlot <- function(fitObj, title, obsID ){
  df <- data.frame(fitObj$post$prob.test[,(obsID*6 +1):(obsID*6+6)])
  #df <- data.frame(fitObj$post$prob.test[,7:13])
  className <- c("Lying", "Sit", "Walk", "Run3", "Run5", "Run7")
  colnames(df) <- className
  df$MCMCsample <- 1:nrow(df)
  #rename_with(~ gsub("X", "Class", .x, fixed = TRUE))
  tail(df)
  
  #df <- df %>% gather(class, sample, -MCMCsample)
  df <- melt(setDT(df), id.vars =  "MCMCsample", measure.vars = className, variable.name = "class", value.name = "sample")
  
  ggplot(data = df, aes(x = MCMCsample, y = sample)) + 
    geom_line()+ 
    facet_wrap(~class, nrow=1) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank()
    ) + 
    ggtitle(title)
}



funcPlotHz <- function(d, location){
  t <- d[loc == location, .(count = .N), by = c("record_time", "id")]
  table(t$count)
  barplot(table(t$count), main="Frequency (measurement per sec), all 48 people, 
          Identical regardless of device location",
          xlab="Number of counts per person-second")
}





funcPlotListTime <- function(d, title=NA,  target= NA, devLocation = NA, pathImage=NA, h = NA, w = NA){
  
  p <- list()
  locVec <- unique(d$loc)
  idVec <- unique(d$id)

  for(i in 1:length(idVec)){
    a <- d[id == idVec[i] & loc==devLocation]
    b <- melt(setDT(a), id.vars  = 1,  measure.vars = target)
    b[,time := .(seq_len(.N)), by = variable]
    
    p[[i]] <- ggplot(b, aes(y = value, x=time))+ 
      geom_point(shape=16, size =0.05) + 
      facet_wrap(~variable, ncol=1, scales="free") + 
      theme_bw() + 
      ylab("raw value") +
      theme(panel.grid.major = element_blank()) +
      ggtitle(paste(title, idVec[i], "_location_",devLocation ,sep =""), subtitle = "Unscaled Y across plots")
    
    if(!is.na(pathImage)) ggsave(paste(pathImage, title, idVec[i], "_location_",devLocation, ".jpg",sep =""), plot = p[[i]], width = w, height = h, units = "in")
    fileName <- paste(pathImage, title, idVec[i], "_location_",devLocation, ".jpg",sep ="")
    p[[i]]$fileName <- fileName
    p[[i]]$personID <- idVec[i]
    p[[i]]$devLocation <- devLocation
    p[[i]]$target <- target
    p[[i]]$title <- title
  }
  names(p) <- idVec
  return(p)
}




# Plot by column, facet participants 
funcPlotFeatureByPerson <- function(ind, d, targetVector, nCol=2){
  targetVariable <- targetVector[ind]
  colNames <- c("id", "time", "loc", targetVariable)
  d <- d[, ..colNames]
  colnames(d)[d[, ncol(d)]] <- "X" 
  loc <- unique(d$loc)

  p <- ggplot(d, aes(y = X, x=time))+ 
    geom_point(shape=16, size =0.05) + 
    facet_wrap(~id, ncol = nCol) + 
    theme_bw() + 
    ylab("Value") +
    theme(panel.grid.major = element_blank()) + 
    ggtitle(paste("Variable_", targetVariable, "_location_", loc, sep = ""))
  return(p)
}


funcPlotSave <- function(p, path, w = 9, h = 4.5){
  f <- paste(path, p$title, "_", p$personID, "_loc_",p$devLocation, ".jpg", sep = "")
  ggsave(filename = f, plot = p, width = w, height = h)
}



#ggsave("~/scratch/BART_test/image/descriptive/pp.jpg", pl, width = 6, height = 3)
#pl <- arrangeGrob(p[[1]], p[[2]], ncol = 1)
#fileName <- paste("30Hz_XYZ_PersonID_", idVec[1], "_location_",devLocation,  ".png",sep ="")
#ggsave(paste("image/descriptive/", fileName, "_location_",devLocation, sep = ""), pl, width = 11, height =4.5)



