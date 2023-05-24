funcNorm <- function(x, y, z) return(sqrt(x^2 + y^2 + z^2))
funcVectorDotFeatures <- function(df, 
                                  removeDiff=FALSE, #remove delta variables or not
                                  deltaGenerate,  #generate delta varaibles or not (make it false if thesse variables are generated already)
                                  ori = c("x", "y", "z")){ #names of axis variables
  if(deltaGenerate){  
    # Norm of first order diff and secon order diff
    diffVar <- paste( ori, "d1", sep="_")
    df[,  (diffVar) := data.table::shift(.SD, 1, fill = NA, "lead") - .SD, .SDcols=ori, by = c("id", "loc")]
    
    diffVar2 <- paste(ori, "d2", sep="_")
    df[,  (diffVar2) := data.table::shift(.SD, 1, fill = NA, "lead") - .SD, .SDcols=diffVar, by = c("id", "loc")]
  }
  
  # Get l2 norms 
  df[, en := funcNorm(x, z, y)]
  
  df[, en_d1 := funcNorm(x_d1, y_d1, z_d1)]
  df[, en_d2 := funcNorm(x_d2, y_d2, z_d2)]
  
  # Remove variables
  if(removeDiff){
    df[, c(diffVar, diffVar2):= NULL]
  }
  
  # ANGLE - two ways to calculate
  leadVar <- paste(ori, "lead", sep="_")
  #df[, ang1 := (x*x_d1 + y*y_d1 + z*z_d1)/(en*en_d1), ]
  df[, (leadVar) := data.table::shift(.SD, 1, fill = NA, "lead"), .SDcols=ori, by = c("id", "loc")]
  df[, en_lead := data.table::shift(en, 1, fill = NA, "lead"), by = c("id", "loc")]
  head(df)
  df[, ang1 := acos((x*x_lead + y*y_lead + z*z_lead)/(en*en_lead)), by = c("id", "loc")]
  df[, ang2 := acos((x/en)*(x_lead/en_lead) + (y/en)*(y_lead/en_lead) + (z/en)*(z_lead/en_lead)), by = c("id", "loc")]

  # Invalid values of ang due to zero denominator, to be imputated later 
  df[, sum(is.na(ang1))]
  df[, sum(is.na(ang2))]
  #df<-na.omit(d, cols = c("ang1", "ang2"))
  nrow(df)
  if(removeDiff){
    df[, c(leadVar, "en_lead"):= NULL]
  }
  return(df)
  ### End of Biosensor feature prep on raw (disaggregated, resampled) data 
}
