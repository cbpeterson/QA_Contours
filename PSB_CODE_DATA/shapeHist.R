#######################################
#Shape Features Histogram Function
#######################################
#### Good


library(R.matlab)
library(EBImage)
library(oro.dicom)
library(tidyverse)
library(radiomics)
library(glcm) 
library(shapes)
library(splancs)
library(concaveman)

folderName <- "Good_Bladder"
g_Files <- list.files("Good_Bladder")

shapeHist <- function(g_Files, folderName, sliceThickness, pixelSpace){

m.feat <- c()
feat.2d <- c()

for (k in 1:length(g_Files)) {
  dataG <- readMat(paste0(paste0(folderName,"/"),g_Files[k] ))
  
  dG <- dataG$mask[[1]]
  vG <- c()
  for (j in 1:dim(dG)[1]) {
    vG[j] <- sum(dG[j,,]) > 0
  }
  
  a <- 1:dim(dG)[1] 
  

  
  #################################################
  # Calculates perimeter better
  
  
  
  # Shift the image in one direction
  s1 <- function(z) cbind(rep(0,nrow(z)), z[,-ncol(z)] )
  s2 <- function(z) cbind(z[,-1], rep(0,nrow(z)) )
  s3 <- function(z) rbind(rep(0,ncol(z)), z[-nrow(z),] )
  s4 <- function(z) rbind(z[-1,], rep(0,ncol(z)) )
  
  # Area, perimeter and circularity
  area <- function(z) sum(z)
  perimeter <- function(z) sum( z != 0 & s1(z)*s2(z)*s3(z)*s4(z) == 0)
  circularity <- function(z) 4*pi*area(z) / perimeter(z)^2
  
  edge <- function(z) z & !(s1(z)&s2(z)&s3(z)&s4(z))
  perimeter <- function(z) {
    e <- edge(z)
    ( 
      # horizontal and vertical segments
      sum( e & s1(e) ) + sum( e & s2(e) ) + sum( e & s3(e) ) + sum( e & s4(e) ) + 
        # diagonal segments
        sqrt(2)*( sum(e & s1(s3(e))) + sum(e & s1(s4(e))) + sum(e & s2(s3(e))) + sum(e & s2(s4(e))) )
    ) / 2  # Each segment was counted twice, once for each end
  }
  #####################################################################################
  
  
  
  
  
  #Feature Vector
  name1 <- c()
  slice1 <- c()
  area1 <- c()
  peri1 <- c()
  rMean1 <- c()
  rMin1 <- c()
  rMax1 <- c()
  cent1 <- c()
  compact1 <- c()
  sphericity1 <- c()
  convexity1 <- c()
  solidity1 <- c()
  roundness1 <- c()
  
  #Histogram features for 2-D shape features
  for (i in 1:length(a[vG])) {
    
    name1[i] <- g_Files[k]
    slice1[i] <- paste0(a[vG][i], "/", max(a))
    #area1[i] <-  area(dG[a[vG][i],,])*pixelSpace*pixelSpace #Pixel spacing for area
    peri1[i] <-  perimeter(dG[a[vG][i],,])*pixelSpace #Pixel spacing
    rMean1[i] <-  computeFeatures.shape(dG[a[vG][i],,])[3]
    rMin1[i] <-  computeFeatures.shape(dG[a[vG][i],,])[5]
    rMax1[i] <-  computeFeatures.shape(dG[a[vG][i],,])[6]
    cent1[i] <- centroid.size(dG[a[vG][i],,])
    #compact1[i] <- (4*pi*area1[i])/(peri1[i]^2)
    sphericity1[i] <-  rMin1[i]/rMax1[i]
    
    #New Shape statistics
    
    # Get X, Y corrdinates
    xyCo <- which(dG[a[vG][i],,]==1, arr.ind=TRUE)
    hpts <- chull(xyCo)
    hpts <- c(hpts, hpts[1])
    
    #calculate Perimeter of convex shape for convexity
    pDist <- c()
    for(j in 1:length(xyCo[hpts, ][,1])-1 ){
      pDist[j] <- sqrt( (xyCo[hpts, ][j,1] -xyCo[hpts, ][j+1,1])^2 +  (xyCo[hpts, ][j,2] -xyCo[hpts, ][j+1,2])^2)
    }
    p_vex <- sum(pDist)*pixelSpace
    convexity1[i] <- p_vex/peri1[i]
    
    if(is.na(convexity1[i])){
      convexity1[i] <-  mean(convexity1 , na.rm=TRUE)
    }
    
    
    #calculate Solidity
    caveXY<- concaveman(xyCo)
    convexArea <- as.matrix(xyCo[hpts, ]) %>% areapl
    regularArea <- as.matrix(caveXY) %>% areapl
    
    area1[i] <- regularArea*pixelSpace*pixelSpace
    compact1[i] <- (4*pi*area1[i])/(peri1[i]^2)
    solidity1[i] <- regularArea/convexArea
    
    #calculate Roundness
    roundness1[i] <- (4*pi*area1[i])/(p_vex^2)
    if(roundness1[i] == Inf){
      roundness1[i] <- mean(roundness1 , inf.rm=TRUE)
    }
    if(is.na(roundness1[i])){
      roundness1[i] <- mean(roundness1 , na.rm=TRUE)
    }
  }
  
  feat1 <- cbind(area1, peri1, rMean1, rMin1, rMax1, cent1, compact1, sphericity1, convexity1, solidity1, roundness1)
  feat.2d <- rbind(feat.2d, cbind(name1, slice1, feat1))
  
  areaH <- c(summary(area1),sd(area1))
  periH <- c(summary(peri1),sd(peri1))
  rMean1H <- c(summary(rMean1),sd(rMean1))
  rMin1H <- c(summary(rMin1),sd(rMin1))
  rMax1H <- c(summary(rMax1),sd(rMax1))
  cent1H <- c(summary(cent1),sd(cent1))
  compact1H <- c(summary(compact1),sd(compact1))
  sphericity1H <- c(summary(sphericity1),sd(sphericity1))
  covexity1H <-  c(summary(convexity1),sd(convexity1))
  solidity1H <- c(summary(solidity1),sd(solidity1))
  roundness1H <- c(summary(roundness1),sd(roundness1))
  
  ####################################################
  feat2  <-  c(areaH, periH,rMean1H, rMin1H ,rMax1H, cent1H, compact1H, sphericity1H, covexity1H, solidity1H, roundness1H)
  
  feat2[78] <- sum(feat1[,2])*sliceThickness
  feat2[79] <- sum(dG)*sliceThickness*pixelSpace*pixelSpace
  feat2[80] <- feat2[78]/feat2[79]
  
  
  
  names(feat2) <-   c(paste0(rep("area", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("peri", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("rMean", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("rMin", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("rMax", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("centroid", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("compact", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("sphericity", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("convexity", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("solidity", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      paste0(rep("roundness", 7), c("Min","1Q", "Median", "Mean","3Q","Max","SD")),
                      "surface.area","volume", "SA_to_V")
  m.feat <- rbind(m.feat,feat2)
  print(k)
  listFeatures <- list(m.feat,feat.2d )
}
return(listFeatures)
}



