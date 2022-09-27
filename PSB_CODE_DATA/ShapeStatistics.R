# Streamlined code for PSB Sumbission

### 3d contour meshes
library(R.matlab)
library(EBImage)
library(oro.dicom)
library(tidyverse)
library(radiomics)
library(glcm) 
library(shapes)
library(misc3d)
library(rgl)
library(alphashape3d)
library(plot3D)
library(precrec)
library(ranger)
library(shapr)

## Random Forest Function
library(randomForest)
library(ggplot2)
library(tidyverse)
library(glmnet)
library(caret)
library(naivebayes)
library(xgboost)

##########################################################################################################################################################################
##########################################################################################################################################################################
# Necessary Functions


#Plot 3D organ mesh, input data
organ3D <- function(data, color="black"){
  d1 <- dim(data$mask[[1]])[1]
  d2 <- data$mask[[1]]
  open3d()
  contour3d(d2, 1.5*(1:d1), 1:512, 1:512, lev=c(0, 1),
            alpha = c(0.2,  1), color = c("white", color))
}


# Comparing different Models:

CVshapeRF <- function(gkidney,bkidney, seed=314){
  
  
  gkidney[is.na(gkidney)]=0
  sum(is.na(gkidney))
  sum(is.na(bkidney))
  bkidney[is.na(bkidney)]=0
  
  
  gkidney[gkidney == Inf] = 0
  bkidney[bkidney == Inf] = 0
  
  
  gg <- cbind(gkidney, 0)
  bb <- cbind(bkidney,1)
  
  
  
  names(gg)[names(gg) != names(bb)] <- "class"
  names(bb)[names(gg) != names(bb)] <- "class"
  
  d <- as.tibble(rbind(gg,bb))
  d$class <- as.factor(d$class)
  
  #######################################################
  #######################################################
  ## New random Sample Here
  ##################################################
  
  
  
  pred_vals <- NULL
  true_vals <- NULL
  
  set.seed(seed)
  smp_size <- floor(.8*nrow(gg))
  s <- sample(seq_len(nrow(gg)), size=nrow(gg))
  s2 <- sample(seq_len(nrow(bb)), size=nrow(bb))
  
  it <- floor(nrow(gg)*.2)
  it2 <- floor(nrow(bb)*.2)
  
  sMat <- c()
  s2Mat <- c()
  for (i in 1:4) {
    sMat <-  rbind(sMat,s[(1+(i-1)*it):(it*i)])
    s2Mat <- rbind(s2Mat,s2[(1+(i-1)*it2):(it2*i)])
  }
  
  
  for (i in 1:4) {
    test <- as.tibble(rbind(gg[sMat[i,],], bb[s2Mat[i,],]))
    train <- as.tibble(rbind(gg[-sMat[i,],], bb[-s2Mat[i,],]))
    
    m1 <- randomForest(class ~ ., data = train, ntree=500,
                       mtry=16, importance = TRUE)
    
    #plot(m1, main="Random Forest")
    
    p1 <- predict(m1, test, type="class" )
    
    
    
    pred_vals = c(pred_vals, p1)
    true_vals = c(true_vals, test$class)
  }
  
  
  
  test <- as.tibble(rbind(gg[-s[1:(4*it)],], bb[-s2[1:(4*it2)],]))
  train <- as.tibble(rbind(gg[s[1:(4*it)],], bb[s2[1:(4*it2)],]))
  
  
  m1 <- randomForest(class ~ ., data = train, ntree=500,
                     mtry=16, importance = TRUE)
  
  
  
  p1 <- predict(m1, test, type="class" )
  
  
  pred_vals = c(pred_vals, p1)
  true_vals = c(true_vals, test$class)
  totvals=cbind(pred_vals, true_vals)
  
  return(totvals)
  
  
}



CVshapeLogReg <- function(gkidney,bkidney, seed=314){
  
  str(gkidney)
  gkidney[is.na(gkidney)]=0
  sum(is.na(gkidney))
  sum(is.na(bkidney))
  bkidney[is.na(bkidney)]=0
  
  sum(bkidney == Inf)
  sum(gkidney == Inf)
  
  gkidney[gkidney == Inf] = 0
  bkidney[bkidney == Inf] = 0
  
  
  gg <- cbind(gkidney, 0)
  bb <- cbind(bkidney,1)
  
  
  
  names(gg)[names(gg) != names(bb)] <- "class"
  names(bb)[names(gg) != names(bb)] <- "class"
  
  d <- as.tibble(rbind(gg,bb))
  d$class <- as.factor(d$class)
  
  #######################################################
  #######################################################
  ## New random Sample Here
  ##################################################
  
  pred_vals <- NULL
  true_vals <- NULL
  
  set.seed(seed)
  smp_size <- floor(.8*nrow(gg))
  s <- sample(seq_len(nrow(gg)), size=nrow(gg))
  s2 <- sample(seq_len(nrow(bb)), size=nrow(bb))
  
  it <- floor(nrow(gg)*.2)
  it2 <- floor(nrow(bb)*.2)
  
  sMat <- c()
  s2Mat <- c()
  for (i in 1:4) {
    sMat <-  rbind(sMat,s[(1+(i-1)*it):(it*i)])
    s2Mat <- rbind(s2Mat,s2[(1+(i-1)*it2):(it2*i)])
  }
  
  
  for (i in 1:4) {
    test <- as.tibble(rbind(gg[sMat[i,],], bb[s2Mat[i,],]))
    train <- as.tibble(rbind(gg[-sMat[i,],], bb[-s2Mat[i,],]))
    
    m1 <- glm(class~., data = train, family = "binomial")
    
    p1 <- predict(m1, test)
    p1 <- exp(p1)/(1+exp(p1))
    
    
    
    ###############################
    # Threshold value for ROC Curve
    ## Warning: package 'ROCit' was built under R version 3.5.2
    
    
    pred_vals = c(pred_vals, p1)
    true_vals = c(true_vals, test$class)
  }
  
  
  
  test <- as.tibble(rbind(gg[-s[1:(4*it)],], bb[-s2[1:(4*it2)],]))
  train <- as.tibble(rbind(gg[s[1:(4*it)],], bb[s2[1:(4*it2)],]))
  
  #####################################################
  m1 <- glm(class~., data = train, family = "binomial")
  
  p1 <- predict(m1, test)
  p1 <- exp(p1)/(1+exp(p1))
  
  
  
  
  pred_vals = c(pred_vals, p1)
  true_vals = c(true_vals, test$class)
  totvals=cbind(pred_vals, true_vals)
  
  return(totvals)
}


CVshapeLasso <- function(gkidney,bkidney, seed=314){
  
  str(gkidney)
  gkidney[is.na(gkidney)]=0
  sum(is.na(gkidney))
  sum(is.na(bkidney))
  bkidney[is.na(bkidney)]=0
  
  sum(bkidney == Inf)
  sum(gkidney == Inf)
  
  gkidney[gkidney == Inf] = 0
  bkidney[bkidney == Inf] = 0
  
  
  gg <- cbind(gkidney, 0)
  bb <- cbind(bkidney,1)
  
  
  
  names(gg)[names(gg) != names(bb)] <- "class"
  names(bb)[names(gg) != names(bb)] <- "class"
  
  d <- as.tibble(rbind(gg,bb))
  d$class <- as.factor(d$class)
  
  #######################################################
  #######################################################
  ## New random Sample Here
  ##################################################
  
  pred_vals <- NULL
  true_vals <- NULL
  
  set.seed(seed)
  smp_size <- floor(.8*nrow(gg))
  s <- sample(seq_len(nrow(gg)), size=nrow(gg))
  s2 <- sample(seq_len(nrow(bb)), size=nrow(bb))
  
  it <- floor(nrow(gg)*.2)
  it2 <- floor(nrow(bb)*.2)
  
  sMat <- c()
  s2Mat <- c()
  for (i in 1:4) {
    sMat <-  rbind(sMat,s[(1+(i-1)*it):(it*i)])
    s2Mat <- rbind(s2Mat,s2[(1+(i-1)*it2):(it2*i)])
  }
  
  
  for (i in 1:4) {
    test <- as.tibble(rbind(gg[sMat[i,],], bb[s2Mat[i,],]))
    train <- as.tibble(rbind(gg[-sMat[i,],], bb[-s2Mat[i,],]))
    
    lambdas <- 10^seq(2, -3, by = -.1)
    
    x <- as.matrix(train[,-length(names(train))])
    y <- as.matrix(train[,length(names(train))])
    
    
    # Setting alpha = 1 implements lasso regression
    lasso_reg <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5, family="binomial")
    
    # Best 
    lambda_best <- lasso_reg$lambda.min 
    
    lasso_model <- glmnet(x, y, alpha = 1, lambda = lambda_best, standardize = TRUE)
    #summary(lasso_model$beta)
    
    p1 <-  predict(lasso_model, s = lambda_best, newx =  as.matrix(test[,-length(names(train))]))
    p1 <- exp(p1)/(1+exp(p1))
    
    
    
    
    pred_vals = c(pred_vals, p1)
    true_vals = c(true_vals, test$class)
  }
  
  
  
  test <- as.tibble(rbind(gg[-s[1:(4*it)],], bb[-s2[1:(4*it2)],]))
  train <- as.tibble(rbind(gg[s[1:(4*it)],], bb[s2[1:(4*it2)],]))
  
  #####################################################
  
  #####################
  ## Lasso
  
  lambdas <- 10^seq(2, -3, by = -.1)
  
  x <- as.matrix(train[,-length(names(train))])
  y <- as.matrix(train[,length(names(train))])
  
  
  # Setting alpha = 1 implements lasso regression
  lasso_reg <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5, family="binomial")
  
  # Best 
  lambda_best <- lasso_reg$lambda.min 
  
  lasso_model <- glmnet(x, y, alpha = 1, lambda = lambda_best, standardize = TRUE)
  #summary(lasso_model$beta)
  
  p1 <-  predict(lasso_model, s = lambda_best, newx =  as.matrix(test[,-length(names(train))]))
  p1 <- exp(p1)/(1+exp(p1))
  
  
  
  
  pred_vals = c(pred_vals, p1)
  true_vals = c(true_vals, test$class)
  totvals=cbind(pred_vals, true_vals)
  
  return(totvals)
}


CVshapeNaiveBayes <- function(gkidney,bkidney, seed=314){
  
  gkidney[is.na(gkidney)]=0
  sum(is.na(gkidney))
  sum(is.na(bkidney))
  bkidney[is.na(bkidney)]=0
  
  sum(bkidney == Inf)
  sum(gkidney == Inf)
  
  gkidney[gkidney == Inf] = 0
  bkidney[bkidney == Inf] = 0
  
  
  gg <- cbind(gkidney, 0)
  bb <- cbind(bkidney,1)
  
  
  
  names(gg)[names(gg) != names(bb)] <- "class"
  names(bb)[names(gg) != names(bb)] <- "class"
  
  d <- as.tibble(rbind(gg,bb))
  d$class <- as.factor(d$class)
  
  #######################################################
  #######################################################
  ## New random Sample Here
  ##################################################
  
  pred_vals <- NULL
  true_vals <- NULL
  
  set.seed(seed)
  smp_size <- floor(.8*nrow(gg))
  s <- sample(seq_len(nrow(gg)), size=nrow(gg))
  s2 <- sample(seq_len(nrow(bb)), size=nrow(bb))
  
  it <- floor(nrow(gg)*.2)
  it2 <- floor(nrow(bb)*.2)
  
  sMat <- c()
  s2Mat <- c()
  for (i in 1:4) {
    sMat <-  rbind(sMat,s[(1+(i-1)*it):(it*i)])
    s2Mat <- rbind(s2Mat,s2[(1+(i-1)*it2):(it2*i)])
  }
  
  
  for (i in 1:4) {
    test <- as.tibble(rbind(gg[sMat[i,],], bb[s2Mat[i,],]))
    train <- as.tibble(rbind(gg[-sMat[i,],], bb[-s2Mat[i,],]))
    
    x <- as.matrix(train[,-length(names(train))])
    y <- as.matrix(train[,length(names(train))])
    y <- as.factor(y)
    m1 <- naive_bayes(x,y, usekernel = T)
    p1 <- predict(m1, test, type="prob")
    p1 <- p1[,2]
    
    
    
    ###############################################################
    pred_vals = c(pred_vals, p1)
    true_vals = c(true_vals, test$class)
  }
  
  
  
  test <- as.tibble(rbind(gg[-s[1:(4*it)],], bb[-s2[1:(4*it2)],]))
  train <- as.tibble(rbind(gg[s[1:(4*it)],], bb[s2[1:(4*it2)],]))
  
  #####################################################
  x <- as.matrix(train[,-length(names(train))])
  y <- as.matrix(train[,length(names(train))])
  y <- as.factor(y)
  m1 <- naive_bayes(x,y, usekernel = T)
  p1 <- predict(m1, test, type="prob")
  p1 <- p1[,2]
  
  
  pred_vals = c(pred_vals, p1)
  true_vals = c(true_vals, test$class)
  totvals=cbind(pred_vals, true_vals)
  
  return(totvals)
}


CVxgboost <- function(gkidney,bkidney, seed=314){
  
  gkidney[is.na(gkidney)]=0
  sum(is.na(gkidney))
  sum(is.na(bkidney))
  bkidney[is.na(bkidney)]=0
  
  sum(bkidney == Inf)
  sum(gkidney == Inf)
  
  gkidney[gkidney == Inf] = 0
  bkidney[bkidney == Inf] = 0
  
  
  gg <- cbind(gkidney, 0)
  bb <- cbind(bkidney,1)
  
  
  
  names(gg)[names(gg) != names(bb)] <- "class"
  names(bb)[names(gg) != names(bb)] <- "class"
  
  d <- as.tibble(rbind(gg,bb))
  d$class <- as.factor(d$class)
  
  #######################################################
  #######################################################
  ## New random Sample Here
  ##################################################
  
  pred_vals <- NULL
  true_vals <- NULL
  
  set.seed(seed)
  smp_size <- floor(.8*nrow(gg))
  s <- sample(seq_len(nrow(gg)), size=nrow(gg))
  s2 <- sample(seq_len(nrow(bb)), size=nrow(bb))
  
  it <- floor(nrow(gg)*.2)
  it2 <- floor(nrow(bb)*.2)
  
  sMat <- c()
  s2Mat <- c()
  for (i in 1:4) {
    sMat <-  rbind(sMat,s[(1+(i-1)*it):(it*i)])
    s2Mat <- rbind(s2Mat,s2[(1+(i-1)*it2):(it2*i)])
  }
  
  
  for (i in 1:4) {
    test <- as.tibble(rbind(gg[sMat[i,],], bb[s2Mat[i,],]))
    train <- as.tibble(rbind(gg[-sMat[i,],], bb[-s2Mat[i,],]))
    
    train2 <- as.matrix(train)
    test2 <- as.matrix(test)
    
    xgb.train = xgb.DMatrix(data=train2[,-dim(train)[2]], label=train2[,dim(train)[2]])
    xgb.test = xgb.DMatrix(data=test2[,-dim(test)[2]], label=test2[,dim(test)[2]])
    
    params = list(
      booster="gbtree",
      eta=0.01,
      max_depth=12,
      gamma=3,
      subsample=0.75,
      colsample_bytree=1,
      objective="multi:softprob",
      eval_metric="auc",
      num_class=2
    )
    
    xgb.fit=xgb.train(
      params=params,
      data=xgb.train,
      nrounds=1000,
      nthreads=1,
      early_stopping_rounds=10,
      watchlist=list(val1=xgb.train,val2=xgb.test),
      verbose=0
    )
    
    xgb.pred = predict(xgb.fit,xgb.test,reshape=T)
    
    p1 <- xgb.pred[,2]
    
    
    ###############################################################
    pred_vals = c(pred_vals, p1)
    true_vals = c(true_vals, test$class)
  }
  
  
  
  test <- as.tibble(rbind(gg[-s[1:(4*it)],], bb[-s2[1:(4*it2)],]))
  train <- as.tibble(rbind(gg[s[1:(4*it)],], bb[s2[1:(4*it2)],]))
  train2 <- as.matrix(train)
  test2 <- as.matrix(test)
  
  xgb.train = xgb.DMatrix(data=train2[,-dim(train)[2]], label=train2[,dim(train)[2]])
  xgb.test = xgb.DMatrix(data=test2[,-dim(test)[2]], label=test2[,dim(test)[2]])
  
  params = list(
    booster="gbtree",
    eta=0.01,
    max_depth=12,
    gamma=3,
    subsample=0.75,
    colsample_bytree=1,
    objective="multi:softprob",
    eval_metric="auc",
    num_class=2
  )
  
  xgb.fit=xgb.train(
    params=params,
    data=xgb.train,
    nrounds=1000,
    nthreads=1,
    early_stopping_rounds=10,
    watchlist=list(val1=xgb.train,val2=xgb.test),
    verbose=0
  )
  
  xgb.pred = predict(xgb.fit,xgb.test,reshape=T)
  
  p1 <- xgb.pred[,2]
  
  
  pred_vals = c(pred_vals, p1)
  true_vals = c(true_vals, test$class)
  totvals=cbind(pred_vals, true_vals)
  
  return(totvals)
}



##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################



# Read in Shape Feature Histogram data files
gkidney <-  read.csv("good_kidney_hist_PSB.csv")
bkidney <-  read.csv("bad_kidney_hist_PSB.csv")


# Data Cleaning to replace NA's with mean value

which(is.na(gkidney))
colSums(is.na(gkidney))
gkidney[is.na(gkidney)]
which(is.na(bkidney))
gkidney$convexityMin[is.na(gkidney$convexityMin)] <- mean(gkidney$convexityMin , na.rm=TRUE)


# Add class vlaues, 0 for acceptable contours, 1 for unacceptable contours
gg <- cbind(gkidney, 0)
bb <- cbind(bkidney,1)
names(gg)[names(gg) != names(bb)] <- "class"
names(bb)[names(gg) != names(bb)] <- "class"

# Run 5-fold cross validation to compare models

# Random Forest
CV1 = CVshapeRF(gkidney,bkidney, seed=1)

# Logistic Regression
CV2 = CVshapeLogReg(gkidney,bkidney, seed=1)

# Lasso Logistic Regression
CV3 = CVshapeLasso(gkidney,bkidney, seed=1)

# Naive Bayes
CV4 = CVshapeNaiveBayes(gkidney,bkidney, seed=1)

# Extreme Gradient Boosting
CV5 = CVxgboost(gkidney,bkidney, seed=1)


# Plot the ROC and PR curves
msmdat2 <- mmdata(cbind(CV1[,1], CV2[,1], CV3[,1], CV4[,1], CV5[,1]),
                  c(CV1[,2]), modnames = c("Random Forest", "Logistic Regression", "Lasso", "Naive Bayes", "XGBoost"))
mscurves <- evalmod(msmdat2)
autoplot(mscurves)

# Get the AUC values
aa <- auc(mscurves)
aa$aucs

#####################################################################################################################################################################
####################################################################################################################################
# Training full model for Shapely Value Analysis
fullData <- as.tibble(rbind(gg,bb))
fullData <- fullData[, 1:81]


#Train Model
set.seed(314)
m1 <- randomForest(class ~ ., data = fullData, ntree=500,
                   mtry=16, importance = TRUE)

#Plot showing full model prediciton
plot(fullData$class, main=paste("Random Forest Flagging Probability"),
     ylab="Class (Probability)", xlab = "Contour Index", pch=16, cex.lab = 1.4,
     cex.main = 2,
     cex.sub = 2)
abline(h=0.5,lwd=3, col="darkgrey", lty=2)
abline(h=0.25,lwd=4, col="darkgrey")
lines(m1$predicted,  type="p", col="dodgerblue3", pch=1, cex=1.5)
legend("topleft",legend=c("True Class","Predicted Class"),
       pch = c(16,1 ), col=c("black","dodgerblue3"), cex=1)
legend(-11.5, 0.875,legend=c("0.50 Threshold","0.25 Threshold"),
       lty = c(2, 1), col=c("darkgrey","darkgrey"), cex=1, lwd=3)



# Read in unlabeled data set and plot the predicted probabilites
unlabeled <-  read.csv("unlab_kidney_hist_PSB.csv")


# build top ten model for shapely value analysis
topten <- sort(m1$importance[,2], decreasing = TRUE)[1:10]
topTen2 <- which(m1$importance[,2] >= min(topten)) # gets top variables by importance
class <- fullData$class
data2 <- cbind(fullData[topTen2], class )



m3 <- ranger(class ~ ., data = data2, num.trees = 500, mtry=8, probability = TRUE)

p3 <-  predict(m3, unlabeled)

plot(p3$predictions[,2], main=paste("Unlabeled Contour Flagging Probability"),
     ylab="Class (Probability)", col="dodgerblue3", cex=1.5,
     xlab = "Contour Index", pch=16, cex.lab = 1.4,
     cex.main = 2,
     cex.sub = 2, ylim=c(0,1))
abline(h=0.5,lwd=3, col="darkgrey", lty=2)
lines(p3$predictions[,2],  type="p", col="dodgerblue3", pch=16, cex=1.5)
legend("topleft",legend=c("Predicted Class"),
       pch = c(16 ), col=c("dodgerblue3"), cex=1)
legend(-.395, 0.925,legend=c("0.50 Threshold"),
       lty = c(2), col=c("darkgrey"), cex=1, lwd=3)

# Find which unlabeled contours get flagged

flaggedID <- which(p3$predictions[,2]>0.5)


#Find shapley values now
explainer <- shapr(data2[,-11], m3)

unlabeled2 <- unlabeled[topTen2] # Get unlabled dataset with only top variables

explain4 <- explain(unlabeled2[flaggedID, ], explainer,
                    approach = "empirical",
                    prediction_zero = 0, n_samples = 1e2
)
pp4 <- plot(explain4)
pp4 +
  ggtitle("")+
  theme(text=element_text(size=16), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=16), #change font size of axis titles
        plot.title=element_text(size=30), #change font size of plot title
        legend.text=element_text(size=16), #change font size of legend text
        legend.title=element_text(size=16)) #change font size of legend title   



#################################################################################################################################################
#################################################################################################################################################
#Histograms of shape features

badk <- read_csv("bad_kidney_2D_PSB.csv")
goodk <- read_csv("good_kidney_2D_PSB.csv")
#view(goodk)

#Data cleaning
colSums(is.na(goodk))
goodk[is.na(goodk)] <- mean(goodk$sphericity1, na.rm=TRUE)
colSums(is.na(badk))

set_good <- split(goodk,goodk$name1)
set_bad <- split(badk,badk$name1)


###### Comparing Area between good and bad contour
x1 <- set_bad[52]$BAD_Kidney_R_mask_DJCervix_CervixFinalST0140.mat$area1
x2 <- set_good[260]$Kidney_R_mask_GT_CervixFinalST0140.mat$area1

t1 <- as_tibble(cbind(c(x1,x2), c(rep(1,length(x1)), rep(2,length(x2)))))
t1$V2 <- as.factor(t1$V2)
names(t1) <- c("Area", "Class")

ggplot(t1, aes(x = Area, fill = Class, colour = Class)) + 
  geom_histogram(aes(y = ..density..),
                 alpha = 0.25,  bins=20, position = "dodge") +
  geom_density(alpha = 0.5, lwd = 1.2,
               linetype = 2)+
  scale_fill_manual(values=c("red1", "green1", "#56B4E9"),
                    labels=c("Unacceptable", "Acceptable"))+
  scale_color_manual(values=c("red4", "green4", "#56B4E9"),
                     labels=c("Unacceptable", "Acceptable"))+
  ggtitle("Area distributions")+
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size=18),
        title = element_text(size=16), 
        legend.text = element_text(size=12))












