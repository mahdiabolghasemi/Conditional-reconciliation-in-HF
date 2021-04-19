library(hts)
library(plyr)
library(doSNOW)
library(matlib)
library(ggplot2)
library(randomForest)
library(rpart)
library(kernlab)
library(gbm)
library(seer)
library(tsfeatures)
library(forecast)
library(foreach)
library(dplyr)
library(xgboost)
library(tsutils)

opt_level <- "level_2"
setwd(paste0("C:/Users/vangelis spil/Desktop/Prison_CHF/",opt_level))

error_cal <- function(insample, outsample, forecasts){
  
  ppy=frequency(insample)
  masep<-mean(abs(diff(insample, lag = ppy)))
  
  outsample <- as.numeric(outsample)
  forecasts <- as.numeric(forecasts)
  
  mase <- mean(abs(outsample-forecasts))/masep
  amse <- abs(mean(outsample-forecasts))/masep
  rmsse <- sqrt(mean(abs(outsample-forecasts)^2)/(masep^2))
  
  return(mase)
}

Data <- read.csv("prison_population.csv", stringsAsFactors =  F)
Data$Legal = Data$Indigenous <- NULL
Data <- ddply(Data, .(Date,State,Gender), colwise(sum))
Data$StateOR <- Data$State ; Data$GenderOR <- Data$Gender
map <- c("A","B","C","D","E","F","G","H")
map_r <- unique(Data$State)
for (i in 1:length(map_r)){
  Data[Data$State==map_r[i],]$State <- map[i]
}
map2 <- c("0","1")
map2_r <- unique(Data$Gender)
for (i in 1:length(map2_r)){
  Data[Data$Gender==map2_r[i],]$Gender <- map2[i]
}
Data$Flag <- paste0(Data$State,Data$Gender)
series_ids <- unique(Data$Flag)

DataFor <- Data[Data$Flag==series_ids[1],c("Date","Count")]
colnames(DataFor)[2] <- series_ids[1] ; row.names(DataFor) <- NULL
for (i in 2:length(series_ids)){
  DataFor_t <- Data[Data$Flag==series_ids[i],c("Date","Count")] 
  colnames(DataFor_t)[2] <- series_ids[i] ; row.names(DataFor_t) <- NULL
  DataFor <- cbind(DataFor, DataFor_t)
  DataFor <- DataFor[,-(ncol(DataFor)-1)]
}
DataFor$Date <- NULL
for (i in 1:length(series_ids)){
  DataFor[,i] <- ts(DataFor[,i], frequency = 4, start = c(2005,1))
}
rm(map,map2,map2_r,map_r,i,DataFor_t,series_ids)

###################################################################################################
#Tourism data set
###################################################################################################
DataForHierarchy <- hts(ts(DataFor, frequency = 4, start = c(2005,1)), characters = c(1,1))
#plot(DataForHierarchy)

#######################################################################
#Reconcile 
########################################################################

#####first phase ######
data <- DataForHierarchy
m_series <- nrow(smatrix(data)) #Total number of series at all levels
b_series <- length(data$labels[[(length(data$labels)-1)+1]]) #Number of series at the bottom level

#Get base forecasts for training and forecast purposes
Forecast_file = Residuals_file <- NULL
Forecast_file_train <- NULL
Summary_error_train <- NULL
chr_train_stack <- ML_data_train <- NULL
Selected_model <- NULL
chr_tarin_all <- NULL

counter <- 0
origin <- 24 #Starting origin (6 years)
tslength <- 48 #Available observations for creating the ML train set
fh = 1 #Forecasting horizon considered

while ((origin+counter+fh) <= tslength){
  
  #Fetch data
  ts_sample <- head(aggts(data), origin+counter+fh)
  train_sample <- head(ts_sample, length(ts_sample[,1])-fh)
  test_sample <- tail(ts_sample, fh+1) #! Bug when having a single observation
  
  hts.data.scan.train <- hts(train_sample[,(ncol(train_sample)-b_series+1):ncol(train_sample)], characters = c(1,1))
  hts.data.scan.test <- hts(test_sample[,(ncol(test_sample)-b_series+1):ncol(test_sample), drop=FALSE], characters = c(1,1))
  
  fore.mint.shr.ro <- forecast(hts.data.scan.train, h=fh, method = "comb", 
                                fmethod = "ets", weights = "mint", covariance = "shr")
  
  fore.bu.ro <- forecast(hts.data.scan.train, h=fh, method = "bu", fmethod = "ets")
  
  
  fore.tdgsf.ro <- forecast(hts.data.scan.train, h=fh, method = "tdgsf", fmethod = "ets")
  
  ##Time series features ####
  
  ##### Forecasts accuracy #####
  test_sample_c <- tail(as.data.frame(allts(hts.data.scan.test)),fh)
  colnames(test_sample_c) <- paste0("A_", as.character(colnames(test_sample_c)))
  ### All time series forecasts #### 
  Min_Error <- 999999.99
  chrtrain <- t(cbind((tsfeatures(allts((hts.data.scan.train)), 
                                  features = c("fluctanal_prop_r1","frequency", "entropy",
                                               "acf_features","lumpiness","stability",
                                               "max_level_shift","max_var_shift","max_kl_shift",
                                               
                                               "hurst", "unitroot_kpss", "heterogeneity",
                                               "nonlinearity", "pacf_features"))), 
                      (tsfeatures(allts((hts.data.scan.train)), 
                                  "stl_features", s.window="periodic", robust= TRUE))) %>%
                  select(unitroot_kpss, trend, linearity, nonlinearity,
                         curvature, hurst, lumpiness, stability,
                         fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
                         max_level_shift, max_var_shift, max_kl_shift,
                         entropy, arch_acf, garch_acf, arch_r2, garch_r2,
                         diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10))
  
  chr.train <- cbind(t(chrtrain[,1]),t(chrtrain[,2]),t(chrtrain[,3]),t(rowMeans(chrtrain[,4:15])))
  
  
  repML=counter+1
  originML=length(train_sample[,1])
  
  chr_train_stack <- rbind(cbind((chr.train),repML,originML))
  
  chr_train_all <- rbind(chr_tarin_all, chr_train_stack)
  
  for (model_id in 1:3){
    
    if (model_id==1){
      allf <- data.frame(allts(fore.mint.shr.ro))
      
    }else if (model_id==2){
      allf <- data.frame(allts(fore.bu.ro))
      
    }else{
      allf <- data.frame(allts(fore.tdgsf.ro))
    }
    
    colnames(allf) <- paste0("F_", as.character(colnames(allts(hts.data.scan.test))))
    
    #Combine with actual values
    all <- cbind(allf, test_sample_c) 
    all$model <- model_id
    all$rep <- counter +1
    all$fh <- c(1:fh)
    all$period <- c((length(train_sample[,1])+1) : (length(train_sample[,1])+fh))
    all$origin <- length(train_sample[,1])
    
    Forecast_file_train <- rbind(Forecast_file_train, all)
    
    error_list_train <- NULL
    for (j in 1:ncol(allts(hts.data.scan.test))){
      error_list_train <- rbind(error_list_train,error_cal(allts(hts.data.scan.train)[,j],all[,j],tail(allts(hts.data.scan.test)[,j],fh)))
    }
    
    Errors <- data.frame(error_list_train) ; colnames(Errors) <- "MASE"
    
    Errors$model <- model_id
    Errors$origin <- origin
    Errors$hs <- c(1:j)
    Summary_error_train <- rbind(Summary_error_train, Errors)
    
    ############ Optimization criterion ################
    if (opt_level=="all_levels"){
      Min_Error_MASE <- (Errors[1,1]+mean(Errors[2:9,1])+mean(Errors[10:25,1]))/3 
    }else if (opt_level=="all_series"){
      Min_Error_MASE <- mean(Errors[,1]) 
    }else if (opt_level=="level_0"){
      Min_Error_MASE <- Errors[1,1]
    }else if (opt_level=="level_1"){
      Min_Error_MASE <- mean(Errors[2:9,1])
    }else if (opt_level=="level_2"){
      Min_Error_MASE <- mean(Errors[10:25,1])
    }
    ##################################################
    
    if (Min_Error_MASE < Min_Error){
      Min_Error <- Min_Error_MASE 
      Selected_model <- model_id} else {
        Min_Error <- Min_Error
        Selected_model <- Selected_model
      }
  }
  ML_data <- data.frame((cbind(t(data.frame(t(chr.train))), Selected_model, repML, originML)))
  ML_data_train <- rbind(ML_data_train, ML_data)
  
  counter <- counter + 1
} 

res <- ddply(Summary_error_train,.(model,hs), colwise(mean))
res$level <- 2
res[res$hs %in% c(2:9),]$level <- 1
res[res$hs %in% 1,]$level <- 0
ddply(res,.(level,model), colwise(mean)) #Error per level
ddply(ddply(res,.(level,model), colwise(mean)),.(model), colwise(mean)) #Average error

######## second phase ########
models_selected <- c( "Shr", "BU","TDGSF")
#This is what I chose for the good XGb results
#ML_data_train.dataframe <- as.data.frame(ML_data_train)[,c(1:25,76:104)]
ML_data_train.dataframe <- as.data.frame(cbind(ML_data_train[,1:25],(ML_data_train[,26:50] +ML_data_train[,51:75])/2, ML_data_train[,76:100], ML_data_train[,101:103]))
table(ML_data_train.dataframe$Selected_model) #Flags per method
##I excluded some other features belwo
write.csv(ML_data_train.dataframe,"MASE_ML_data_train.dataframe.csv")
write.csv(ML_data_train,"MASE_ML_data_train.csv")
write.csv(Forecast_file_train,"MASE_Forecast_file_train.csv")
write.csv(Summary_error_train,"MASE_Summary_error_train.csv")
save.image("Prison_train.Rdata")

load("Prison_train.Rdata")
Forecast_file_test <- NULL
chr_test_stack <- NULL
Errors_test1 <- NULL
error_list_test1 <- NULL
Summary_error_test1 <- NULL
Summary_error_test_f <- NULL
recon_bottom_f_all <- NULL
chr_test <- NULL 
recon_f <- NULL
recon_bottom_f <- NULL
selected_vars_all <- NULL
selected_vars_all_name <- NULL

counter <- 0
origin_id <- 40 #Starting origin (26weeks)
tslength_id <- 48 #Available observations
fh = 1 #Forecasting horizon considered
svmC = xgvC = RFc <- c()

while ((origin_id+counter+fh) <= tslength_id){
  
  #Fetch data
  ts_sample <- head(aggts(data), origin_id+counter+fh)
  train_sample <- head(ts_sample, length(ts_sample[,1])-fh)
  test_sample <- tail(ts_sample, fh+1) #! Bug when having a single observation
  
  hts.data.scan.train <- hts(train_sample[,(ncol(train_sample)-b_series+1):ncol(train_sample)], characters = c(1,1))
  hts.data.scan.test <- hts(test_sample[,(ncol(test_sample)-b_series+1):ncol(test_sample), drop=FALSE], characters = c(1,1))
  
  fore.mint.shr.ro <- forecast(hts.data.scan.train, h=fh, method = "comb", 
                               fmethod = "ets", weights = "mint", covariance = "shr")
  
  fore.bu.ro <- forecast(hts.data.scan.train, h=fh, method = "bu", fmethod = "ets")
  
  
  fore.tdgsf.ro <- forecast(hts.data.scan.train, h=fh, method = "tdgsf", fmethod = "ets")
  
  Errors_test1 <- NULL
  Errors_test <- NULL
  
  error_list_test1 <- NULL
  
  for (model_id in 1:3){
    
    if (model_id==1){
      all_test_f <- data.frame(allts(fore.mint.shr.ro))
      
    }else if (model_id==2){
      all_test_f <- data.frame(allts(fore.bu.ro))
      
    }else{
      all_test_f <- data.frame(allts(fore.tdgsf.ro))
    }
    
    error_list_test1 <- NULL
    for (j in 1:ncol(allts(hts.data.scan.test))){
      error_list_test1 <- rbind(error_list_test1,error_cal(allts(hts.data.scan.train)[,j],(all_test_f)[,j],tail(allts(hts.data.scan.test)[,j],fh)))
    }
    Errors_test1 <- cbind(Errors_test1,error_list_test1)
  }
  
  Errors_test1 <- as.data.frame(Errors_test1)  
  colnames(Errors_test1) <- models_selected
  
  ############ Optimization criterion ################
  if (opt_level=="all_levels"){
    Errors_test1$best <- which.min(c((Errors_test1[1,1]+mean(Errors_test1[2:9,1])+mean(Errors_test1[10:25,1]))/3,
                                     (Errors_test1[1,2]+mean(Errors_test1[2:9,2])+mean(Errors_test1[10:25,2]))/3,
                                     (Errors_test1[1,3]+mean(Errors_test1[2:9,3])+mean(Errors_test1[10:25,3]))/3))
  }else if (opt_level=="all_series"){
    Errors_test1$best <- which.min(colMeans(Errors_test1[,1:3]))
  }else if (opt_level=="level_0"){
    Errors_test1$best <- which.min(c(Errors_test1[1,1],
                                     Errors_test1[1,2],
                                     Errors_test1[1,3]))
  }else if (opt_level=="level_1"){
    Errors_test1$best <- which.min(c(mean(Errors_test1[2:9,1]),
                                     mean(Errors_test1[2:9,2]),
                                     mean(Errors_test1[2:9,3])))
  }else if (opt_level=="level_2"){
    Errors_test1$best <- which.min(c(mean(Errors_test1[10:25,1]),
                                     mean(Errors_test1[10:25,2]),
                                     mean(Errors_test1[10:25,3])))
  }
  ##################################################
  
  ######################################################################
  Errors_test1$counter <- counter+2
  Errors_test1$level <- c(1,rep(2,8),rep(3,16))
  Summary_error_test1 <- rbind(Summary_error_test1, Errors_test1) 
  
  ## forecasting the bottom series with ML ###
  chr_test <- chrtest <- NULL
  for (hs in 1:m_series) {
    chr.test <- cbind(t(cbind((tsfeatures(allts((hts.data.scan.train))[,hs], 
                                          features = c("fluctanal_prop_r1","frequency", "entropy",
                                                       "acf_features","lumpiness","stability",
                                                       "max_level_shift","max_var_shift","max_kl_shift",
                                                       
                                                       "hurst", "unitroot_kpss", "heterogeneity",
                                                       "nonlinearity", "pacf_features"))), 
                              (tsfeatures(allts((hts.data.scan.train))[,hs], 
                                          "stl_features", s.window="periodic", robust= TRUE)))%>%
                          select( unitroot_kpss,
                                  trend,linearity, nonlinearity,curvature,
                                  hurst, lumpiness, stability,
                                  fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
                                  max_level_shift, max_var_shift, max_kl_shift,entropy,
                                  arch_acf,garch_acf, arch_r2, garch_r2,
                                  diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)))
    
    
    chrtest <- rbind(chrtest,(t(chr.test)))
  }
  #This one used features of series at all levels
  chr_test <- data.frame(cbind(t(chrtest[1,]),(t(chrtest[2,])+t(chrtest[3,])/2),t(colMeans(chrtest[4:15,])),t(colMeans(chrtest[4:15,]))))
  chr_test <- chr_test[,-c(51:75)]
  
  # This one for just bottom level
  repML=counter+1
  originML=length(train_sample[,1])
  
  chr_test_stack <- rbind(chr_test_stack,cbind(chr_test,repML,originML))
  
  
  ### ML model to classify the appropriate reconciliation method ####
  Reconcile_data <-  ML_data_train.dataframe  
  #Reconcile_data <- Reconcile_data[Reconcile_data$originML<origin_id,]
  Reconcile_data <- Reconcile_data[Reconcile_data$originML<length(allts(hts.data.scan.train)[,1]),]
  
  
  y_value <- t(as.factor(Reconcile_data$Selected_model))
  
  set.seed(100)
  ML_model <- randomForest(x = ((Reconcile_data)[,-c(dim(Reconcile_data)[2]:(dim(Reconcile_data)[2]-3))]) , y= y_value, ntree = 500)
  selected_vars <- sort((ML_model$importance[,1]), decreasing = TRUE)[1:25]
  selected_vars_name <- names(selected_vars)
  
  selected_vars_all <- rbind(selected_vars_all,selected_vars)
  selected_vars_all_name <- rbind(selected_vars_all_name,selected_vars_name)
  
  chr_test_selected <- cbind(chr_test[which(names(chr_test) %in% c(names(selected_vars)))],y=Errors_test1$best[1])
  
  set.seed(100)
  ML_model_RF <- randomForest(x = Reconcile_data[which(names(Reconcile_data) %in% c(names(selected_vars)))] ,
                              y= y_value, ntree = 250)
  
  
  RF_output <- predict(ML_model_RF, chr_test_selected[1:25])
  
  
  
  dat= as.data.frame(cbind(Reconcile_data[which(names(Reconcile_data) %in% c(names(selected_vars)))] ,
                           y=t(y_value)))
  set.seed(100)
  kernfit <- ksvm(as.matrix(dat[,-dim(dat)[2]]),dat$y, type = "C-svc", kernel = 'rbfdot', 
                  C = 100, scaled = TRUE)
  
  SVM_output <- predict(kernfit, chr_test_selected[1:25])
  
  #### XGBOOst
  train_x = data.matrix(dat[,-dim(dat)[2]])
  train_y = dat[,dim(dat)[2]]
  
  test_x = data.matrix(chr_test_selected[1:25])
  test_y = chr_test_selected[length(chr_test_selected)]
  
  train.y = (as.numeric(dat[,dim(dat)[2]])-1)
  xgbtrain = xgb.DMatrix(data=train_x, label=train.y)
  
  testy = as.numeric(chr_test_selected[length(chr_test_selected)])-1
  xgbtest = xgb.DMatrix(data=test_x, label=testy) 
  set.seed(100)
  XGB_class <- xgb.train(booster="gbtree",
                         eta=0.3,
                         max_depth=5,
                         gamma=3,
                         subsample=0.4,
                         colsample_bytree=1,
                         objective="multi:softprob",
                         eval_metric="mlogloss",
                         num_class=3, 
                         data=xgbtrain,
                         nrounds=1000,
                         nthreads=3,
                         early_stopping_rounds=10,
                         watchlist=list(val1=xgbtrain,val2=xgbtest),
                         verbose=0)
  
  xgb.pred = predict(XGB_class,xgbtest,reshape=T)
  xgb.pred = as.data.frame(xgb.pred);colnames(xgb.pred) <- c(1,2,3)
  xgb.pred$prediction = apply(xgb.pred,1,function(x) colnames(xgb.pred)[which.max(x)])
  
  svmC <- c(svmC, as.numeric(SVM_output)) 
  xgvC <- c(xgvC, as.numeric(xgb.pred$prediction))
  RFc <- c(RFc, as.numeric(RF_output))
  
  ML_output <- cbind(RF_output,SVM_output, xgb.pred$prediction)
  
  for (ML in 1:3) {
    
    if (ML_output[ML]==1){
      recon_f <- data.frame(allts(fore.mint.shr.ro))
      
    }else if (ML_output[ML]==2){
      recon_f <- data.frame(allts(fore.bu.ro))
      
    }else{
      recon_f <- data.frame(allts(fore.tdgsf.ro))
    }
    
    
    recon_bottom_f_all <- rbind(recon_bottom_f_all, cbind(recon_f, counter, ML, as.numeric(ML_output[ML]),   Errors_test1$best[1])) 
    ##Time series features ####
    
    ##### Forecasts accuracy #####
    test_sample <- tail(as.data.frame(allts(hts.data.scan.test)),fh)
    colnames(test_sample) <- paste0("A_", as.character(colnames(test_sample)))
    ### All time series forecasts #### 
    error_list_test <- NULL
    
    for (j in 1:ncol(allts(hts.data.scan.test))){
      error_list_test <- rbind(error_list_test,error_cal(allts(hts.data.scan.train)[,j],recon_f[,j],tail(allts(hts.data.scan.test)[,j],fh)))
    }
    Errors_test <- data.frame(cbind(Errors_test,(error_list_test))) 
  }
  colnames(Errors_test) <- c("CR_RF","CR_SVM", "XGB_class")
  
  Summary_error_test_f <- rbind(Summary_error_test_f, cbind(Errors_test,Errors_test1))
  counter <- counter + 1
}

correct_class <- which(recon_bottom_f_all$`as.numeric(ML_output[ML])`== recon_bottom_f_all$`Errors_test1$best[1]`)

write.csv(correct_class,"MASE_correct_class_test.csv")
write.csv(Summary_error_test_f, "MASE_Summary_error_test.csv")
write.csv(recon_bottom_f_all,"MASE_recon_bottom.csv")
write.csv(chr_test_stack, "MASE_chr_test_stack.csv")
write.csv(chr_test_selected, "MASE_chr_test_selected.csv")
write.csv(selected_vars_all, "MASE_selected_vars_all.csv")

res <- ddply(Summary_error_test_f,.(level), colwise(mean)) #Error per level
res
colMeans(res) #Average error

dataMCB <- Summary_error_test_f
dataMCB$X = dataMCB$best <- NULL
colnames(dataMCB)[1] <- "CHF-RF"
colnames(dataMCB)[2] <- "CHF-SVM"
colnames(dataMCB)[3] <- "CHF-XGB"
if (opt_level=="all_levels"){
  dataMCB_3 <- dataMCB
  dataMCB_3 <- ddply(dataMCB_3, .(counter,level), colwise(mean))
  dataMCB_3$level <- NULL
  dataMCB_3 <- ddply(dataMCB_3, .(counter), colwise(mean))
  dataMCB_3$counter <- NULL
  nemenyi(as.matrix(dataMCB_3), conf.level = 0.95, plottype = "vmcb", main="Average per level")
}else if (opt_level=="all_series"){
  dataMCB_3 <- dataMCB
  dataMCB_3$counter = dataMCB_3$level <- NULL
  nemenyi(as.matrix(dataMCB_3), conf.level = 0.95, plottype = "vmcb", main="All series")
}else if (opt_level=="level_0"){
  dataMCB_3 <- dataMCB[dataMCB$level==1,]
  dataMCB_3$counter = dataMCB_3$level <- NULL
  nemenyi(as.matrix(dataMCB_3), conf.level = 0.95, plottype = "vmcb", main="Level 0")
}else if (opt_level=="level_1"){
  dataMCB_3 <- dataMCB[dataMCB$level==2,]
  dataMCB_3$counter = dataMCB_3$level <- NULL
  nemenyi(as.matrix(dataMCB_3), conf.level = 0.95, plottype = "vmcb", main="Level 1")
}else if (opt_level=="level_2"){
  dataMCB_3 <- dataMCB[dataMCB$level==3,]
  dataMCB_3$counter = dataMCB_3$level <- NULL
  nemenyi(as.matrix(dataMCB_3), conf.level = 0.95, plottype = "vmcb", main="Level 2")
}

k <- data.frame(RFc,svmC,xgvC,
                unique(Summary_error_test_f[,c("counter","best")])$best)
colnames(k)[4] <- "actual"
nrow(k[k$RFc==k$actual,])/nrow(k)
nrow(k[k$svmC==k$actual,])/nrow(k)
nrow(k[k$xgvC==k$actual,])/nrow(k)
k


