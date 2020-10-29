library(hts)
library(plyr)
library(foreach)
library(doSNOW)
library(matlib)
library(ggplot2)
library(randomForest)
library(rpart)
library(kernlab)
library(gbm)
library(seer)
library(tsfeatures)
library(caret)
library(mcb)
library(randomForest)
library(xgboost)
library(tstools)

setwd("C:/Temp/Model selection")

error_cal <- function(insample, outsample, forecasts, ppy){
  
  #insample=Insample[,i];outsample=Forecasts_temp[,i];forecasts=Outsample[,i];ppy=12
  masep<-mean(abs(diff(insample, lag = ppy)))
  
  outsample <- as.numeric(outsample)
  forecasts <- as.numeric(forecasts)
  
  #mase <- mean(abs(outsample-forecasts))/masep
  #amse <- abs(mean(outsample-forecasts))/masep
  rmsse <- sqrt(mean(abs(outsample-forecasts)^2)/(masep^2))
  
  #output <- c(mase,amse,rmsse) ; names(output) <- c("MASE","AMSE","RMSSE")
  
  return(rmsse)
}



###################################################################################################
#Read data & Create validation set
###################################################################################################
#Read file (from 01/1998 to 12/2017)
input <- read.csv("datanew.csv", stringsAsFactors = F)


library(janitor)
data = janitor::clean_names(input)

library(tidyverse)
data = data %>% dplyr::rename(
  date = week_commencing,
  category = category1,
  brand = brand1,
  flavour = flavour1,
  size = unit_size1,
  cbfs = name1,
  scan_units = sum_of_measures_actual_scan_cases,
  scan_cases = sales,
  gross_sales = sum_of_measures_actual_scan_retailer_gross_sales,
  prices = actual_scan_sales_price_per_unit,
  sku = namels,
  sl = ls#,
) %>% 
  dplyr::mutate(date = lubridate::dmy(date))

# group by customer state
sl.names <- data %>% group_by(sl) %>% 
  nest()


# group_by both
hr.data <- data %>% 
  # subset so that all obercations have a date greater than some value
  filter(date >= lubridate::dmy("28-08-2016")) %>% 
  # group by to establish the categories
  group_by(sl,cbfs) %>% 
  # subset so that you have a minimum number of observations in each category
  filter(n()>112) %>% 
  # filter to only those products that have 12 time series
  group_by(cbfs) %>% 
  filter(n_distinct(sl) == 12)  %>% nest()


hr.data.n <- hr.data.n.all <- list()
for (i in 1:length(hr.data$data)) {
  for (j in 1:12) {
    hr.data.n[[j]] <- hr.data$data[[i]] %>% filter(sl==sl.names$sl[j])
  }
  hr.data.n.all[[i]] <- hr.data.n 
}

## extracting scan for 12 sl(states customers) and 62 SKU
scan = scan1 = list()
for (i in 1:length(hr.data.n.all)) {
  for (j in 1:12) {
    scan1[[j]] <- ts(hr.data.n.all[[i]][[j]]$scan_cases[10:129], frequency = 52, start = c(2016,1))
  }
  scan[[i]] <- scan1
}

## extracting prices for 12 sl and 62 SKU
price = price1 = list()
for (i in 1:length(hr.data.n.all)) {
  for (j in 1:12) {
    price1[[j]] <- hr.data.n.all[[i]][[j]]$prices[10:129]
  }
  price[[i]] <- price1
}

outliers <- c(49,50,54, 62,11,37,38)
scan=scan[-outliers]
price=price[-outliers]

#### hierarchical time series model for total Scan ####
hts.scan <- list()
nodes <- list(2, c(6, 6))
bnames = c("ANSW","AQLD","ASAT",
           "ATAS","AVIC","AWAT",   
           "BNSW","BQLD","BSAT",
           "BTAS","BVIC","BWAT")

scan.matrix <- matrix(NA, ncol = length(bnames), nrow = 120)

for (i in 1:length(scan)) {
  scan.matrix <- bind_cols(scan[[i]])
  colnames(scan.matrix) <- bnames
  hts.scan[[i]] <- hts(log(scan.matrix), nodes=nodes)
}

#### hierarchical series of price
#### hierarchical time series model for total Scan ####
hts.price <- list()
nodes <- list(2, c(6, 6))
bnames = c("ANSW","AQLD","ASAT",
           "ATAS","AVIC","AWAT",   
           "BNSW","BQLD","BSAT",
           "BTAS","BVIC","BWAT")

price.matrix <- matrix(NA, ncol = length(bnames), nrow = 120)

for (i in 1:length(price)) {
  price.matrix <- bind_cols(price[[i]])
  colnames(price.matrix) <- bnames
  hts.price[[i]] <- hts(price.matrix, nodes=nodes)
}



#######################################################################
#Reconcile 
########################################################################

#####first phase ######
data <- hts.scan
data.price <- hts.price

structure.n <-  window(data, start = 1)[[1]]
m_series <- nrow(smatrix(structure.n)) #Total number of series at all levels
b_series <- length(structure.n$labels[[(length(structure.n$labels)-1)+1]]) #Number of series at the bottom level

#Get base forecasts for training and forecast purposes
Forecast_file = Residuals_file <- NULL
Forecast_file_train <- NULL
Summary_error_train <- NULL
chr_train_stack <- ML_data_train <- NULL
Selected_model <- NULL
chr_tarin_all <- NULL

for (tsid in 1:55) {
  
  counter <- 0
  origin <- 26 #Starting origin (26weeks)
  tslength <- 85 #Available observations
  #tslength <- length(data[[tsid]]$bts[,1]) #Available observations
  fh = 2 #Forecasting window considered
  
  while ((origin+counter+fh) <= tslength){
    
    #Fetch data
    ts_sample <- head(aggts(data[[tsid]]), origin+counter+fh)
    train_sample <- head(ts_sample, length(ts_sample[,1])-fh)
    test_sample <- tail(ts_sample, fh)
    
    ## Include price
    ts_sample_p <- head(aggts(data.price[[tsid]]), origin+counter+fh)
    train_sample_p <- head(ts_sample_p, length(ts_sample_p[,1])-fh)
    test_sample_p <- tail(ts_sample_p, fh)
    
    
    hts.data.scan.train <- hts(((train_sample)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
    hts.data.scan.test <- hts(((test_sample)[,4:15, drop=FALSE]), nodes=list(2, c(6,6)), bnames = bnames)
    
    hts.data.price.train <- hts(((train_sample_p)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
    hts.data.price.test <- hts(((test_sample_p)[,4:15, drop=FALSE]), nodes=list(2, c(6,6)), bnames = bnames)
    
    
    fore.mint.shr.ro <- (forecast(hts.data.scan.train, h=fh, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                                  , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)), weights = c("mint"), covariance = "shr"))
    
    fore.bu.ro <- forecast(hts.data.scan.train, h=fh, method = "bu", fmethod = "arima", xreg = as.numeric(aggts(hts.data.price.train, levels = 0))
                           , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
   
    fore.tdgsf.ro <- forecast(hts.data.scan.train, h=fh, method = "tdgsf", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                              , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
    ##Time series features ####
    
    ##### Forecasts accuracy #####
    
    #models <- c("WLS", "Shr", "Sam", "BU", "MO", "TDFP", "TDGSA", "TDGSF")
    test_sample_c <- as.data.frame(allts(hts.data.scan.test))
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
                    #  (tsfeatures(allts((hts.data.scan.train))[,hs],
                    #             c("mean","var"), scale=FALSE, na.rm=TRUE)))%>%
                    select(unitroot_kpss, trend, linearity, nonlinearity,
                           curvature, hurst, lumpiness, stability,
                           fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
                           max_level_shift, max_var_shift, max_kl_shift,
                           entropy, arch_acf, garch_acf, arch_r2, garch_r2,
                           diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10))
    
    
    #train_cor <- rcorr(as.matrix(allts(hts.data.scan.train)))$r
    #train_cor <- cbind(min(train_cor[1,2],train_cor[1,3]),min(train_cor[1,4:15]), min(train_cor[2,4:15],train_cor[3,4:15]),
     #                  min(train_cor[4:15,4:15]), train_cor[2,3]) 
    #colnames(train_cor) <- c("cor01","cor02","cor12", "corL2", "corL1")
    
    
    chr.train <- cbind(t(chrtrain[,1]),t(chrtrain[,2]),t(chrtrain[,3]),t(rowMeans(chrtrain[,4:15])))
    
    repML=counter+1
    originML=length(train_sample[,1])
    
    chr_train_stack <- rbind(cbind((chr.train),repML,originML,tsid))
    
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
      all$tsid <- tsid
      all$rep <- counter +1
      all$fh <- c(1:fh)
      all$period <- c((length(train_sample[,1])+1) : (length(train_sample[,1])+fh))
      all$origin <- length(train_sample[,1])
      
      Forecast_file_train <- rbind(Forecast_file_train, all)
      
      
      
      error_list_train <- NULL
      for (j in 1:ncol(allts(hts.data.scan.test))){
        error_list_train <- rbind(error_list_train,error_cal(allts(hts.data.scan.train)[,j],all[,j],allts(hts.data.scan.test)[,j],ppy=1))
      }
      
      Errors <- data.frame(error_list_train) ; colnames(Errors) <- ("RMSSE")#,"AMSE","RMSSE")
     
      Errors$model <- model_id
      Errors$tsid <- tsid
      Errors$origin <- origin
      Errors$hs <- c(1:j)
      #Errors$best <- which.min(colMeans(Errors[,1:3]))
      Summary_error_train <- rbind(Summary_error_train, Errors)
      
      
      #Min_Error_MASE <- (Errors[hs,1])
      
      Min_Error_RMSSE <- mean(Errors[,1])
      
      #Min_Error <- Min_Error_MASE
      if (Min_Error_RMSSE < Min_Error){
        Min_Error <- Min_Error_RMSSE 
        Selected_model <- model_id} else {
          Min_Error <- Min_Error
          Selected_model <- Selected_model
        }
    }
    ML_data <- data.frame((cbind(t(data.frame(t(chr.train))), Selected_model, repML, originML, tsid)))
    ML_data_train <- rbind(ML_data_train, ML_data)
  
    counter <- counter + 2
  } 
}



#models <- c("WLS", "Shr", "BU", "MO", "TDFP", "TDGSA", "TDGSF")
models_selected <- c( "Shr", "BU","TDGSF")



write.csv(ML_data_train.dataframe,"RMSSE_ML_data_train.dataframe.csv")
write.csv(ML_data_train,"RMSSE_ML_data_train.csv")
write.csv(Forecast_file_train,"RMSSE_Forecast_file_train.csv")
write.csv(Summary_error_train,"RMSSE_Summary_error_train.csv")


######## second phase ########
#This is what I chose for the good XGb results
#ML_data_train.dataframe <- as.data.frame(ML_data_train)[,c(1:25,76:104)]
ML_data_train.dataframe <- as.data.frame(cbind(ML_data_train[,1:25],(ML_data_train[,26:50] +ML_data_train[,51:75])/2, ML_data_train[,76:100], ML_data_train[,101:104]))


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
#fh = 5
for (tsid in 1:55) {
  
  counter <- 0
  origin_id <- 83 #Starting origin (26weeks)
  tslength_id <- 120 #Available observations
  #tslength <- length(data[[tsid]]$bts[,1]) #Available observations
  fh = 2 #Forecasting horizon considered
  
  while ((origin_id+counter+fh) <= tslength_id){
    
    #Fetch data
    ts_sample <- head(aggts(data[[tsid]]), origin_id+counter+fh)
    train_sample <- head(ts_sample, length(ts_sample[,1])-fh)
    test_sample <- tail(ts_sample, fh)
    
    ## Include price
    ts_sample_p <- head(aggts(data.price[[tsid]]), origin_id+counter+fh)
    train_sample_p <- head(ts_sample_p, length(ts_sample_p[,1])-fh)
    test_sample_p <- tail(ts_sample_p, fh)
    
    
    hts.data.scan.train <- hts(((train_sample)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
    hts.data.scan.test <- hts(((test_sample)[,4:15,drop=FALSE]), nodes=list(2, c(6,6)), bnames = bnames)
    
    hts.data.price.train <- hts(((train_sample_p)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
    hts.data.price.test <- hts(((test_sample_p)[,4:15,drop=FALSE]), nodes=list(2, c(6,6)), bnames = bnames)
    
    
    
    fore.mint.shr.ro <- (forecast(hts.data.scan.train, h=fh, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                                  , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)), weights = c("mint"), covariance = "shr"))
    
     
    fore.bu.ro <- forecast(hts.data.scan.train, h=fh, method = "bu", fmethod = "arima", xreg = as.numeric(aggts(hts.data.price.train, levels = 0))
                           , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
   
    fore.tdgsf.ro <- forecast(hts.data.scan.train, h=fh, method = "tdgsf", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                              , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
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
        error_list_test1 <- rbind(error_list_test1,error_cal(allts(hts.data.scan.train)[,j],(all_test_f)[,j],allts(hts.data.scan.test)[,j],ppy=1))
        #error_list_train <- rbind(error_list_train,error_cal(allts(hts.data.scan.train)[,hs],all[,hs],allts(hts.data.scan.test)[,hs],ppy=1))
      }
      Errors_test1 <- cbind(Errors_test1,error_list_test1)
    }

    Errors_test1 <- as.data.frame(Errors_test1)  
    colnames(Errors_test1) <- models_selected
   
    Errors_test1$tsid <- tsid
    Errors_test1$best <- which.min(colMeans(Errors_test1[,1:3]))
    Errors_test1$counter <- counter+2
    Errors_test1$level <- c(1,2,2,3,3,3,3,3,3,3,3,3,3,3,3)
    Summary_error_test1 <- rbind(Summary_error_test1, Errors_test1) 
    
    ## forecasting the bottom series with ML ###
    chr_test <- chrtest <- NULL
    for (hs in 1:m_series) {
      #Min_Error <- 999999.99
      chr.test <- cbind(t(cbind((tsfeatures(allts((hts.data.scan.train))[,hs], 
                                            features = c("fluctanal_prop_r1","frequency", "entropy",
                                                         "acf_features","lumpiness","stability",
                                                         "max_level_shift","max_var_shift","max_kl_shift",
                                                         
                                                         "hurst", "unitroot_kpss", "heterogeneity",
                                                         "nonlinearity", "pacf_features"))), 
                                (tsfeatures(allts((hts.data.scan.train))[,hs], 
                                            "stl_features", s.window="periodic", robust= TRUE)))%>%
                            #(tsfeatures(allts((hts.data.scan.train))[,hs],
                            #           c("mean","var"), scale=FALSE, na.rm=TRUE)))%>%
                            select( unitroot_kpss,
                                    trend,linearity, nonlinearity,curvature,
                                    hurst, lumpiness, stability,
                                    fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
                                    max_level_shift, max_var_shift, max_kl_shift,entropy,
                                    arch_acf,garch_acf, arch_r2, garch_r2,
                                    diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)))
      
      
      chrtest <- rbind(chrtest,(t(chr.test)))
    }
   # test_cor <- rcorr(as.matrix(allts(hts.data.scan.test)))$r
    #test_cor <- cbind(min(test_cor[1,2],test_cor[1,3]),min(test_cor[1,4:15]), min(test_cor[2,4:15],test_cor[3,4:15]),
     #                  min(test_cor[4:15,4:15]), test_cor[2,3]) 
    #colnames(test_cor) <- c("cor01","cor02","cor12", "corL2", "corL1")
    
    #This one used features of series at all levels
    chr_test <- data.frame(cbind(t(chrtest[1,]),(t(chrtest[2,])+t(chrtest[3,])/2),t(colMeans(chrtest[4:15,])),t(colMeans(chrtest[4:15,]))))
    chr_test <- chr_test[,-c(51:75)]
    #chr_test <- data.frame(cbind(t(chrtest[1,]),t(colMeans(chrtest[4:15,]))))
    
    # This one for just bottom level
    #chr_test <- data.frame(t(colMeans(chrtest)))
    
    repML=counter+1
    originML=length(train_sample[,1])
    
    chr_test_stack <- rbind(chr_test_stack,cbind(chr_test,repML,originML,tsid))
    
   
    ### ML model to classify the appropriate reconciliation method ####
    Reconcile_data <-  ML_data_train.dataframe  
    
    
    y_value <- t(as.factor(Reconcile_data$Selected_model))
    y_value_reg <-  t(as.numeric(Reconcile_data$Selected_model))
    #ML_model <- randomForest(x = ((Reconcile_data)[,2:109]) , y= y_value, ntree = 500)
    # ML_output <- predict(ML_model, newdata = (chr_test))
    
    #y_value <- t(as.factor(Reconcile_data$Selected_model))
    
    ML_model <- randomForest(x = ((Reconcile_data)[,-c(dim(Reconcile_data)[2]:(dim(Reconcile_data)[2]-3))]) , y= y_value, ntree = 500)
    selected_vars <- sort((ML_model$importance[,1]), decreasing = TRUE)[1:25]
    selected_vars_all <- rbind(selected_vars_all,selected_vars)
    chr_test_selected <- cbind(chr_test[which(names(chr_test) %in% c(names(selected_vars)))],y=Errors_test1$best[1])
    
    ML_model_RF <- randomForest(x = Reconcile_data[which(names(Reconcile_data) %in% c(names(selected_vars)))] ,
                                y= y_value, ntree = 150)
    
    
    RF_output <- predict(ML_model_RF, chr_test_selected[1:25])
   
    
    
    dat= as.data.frame(cbind(Reconcile_data[which(names(Reconcile_data) %in% c(names(selected_vars)))] ,
                             y=t(y_value)))
    
    
    kernfit <- ksvm(as.matrix(dat[,-dim(dat)[2]]),dat$y, type = "kbb-svc", kernel = 'rbfdot', 
                    C = 200, scaled = TRUE)
    
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
    
    XGB_class <- xgb.train(booster="gbtree",
                           eta=0.01,
                           max_depth=5,
                           gamma=3,
                           subsample=0.75,
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
    
    ML_output <- cbind(RF_output,SVM_output, xgb.pred$prediction)
    
    for (ML in 1:3) {
      if (ML_output[ML]==1){
        recon_f <- data.frame(allts(fore.mint.shr.ro))#[,hs]
        
      }else if (ML_output[ML]==2){
        recon_f <- data.frame(allts(fore.bu.ro))#[,hs]
        
      }else{
        recon_f <- data.frame(allts(fore.tdgsf.ro))#[,hs]
      }
     
      
      hts.bu.forecast <- hts((recon_f[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
      
      recon_bottom_f_all <- rbind(recon_bottom_f_all, cbind(recon_f, tsid, counter, ML, as.numeric(ML_output[ML]),   Errors_test1$best[1])) 
      ##Time series features ####
      
      ##### Forecasts accuracy #####
      
      #models <- c("WLS", "Shr", "Sam", "BU", "MO", "TDFP", "TDGSA", "TDGSF")
      test_sample <- as.data.frame(allts(hts.data.scan.test))
      colnames(test_sample) <- paste0("A_", as.character(colnames(test_sample)))
      ### All time series forecasts #### 
      error_list_test <- NULL
  
      for (j in 1:ncol(allts(hts.data.scan.test))){
        error_list_test <- rbind(error_list_test,error_cal(allts(hts.data.scan.train)[,j],allts(hts.bu.forecast)[,j],allts(hts.data.scan.test)[,j],ppy=1))
      }
      Errors_test <- data.frame(cbind(Errors_test,(error_list_test))) 
    }
    colnames(Errors_test) <- c("CR_RF","CR_SVM", "XGB_class")
    
    Summary_error_test_f <- rbind(Summary_error_test_f, cbind(Errors_test,Errors_test1))
    counter <- counter + 2
  }
  
}

correct_class <- which(recon_bottom_f_all$`as.numeric(ML_output[ML])`== recon_bottom_f_all$`Errors_test1$best[1]`)
write.csv(correct_class,"RMSSE_correct_class.csv")

write.csv(Summary_error_test_f, "RMSSE_Summary_error.csv")


write.csv(recon_bottom_f_all,"RMSSE_recon_bottom.csv")
#recon_bottom_f_all <- read.csv("EMSSE_recon_bottom.csv")

write.csv(chr_test_stack, "RMSSE_chr_test_stack.csv")
#chr_test_stack <- read.csv("RMSSE_chr_test_stack.csv")

write.csv(chr_test_selected, "RMSSE_chr_test_selected.csv")

