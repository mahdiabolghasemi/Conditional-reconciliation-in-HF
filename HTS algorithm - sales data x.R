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

setwd("C:/Temp")

error_cal <- function(insample, outsample, forecasts, ppy){
  
  #insample=Insample[,i];outsample=Forecasts_temp[,i];forecasts=Outsample[,i];ppy=12
  masep<-mean(abs(diff(insample, lag = ppy)))
  
  outsample <- as.numeric(outsample)
  forecasts <- as.numeric(forecasts)
  
  mase <- mean(abs(outsample-forecasts))/masep
  #amse <- abs(mean(outsample-forecasts))/masep
  #rmsse <- sqrt(mean(abs(outsample-forecasts)^2)/(masep^2))
  
  #output <- c(mase,amse,rmsse) ; names(output) <- c("MASE","AMSE","RMSSE")
 
  return(mase)
}
reconcile <- function(frc_base){
  
  S <- smatrix(structure.n)
  W <- diag(rowSums(S))
  
  fb <- t(as.matrix(as.data.frame(frc_base)))
  
  frcst <- S%*%matlib::inv(t(S)%*% matlib::inv(W)%*%S)%*%t(S)%*% matlib::inv(W)%*%fb
  return(frcst)
}
bu <- function(frc_base){
  
  S <- smatrix(structure.n)
  m <- nrow(S) #total number of series
  K <- length(structure.n$labels)-1 #Numer of levels in the hierarchy
  mK <- length(structure.n$labels[[K+1]]) #Numer of series at the bottom level of the hierarchy
  P <- cbind(matrix(0, nrow = mK,  ncol = m-mK) ,   diag(1, mK))
  
  fb <- t(as.matrix(as.data.frame(frc_base)))
  frcst <- S%*%P%*%fb
  return(frcst)
}
td <- function(frc_base){
  
  S <- smatrix(structure.n)
  m <- nrow(S) #total number of series
  K <- length(structure.n$labels)-1 #Numer of levels in the hierarchy
  mK <- length(structure.n$labels[[K+1]]) #Numer of series at the bottom level of the hierarchy
  
  weigtsTD <- colSums(allts(structure.n))/colSums(allts(structure.n))[1]
  weigtsTD <- matrix(tail(weigtsTD, mK))
  
  P <- cbind( weigtsTD,  matrix(0, nrow = mK,  ncol = m-1))
  
  fb <- t(as.matrix(as.data.frame(frc_base)))
  frcst <- S%*%P%*%fb
  return(frcst)
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
  #cnc = cnc,
  #fff=fff
) %>% 
  dplyr::mutate(date = lubridate::dmy(date))

# group by customer state
sl.names <- data %>% group_by(sl) %>% 
  nest()


# group_by both
hr.data <- data %>% 
  # subset so that all obercations have a date greater than some value
  filter(date >= lubridate::dmy("28-08-2013")) %>% 
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


### Time series features ####
ts.scan.char <- list()
tryCatch(for (i in 1:length(scan)) {
  ts.scan.char[[i]] <- cbind(t(tsfeatures(ts(rowSums(hts.scan[[i]]$bts),frequency = 52))),
                            t(tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,1:6]),frequency = 52))),
                            t(tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,7:12]),frequency = 52))),
                            t(tsfeatures(ts((hts.scan[[i]]$bts[,1:12]),frequency = 52))))
})

chr.scan <- list()
for (i in 1:length(scan)) {
  chr.scan[[i]] <- cbind(t(cbind((tsfeatures(ts(rowSums(hts.scan[[i]]$bts), frequency = 52), 
                                        features = c("fluctanal_prop_r1","frequency", 
                                                     "acf_features","lumpiness","stability","entropy",
                                                     "max_level_shift","max_var_shift","max_kl_shift",
                                                     "crossing_points", 
                                                     "hurst", "unitroot_kpss", "heterogeneity",
                                                     "nonlinearity", "pacf_features"))), 
                             (tsfeatures(ts(rowSums(hts.scan[[i]]$bts), frequency = 52), 
                                        "stl_features", s.window="periodic", robust= TRUE)),
                             (tsfeatures(ts(rowSums(hts.scan[[i]]$bts), frequency = 52),
                                        c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
    select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
           curvature, hurst, lumpiness, stability,seasonal_strength,
           fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
           max_level_shift, max_var_shift, max_kl_shift, crossing_points,
           entropy, arch_acf, garch_acf, arch_r2, garch_r2)), 
    t(cbind((tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,1:6]), frequency = 52), 
                        features = c("fluctanal_prop_r1","frequency",
                                     "acf_features","lumpiness","stability","entropy",
                                     "max_level_shift","max_var_shift","max_kl_shift",
                                     "crossing_points", 
                                     "hurst", "unitroot_kpss", "heterogeneity",
                                     "nonlinearity", "pacf_features"))), 
            (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,1:6]), frequency = 52), 
                        "stl_features", s.window="periodic", robust= TRUE)),
            (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,1:6]), frequency = 52),
                        c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
        select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
               curvature, hurst, lumpiness, stability,seasonal_strength,
               fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
               max_level_shift, max_var_shift, max_kl_shift, crossing_points,
               entropy, arch_acf, garch_acf, arch_r2, garch_r2)),
    
    t(cbind((tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,7:12]), frequency = 52), 
                        features = c("fluctanal_prop_r1","frequency", "entropy",
                                     "acf_features","lumpiness","stability","entropy",
                                     "max_level_shift","max_var_shift","max_kl_shift",
                                     "crossing_points", 
                                     "hurst", "unitroot_kpss", "heterogeneity",
                                     "nonlinearity", "pacf_features"))), 
            (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,7:12]), frequency = 52), 
                        "stl_features", s.window="periodic", robust= TRUE)),
            (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[,7:12]), frequency = 52),
                        c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
        select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
               curvature, hurst, lumpiness, stability,seasonal_strength,
               fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
               max_level_shift, max_var_shift, max_kl_shift, crossing_points,
               entropy_entropy, arch_acf, garch_acf, arch_r2, garch_r2)), 
    
    t(cbind((tsfeatures(ts((hts.scan[[i]]$bts[,1:12]), frequency = 52), 
                        features = c("fluctanal_prop_r1","frequency", "entropy",
                                     "acf_features","lumpiness","stability",
                                     "max_level_shift","max_var_shift","max_kl_shift",
                                     "crossing_points", 
                                     "hurst", "unitroot_kpss", "heterogeneity",
                                     "nonlinearity", "pacf_features"))), 
            (tsfeatures(ts((hts.scan[[i]]$bts[,1:12]), frequency = 52), 
                        "stl_features", s.window="periodic", robust= TRUE)),
            (tsfeatures(ts((hts.scan[[i]]$bts[,1:12]), frequency = 52),
                        c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
        select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
               curvature, hurst, lumpiness, stability,seasonal_strength,
               fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
               max_level_shift, max_var_shift, max_kl_shift, crossing_points,
               entropy, arch_acf, garch_acf, arch_r2, garch_r2)))
  }

saveRDS(chr.scan, "chr.scan.RDS")

Chr_file <- Chr_filen <- NULL
for (i in 1:length(scan)) {
  chr.scan.t <- (t(chr.scan[[i]]))
  #Chr_file <- rbind(Chr_file, cbind(t(chr.scan[[i]]), i))
  rownames(chr.scan.t) = (paste0(colnames(allts(hts.data.scan)), i))
  Chr_filen <- rbind(Chr_filen, chr.scan.t)
}

#write.csv(Chr_file, "Chr_file.csv")
Chr_file <- read.csv("Chr_file.csv")
write.csv(Chr_filen, "Chr_filen.csv")

chr.scan.train <- list()
for (i in 1:length(scan)) {
  chr.scan.train[[i]] <- cbind(t(cbind((tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,]), frequency = 26), 
                                             features = c("fluctanal_prop_r1","frequency", 
                                                          "acf_features","lumpiness","stability","entropy",
                                                          "max_level_shift","max_var_shift","max_kl_shift",
                                                          "crossing_points", 
                                                          "hurst", "unitroot_kpss", "heterogeneity",
                                                          "nonlinearity", "pacf_features"))), 
                                 (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,]), frequency = 26), 
                                             "stl_features", s.window="periodic", robust= TRUE)),
                                 (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,]), frequency = 26),
                                             c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
                             select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
                                    curvature, hurst, lumpiness, stability,seasonal_strength,
                                    fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
                                    max_level_shift, max_var_shift, max_kl_shift, crossing_points,
                                    entropy, arch_acf, garch_acf, arch_r2, garch_r2)), 
                         t(cbind((tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,1:6]), frequency = 26), 
                                             features = c("fluctanal_prop_r1","frequency",
                                                          "acf_features","lumpiness","stability","entropy",
                                                          "max_level_shift","max_var_shift","max_kl_shift",
                                                          "crossing_points", 
                                                          "hurst", "unitroot_kpss", "heterogeneity",
                                                          "nonlinearity", "pacf_features"))), 
                                 (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,1:6]), frequency = 26), 
                                             "stl_features", s.window="periodic", robust= TRUE)),
                                 (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,1:6]), frequency = 26),
                                             c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
                             select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
                                    curvature, hurst, lumpiness, stability,seasonal_strength,
                                    fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
                                    max_level_shift, max_var_shift, max_kl_shift, crossing_points,
                                    entropy, arch_acf, garch_acf, arch_r2, garch_r2)),
                         
                         t(cbind((tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,7:12]), frequency = 26), 
                                             features = c("fluctanal_prop_r1","frequency", "entropy",
                                                          "acf_features","lumpiness","stability","entropy",
                                                          "max_level_shift","max_var_shift","max_kl_shift",
                                                          "crossing_points", 
                                                          "hurst", "unitroot_kpss", "heterogeneity",
                                                          "nonlinearity", "pacf_features"))), 
                                 (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,7:12]), frequency = 26), 
                                             "stl_features", s.window="periodic", robust= TRUE)),
                                 (tsfeatures(ts(rowSums(hts.scan[[i]]$bts[1:55,7:12]), frequency = 26),
                                             c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
                             select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
                                    curvature, hurst, lumpiness, stability,seasonal_strength,
                                    fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
                                    max_level_shift, max_var_shift, max_kl_shift, crossing_points,
                                    entropy_entropy, arch_acf, garch_acf, arch_r2, garch_r2)), 
                         
                         t(cbind((tsfeatures(ts((hts.scan[[i]]$bts[1:55,1:12]), frequency = 26), 
                                             features = c("fluctanal_prop_r1","frequency", "entropy",
                                                          "acf_features","lumpiness","stability",
                                                          "max_level_shift","max_var_shift","max_kl_shift",
                                                          "crossing_points", 
                                                          "hurst", "unitroot_kpss", "heterogeneity",
                                                          "nonlinearity", "pacf_features"))), 
                                 (tsfeatures(ts((hts.scan[[i]]$bts[1:55,1:12]), frequency = 26), 
                                             "stl_features", s.window=26, robust= TRUE)),
                                 (tsfeatures(ts((hts.scan[[i]]$bts[1:55,1:12]), frequency = 26),
                                             c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
                             select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
                                    curvature, hurst, lumpiness, stability,
                                    fluctanal_prop_r1,x_acf1,x_acf10,seas_acf1, e_acf1, e_acf10,
                                    max_level_shift, max_var_shift, max_kl_shift,
                                    entropy, arch_acf, garch_acf, arch_r2, garch_r2,
                                    diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)))
}

saveRDS(chr.scan.train, "chr.scan.train.RDS")

Chr_file_train <- Chr_file_train_n <- NULL
for (i in 1:length(scan)) {
  chr.scan.train.t <- (t(chr.scan.train[[i]]))
  #Chr_file <- rbind(Chr_file, cbind(t(chr.scan[[i]]), i))
  rownames(chr.scan.train.t) = (paste0(colnames(allts(hts.data.scan)), i))
  Chr_file_train_n <- rbind(Chr_file_train_n, chr.scan.train.t)
}


##normalising pos features data

# minmax <- function(x){
#   (x - min(x)) / (max(x) - min(x))
# }
# pos.features.scaled <- sapply(pos.features, minmax)
# 
# 
# bind_rows(pos.features.scaled %>% 
#             as_tibble() %>%
#             select(mean, lumpiness,trend,linearity, hurst, stability, seasonal_strength,nonlinearity, entropy_entropy ) %>% 
#             rename(variance=lumpiness,  entropy=entropy_entropy) %>%
#             mutate(Source="POS"), 
#           exf.features.scaled %>% 
#             as_tibble() %>%
#             select(mean, lumpiness,trend,linearity, hurst, stability, seasonal_strength,nonlinearity ,entropy_entropy) %>% 
#             mutate(Source="Orders") %>% rename(variance=lumpiness, entropy=entropy_entropy)) %>%
#   melt() %>%
  #pivot_longer(names_to = "features", cols = c(POS_trend, Order_trend)) %>%
#  ggplot(aes(value, fill=Source)) + 
 # geom_density(position = "stack") + facet_wrap(~variable, scales = "free_y")




### HF methods ####

#test
# x=scan[[1]]
# y=price[[1]]
# train.obs = 112
# total.obs= 120

train.obs = 55
total.obs= 112

Forecast_file <- NULL
Summary_error <- NULL
error_list <- NULL

hf.arimax <- function(x,y){
  hts.data.scan <- hts(data.frame(x)[1:train.obs,], nodes=list(2, c(6,6)), bnames = bnames)
  hts.data.scan.total <- hts(data.frame(x), nodes=list(2, c(6,6)), bnames =bnames)
  
  hts.data.price <- hts(data.frame(y), nodes=list(2, c(6,6)), bnames = bnames)
  hts.data.scan.test <- hts(data.frame(x)[(1+train.obs) : total.obs, ], nodes=list(2, c(6,6)), bnames = bnames)
  
  fore.wls <- forecast(hts.data.scan, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                       , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0, weights = c("wls"))
  
  fore.mint.shr <- forecast(hts.data.scan, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                            , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0, weights = c("mint"), covariance = "shr")
  
  fore.mint.sam <- forecast(hts.data.scan, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                        , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0, weights = c("mint"), covariance = "sam")
  fore.bu <- forecast(hts.data.scan, h=8, method = "bu", fmethod = "arima", xreg = as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                      , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0)
  fore.mo <- forecast(hts.data.scan, h=8, method = "mo", level = 1, fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                      , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0)
  fore.tdfp <- forecast(hts.data.scan, h=8, method = "tdfp", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                        , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0)
  fore.tdgsa <- forecast(hts.data.scan, h=8, method = "tdgsa", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                      , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0)
  fore.tdgsf <- forecast(hts.data.scan, h=8, method = "tdgsf", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price, levels = 0))[1:train.obs]
                         , newxreg =  as.numeric(aggts(hts.data.price, levels = 0))[(1+train.obs): total.obs], lambda = 0)
  
  ##### Forecasts accuracy #####
  
  models <- c("WLS", "Shr", "Sam", "BU", "MO", "TDFP", "TDGSA", "TDGSF")
  
  ### All time series forecasts #### 
  
  for (model_id in 1:8){
    
    if (model_id==1){
      allf <- data.frame(allts(fore.wls))
    }else if (model_id==2){
      allf <- data.frame(allts(fore.mint.shr))
    }else if (model_id==3){
      allf <- data.frame(allts(fore.mint.sam))
    }else if (model_id==4){
      allf <- data.frame(allts(fore.bu))
    }else if (model_id==5){
      allf <- data.frame(allts(fore.mo))
    }else if (model_id==6){
      allf <- data.frame(allts(fore.tdfp))
    }else if (model_id==7){
      allf <- data.frame(allts(fore.tdgsa))
    }else{
      allf <- data.frame(allts(fore.tdgsf))
    }
    
    colnames(allf) <- paste0("F_", as.character(colnames(allts(hts.data.scan.test))))
    test_sample_c <- as.data.frame(allts(hts.data.scan.test))
    colnames(test_sample_c) <- paste0("A_", as.character(colnames(test_sample_c)))
    #Combine with actual values
    allf <- cbind(allf, test_sample_c) 
    allf$model <- model_id
    allf$tsid <- tsid
    Forecast_file <- rbind(Forecast_file, allf)



    structure <- window(hts.data.scan)
    m_series <- nrow(smatrix(structure)) #Total number of series at all levels
    b_series <- length(structure$labels[[(length(structure$labels)-1)+1]]) #Number of series at the bottom level
    error_list <- NULL
  for (j in 1:ncol(allts(hts.data.scan.test))){
    error_list <- rbind(error_list,error_cal(allts(hts.data.scan)[,j],allf[,j],allts(hts.data.scan.test)[,j],ppy=1))
  }
  Errors <- data.frame(error_list) ; colnames(Errors) <- c("MASE","AMSE","RMSSE")
  Errors$Tags <- as.character(unlist(colnames(allts(hts.data.scan.test))))
  for (z in 1:nrow(Errors)){
    Errors$Tags[z] <- substr(Errors$Tags[z],1, nchar(Errors$Tags[z]))
  }
  Errors$Level <- NA
  for (idd in 1:nrow(Errors)){
    for (ldd in 1:length(structure$labels)){
      if (Errors$Tags[idd] %in% structure$labels[[ldd]]){
        Errors$Level[idd] <- ldd
      }
    }
  }
  Errors$model <- model_id
  Errors$tsid <- tsid
  Summary_error <- rbind(Summary_error, Errors)
  }
  mylist <- list(Summary_error, Forecast_file)
  return(mylist)
}

all_results <- list()

for (tsid in 1:length(scan)) {
  all_results [[tsid]] <- hf.arimax(scan[[tsid]], price[[tsid]])
}

saveRDS(all_results, "all_results.RDS")

F_file <- NULL

# Chr_file <- Chr_filen <- NULL
# for (i in 1:length(scan)) {
#   F_file.t <- (t(chr.scan[[i]][[2]]))
#   #Chr_file <- rbind(Chr_file, cbind(t(chr.scan[[i]]), i))
#   rownames(chr.scan.t) = (paste0(colnames(allts(hts.data.scan)), i))
#    F_file.t <- rbind(Chr_filen, chr.scan.t)
# }


for (tsid in 1:length(scan)) {
  all_results[[tsid]][[2]][,c(1:15,31,32)]
  F_file <- rbind(F_file, all_results[[tsid]][[2]])
}
write.csv(F_file[,c(1:15,31,32)], "F_file.csv")
Actual_file <- F_file[,16:32]
write.csv(F_file[,16:32], "Actual_file.csv")

Acc_file <- NULL
for (tsid in 1:length(scan)) {
  Acc_file <- rbind(Acc_file, all_results[[tsid]][[1]])
}
write.csv(Acc_file, "Acc_file.csv")

Acc_file_n <- Acc_file %>% as_tibble() %>%
  mutate(name= paste0(Tags,tsid))




###3 training model and validation set
train.obs = 55
total.obs= 112

Forecast_file_train <- NULL
Summary_error_train <- NULL
error_list_train <- NULL


data <- hts.scan
data.price <- hts.price
#structure <- window(data, start = c(1998,1))

structure.n <-  window(data, start = 1)
m_series <- nrow(smatrix(structure.n)) #Total number of series at all levels
b_series <- length(structure.n$labels[[(length(structure.n$labels)-1)+1]]) #Number of series at the bottom level
input <- NULL

#Get base forecasts for training and forecast purposes
Forecast_file = Residuals_file <- NULL
Forecast_file_train <- NULL
Summary_error_train <- NULL
chr_train_stack <- ML_data_train <- NULL
Selected_model <- NULL
tsid.outlier <- c(3,5:10)
for (tsid in 6:55) {
  
  counter <- 0
  origin <- 26 #Starting origin (26weeks)
  tslength <- 85 #Available observations
  #tslength <- length(data[[tsid]]$bts[,1]) #Available observations
  fh = 8 #Forecasting horizon considered
  
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
   hts.data.scan.test <- hts(((test_sample)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
   
   hts.data.price.train <- hts(((train_sample_p)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
   hts.data.price.test <- hts(((test_sample_p)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
   
  # hts.data.scan <- hts(data.frame(x)[1:train.obs,], nodes=list(2, c(6,6)), bnames = bnames)
  # hts.data.scan.total <- hts(data.frame(x), nodes=list(2, c(6,6)), bnames =bnames)
  # 
  # hts.data.price <- hts(data.frame(y), nodes=list(2, c(6,6)), bnames = bnames)
  # hts.data.scan.test <- hts(data.frame(x)[(1+train.obs) : total.obs, ], nodes=list(2, c(6,6)), bnames = bnames)
  # 
 # fore.wls.ro <- forecast(hts.data.scan.train, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
  #                     , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)), lambda = 0, weights = c("wls"))
  
  fore.wls.ro <- tryCatch(forecast(hts.data.scan.train, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                       , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)), weights = c("wls")))
  
  fore.mint.shr.ro <- tryCatch(forecast(hts.data.scan.train, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                            , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)), weights = c("mint"), covariance = "shr"))
  
  # fore.mint.sam.ro <- tryCatch(forecast(hts.data.scan.train, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0)+0.01)
   #                          , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)+.01), weights = c("mint"), covariance = "sam"))
  fore.mint.sam.ro <- fore.mint.shr.ro
  fore.bu.ro <- forecast(hts.data.scan.train, h=8, method = "bu", fmethod = "arima", xreg = as.numeric(aggts(hts.data.price.train, levels = 0))
                      , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
  
  fore.mo.ro <- forecast(hts.data.scan.train, h=8, method = "mo", level = 1, fmethod = "arima",xreg =as.numeric(aggts(hts.data.price.train, levels = 0))
                      , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
  
  fore.tdfp.ro <- forecast(hts.data.scan.train, h=8, method = "tdfp", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                        , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
  
  fore.tdgsa.ro <- forecast(hts.data.scan.train, h=8, method = "tdgsa", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                         , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
  
  fore.tdgsf.ro <- forecast(hts.data.scan.train, h=8, method = "tdgsf", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                         , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
  
  ##Time series features ####
  
  ##### Forecasts accuracy #####
  
  #models <- c("WLS", "Shr", "Sam", "BU", "MO", "TDFP", "TDGSA", "TDGSF")
  test_sample_c <- as.data.frame(allts(hts.data.scan.test))
  colnames(test_sample_c) <- paste0("A_", as.character(colnames(test_sample_c)))
  ### All time series forecasts #### 
  for (hs in 1:b_series) {
    Min_Error <- 999999.99
    chr.train <- cbind(t(cbind((tsfeatures(ts((hts.data.scan.train$bts[,hs]), frequency = 1), 
                                           features = c("fluctanal_prop_r1","frequency", "entropy",
                                                        "acf_features","lumpiness","stability",
                                                        "max_level_shift","max_var_shift","max_kl_shift",
                                                        
                                                        "hurst", "unitroot_kpss", "heterogeneity",
                                                        "nonlinearity", "pacf_features"))), 
                               (tsfeatures(ts((hts.data.scan.train$bts[,hs]), frequency = 1), 
                                           "stl_features", s.window="periodic", robust= TRUE)),
                               (tsfeatures(ts((hts.data.scan.train$bts[,hs]), frequency = 1),
                                           c("mean","var"), scale=FALSE, na.rm=TRUE)))%>%
                           select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
                                  curvature, hurst, lumpiness, stability,
                                  fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
                                  max_level_shift, max_var_shift, max_kl_shift,
                                  entropy, arch_acf, garch_acf, arch_r2, garch_r2,
                                  diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)))
    
    repML=counter+1
    originML=length(train_sample[,1])
    
    chr_train_stack <- rbind(chr_train_stack, cbind(chr.train,repML,originML,tsid, hs))
    
  for (model_id in 1:8){
    
    if (model_id==1){
      allf <- data.frame(allts(fore.wls.ro))
    }else if (model_id==2){
      allf <- data.frame(allts(fore.mint.shr.ro))
    }else if (model_id==3){
      allf <- data.frame(allts(fore.mint.sam.ro))
    }else if (model_id==4){
      allf <- data.frame(allts(fore.bu.ro))
    }else if (model_id==5){
      allf <- data.frame(allts(fore.mo.ro))
    }else if (model_id==6){
      allf <- data.frame(allts(fore.tdfp.ro))
    }else if (model_id==7){
      allf <- data.frame(allts(fore.tdgsa.ro))
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
    
    
    
    #structure <- window(hts.data.scan)
    #m_series <- nrow(smatrix(structure)) #Total number of series at all levels
    #b_series <- length(structure$labels[[(length(structure$labels)-1)+1]]) #Number of series at the bottom level
    error_list_train <- NULL
    #for (j in 1:ncol(allts(hts.data.scan.test))){
    error_list_train <- error_cal(allts(hts.data.scan.train)[,hs],all[,hs],allts(hts.data.scan.test)[,hs],ppy=1)
    #}
    #for (j in 1:ncol(allts(hts.data.scan.test))){
      #error_list_train <- rbind(error_list_train,error_cal(allts(hts.data.scan.train)[,hs],all[,hs],allts(hts.data.scan.test)[,hs],ppy=1))
    #}
    Errors <- data.frame(error_list_train) ; colnames(Errors) <- ("MASE")#,"AMSE","RMSSE")
    # Errors$Tags <- as.character(unlist(colnames(allts(hts.data.scan.test))))
    # for (z in 1:nrow(Errors)){
    #   Errors$Tags[z] <- substr(Errors$Tags[z],1, nchar(Errors$Tags[z]))
    # }
    # Errors$Level <- NA
    # for (idd in 1:nrow(Errors)){
    #   for (ldd in 1:length(structure$labels)){
    #     if (Errors$Tags[idd] %in% structure$labels[[ldd]]){
    #       Errors$Level[idd] <- ldd
    #     }
    #   }
    # }
    Errors$model <- model_id
    Errors$tsid <- tsid
    Errors$origin <- origin
    Errors$hs <- hs
    Summary_error_train <- rbind(Summary_error_train, Errors)
    
    
    #Min_Error_MASE <- (Errors[hs,1])
    
    Min_Error_MASE <- (Errors[1,])
    
    #Min_Error <- Min_Error_MASE
    if (Min_Error_MASE < Min_Error){
        Min_Error <- Min_Error_MASE 
    Selected_model <- model_id} else {
      Min_Error <- Min_Error
      Selected_model <- Selected_model
    }
  
  
  
   
     #t(cbind((tsfeatures(ts(rowSums(hts.data.scan.train$bts), frequency = 1), 
  #                                                  features = c("fluctanal_prop_r1","frequency", 
  #                                                               "acf_features","lumpiness","stability","entropy",
  #                                                               "max_level_shift","max_var_shift","max_kl_shift",
  #                                                               
  #                                                               "hurst", "unitroot_kpss", "heterogeneity",
  #                                                               "nonlinearity", "pacf_features"))), 
  #                                      (tsfeatures(ts(rowSums(hts.data.scan.train$bts), frequency = 1), 
  #                                                  "stl_features", s.window="periodic", robust= TRUE)),
  #                                      (tsfeatures(ts(rowSums(hts.data.scan.train$bts), frequency = 1),
  #                                                  c("mean","var"), scale=FALSE, na.rm=TRUE))) %>%
  #                        select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
  #                               curvature, hurst, lumpiness, stability,
  #                               fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
  #                               max_level_shift, max_var_shift, max_kl_shift,
  #                               entropy, arch_acf, garch_acf, arch_r2, garch_r2,
  #                               diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)), 
  #                              t(cbind((tsfeatures(ts(rowSums(hts.data.scan.train$bts[,1:6]), frequency = 1), 
  #                                                  features = c("fluctanal_prop_r1","frequency",
  #                                                               "acf_features","lumpiness","stability","entropy",
  #                                                               "max_level_shift","max_var_shift","max_kl_shift",
  #                                                               
  #                                                               "hurst", "unitroot_kpss", "heterogeneity",
  #                                                               "nonlinearity", "pacf_features"))), 
  #                                      (tsfeatures(ts(rowSums(hts.data.scan.train$bts[,1:6]), frequency = 1), 
  #                                                  "stl_features", s.window="periodic", robust= TRUE)),
  #                                      (tsfeatures(ts(rowSums(hts.data.scan.train$bts[,1:6]), frequency = 1),
  #                                                  c("mean","var"), scale=FALSE, na.rm=TRUE)))%>%
  #                                  select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
  #                                         curvature, hurst, lumpiness, stability,
  #                                         fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
  #                                         max_level_shift, max_var_shift, max_kl_shift,
  #                                         entropy, arch_acf, garch_acf, arch_r2, garch_r2,
  #                                         diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)),
  #                              
  #                              t(cbind((tsfeatures(ts(rowSums(hts.data.scan.train$bts[,7:12]), frequency = 1), 
  #                                                  features = c("fluctanal_prop_r1","frequency", "entropy",
  #                                                               "acf_features","lumpiness","stability","entropy",
  #                                                               "max_level_shift","max_var_shift","max_kl_shift",
  #                                                               
  #                                                               "hurst", "unitroot_kpss", "heterogeneity",
  #                                                               "nonlinearity", "pacf_features"))), 
  #                                      (tsfeatures(ts(rowSums(hts.data.scan.train$bts[,7:12]), frequency = 1), 
  #                                                  "stl_features", s.window="periodic", robust= TRUE)),
  #                                      (tsfeatures(ts(rowSums(hts.data.scan.train$bts[,7:12]), frequency = 1),
  #                                                  c("mean","var"), scale=FALSE, na.rm=TRUE)))%>%
  #                                  select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
  #                                         curvature, hurst, lumpiness, stability,
  #                                         fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
  #                                         max_level_shift, max_var_shift, max_kl_shift,
  #                                         entropy_entropy, arch_acf, garch_acf, arch_r2, garch_r2,
  #                                         diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)), 

  }
    ML_data <- cbind(t(rowMeans(chr.train)), Selected_model, repML, originML, tsid, hs)
    ML_data_train <- rbind(ML_data_train, ML_data)
  }
  counter <- counter + 1
  } 
}


# for (tsid in 1:length(scan)) {
#   all_results[[tsid]][[2]][,c(1:15,31,32)]
#   F_file <- rbind(F_file, all_results[[tsid]][[2]])
# }
# write.csv(F_file[,c(1:15,31,32)], "F_file.csv")
# Actual_file <- F_file[,16:32]
# write.csv(F_file[,16:32], "Actual_file.csv")
# 
# Acc_file <- NULL
# for (tsid in 1:length(scan)) {
#   Acc_file <- rbind(Acc_file, all_results[[tsid]][[1]])
# }
# write.csv(Acc_file, "Acc_file.csv")
# 
# Acc_file_n <- Acc_file %>% as_tibble() %>%
#   mutate(name= paste0(Tags,tsid))
# 

# features <- c("mean", "var", "unitroot_kpss", "trend", "linearity", "nonlinearity",
#   "curvature", "hurst", "lumpiness", "stability","seasonal_strength",
#   "fluctanal_prop_r1","x_acf1","x_acf10","seas_acf1", "e_acf1", "e_acf10",
#   "max_level_shift", "max_var_shift", "max_kl_shift", "crossing_points",
#   "entropy_entropy", "arch_acf", "garch_acf", "arch_r2", "garch_r2")


##########################Conditional reconciliation#############################################

ML_data_train.dataframe <- as.data.frame(ML_data_train)
Forecast_file_test <- NULL
chr_test_stack <- NULL
Errors_test1 <- NULL
error_list_test1 <- NULL
Summary_error_test1 <- NULL
Summary_error_test_f <- NULL
recon_bottom_f_all <- NULL
for (tsid in 1:8) {
      
  counter <- 0
  origin_id <- 77 #Starting origin (26weeks)
  tslength_id <- 120 #Available observations
  #tslength <- length(data[[tsid]]$bts[,1]) #Available observations
  fh = 8 #Forecasting horizon considered
  
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
    hts.data.scan.test <- hts(((test_sample)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
    
    hts.data.price.train <- hts(((train_sample_p)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
    hts.data.price.test <- hts(((test_sample_p)[,4:15]), nodes=list(2, c(6,6)), bnames = bnames)
    
    
    ## forecasting the series ###
    fore.wls.ro <- (forecast(hts.data.scan.train, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                             , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)), weights = c("wls")))
    
    fore.mint.shr.ro <- (forecast(hts.data.scan.train, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                                  , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)), weights = c("mint"), covariance = "shr"))
    
    #fore.mint.sam.ro <- (forecast(hts.data.scan.train, h=8, method = "comb", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
     #                         , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)+.01), weights = c("mint"), covariance = "sam"))
    fore.mint.sam.ro <- fore.mint.shr.ro
    
    fore.bu.ro <- forecast(hts.data.scan.train, h=8, method = "bu", fmethod = "arima", xreg = as.numeric(aggts(hts.data.price.train, levels = 0))
                           , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
    fore.mo.ro <- forecast(hts.data.scan.train, h=8, method = "mo", level = 1, fmethod = "arima",xreg =as.numeric(aggts(hts.data.price.train, levels = 0))
                           , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
    fore.tdfp.ro <- forecast(hts.data.scan.train, h=8, method = "tdfp", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                             , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
    fore.tdgsa.ro <- forecast(hts.data.scan.train, h=8, method = "tdgsa", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                              , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
    fore.tdgsf.ro <- forecast(hts.data.scan.train, h=8, method = "tdgsf", fmethod = "arima",xreg =  as.numeric(aggts(hts.data.price.train, levels = 0))
                              , newxreg =  as.numeric(aggts(hts.data.price.test, levels = 0)))
    
    Errors_test1 <- NULL
    error_list_test1 <- NULL
    
    for (model_id in 1:8){
      
      if (model_id==1){
        all_test_f <- data.frame(allts(fore.wls.ro))
      }else if (model_id==2){
        all_test_f <- data.frame(allts(fore.mint.shr.ro))
      }else if (model_id==3){
        all_test_f <- data.frame(allts(fore.mint.sam.ro))
      }else if (model_id==4){
        all_test_f <- data.frame(allts(fore.bu.ro))
      }else if (model_id==5){
        all_test_f <- data.frame(allts(fore.mo.ro))
      }else if (model_id==6){
        all_test_f <- data.frame(allts(fore.tdfp.ro))
      }else if (model_id==7){
        all_test_f <- data.frame(allts(fore.tdgsa.ro))
      }else{
        all_test_f <- data.frame(allts(fore.tdgsf.ro))
      }
      
      # colnames(all_test_f) <- paste0("F_", as.character(colnames(allts(hts.data.scan.test))))
      # 
      # #Combine with actual values
      # all_test_f <- cbind(all_test_f, test_sample_c) 
      # all_test_f$model <- model_id
      # all_test_f$tsid <- tsid
      # all_test_f$rep <- counter +1
      # all_test_f$fh <- c(1:fh)
      # all_test_f$period <- c((length(train_sample[,1])+1) : (length(train_sample[,1])+fh))
      # all_test_f$origin <- length(train_sample[,1])
      # 
      # Forecast_file_test <- rbind(Forecast_file_test, all_test_f)
      error_list_test1 <- NULL
      for (j in 1:ncol(allts(hts.data.scan.test))){
        error_list_test1 <- rbind(error_list_test1,error_cal(allts(hts.data.scan.train)[,j],(all_test_f)[,j],allts(hts.data.scan.test)[,j],ppy=1))
        #error_list_train <- rbind(error_list_train,error_cal(allts(hts.data.scan.train)[,hs],all[,hs],allts(hts.data.scan.test)[,hs],ppy=1))
      }
      Errors_test1 <- cbind(Errors_test1,error_list_test1)
    }
      #rownames(Errors_test1) <-  as.character(unlist(colnames(allts(hts.data.scan.test))))#,"AMSE","RMSSE")
    Errors_test1 <- as.data.frame(Errors_test1)  
    colnames(Errors_test1) <- models
      # Errors_test1$Tags <- as.character(unlist(colnames(allts(hts.data.scan.test))))
      # for (z in 1:nrow(Errors_test1)){
      #   Errors_test1$Tags[z] <- substr(Errors_test1$Tags[z],1, nchar(Errors_test1$Tags[z]))
      # }
      # Errors_test1$Level <- NA
      # for (idd in 1:nrow(Errors_test1)){
      #   for (ldd in 1:length(structure$labels)){
      #     if (Errors_test1$Tags[idd] %in% structure$labels[[ldd]]){
      #       Errors_test1$Level[idd] <- ldd
      #     }
      #   }
      # }
      #Errors_test1$model <- model_id
      Errors_test1$tsid <- tsid
      #Errors_test1$origin_id <- origin_id
      Errors_test1$counter <- counter+1
      Errors_test1$level <- c(1,2,2,3,3,3,3,3,3,3,3,3,3,3,3)
      Summary_error_test1 <- rbind(Summary_error_test1, Errors_test1) 
      
    ## forecasting the bottom series with ML ###
      recon_f <- NULL
      recon_bottom_f <- NULL
    for (hs in 1:b_series) {
      #Min_Error <- 999999.99
      chr.test <- cbind(t(cbind((tsfeatures(ts((hts.data.scan.train$bts[,hs]), frequency = 1), 
                                             features = c("fluctanal_prop_r1","frequency", "entropy",
                                                          "acf_features","lumpiness","stability",
                                                          "max_level_shift","max_var_shift","max_kl_shift",
                                                          
                                                          "hurst", "unitroot_kpss", "heterogeneity",
                                                          "nonlinearity", "pacf_features"))), 
                                 (tsfeatures(ts((hts.data.scan.train$bts[,hs]), frequency = 1), 
                                             "stl_features", s.window="periodic", robust= TRUE)),
                                 (tsfeatures(ts((hts.data.scan.train$bts[,hs]), frequency = 1),
                                             c("mean","var"), scale=FALSE, na.rm=TRUE)))%>%
                             select(mean, var, unitroot_kpss, trend, linearity, nonlinearity,
                                    curvature, hurst, lumpiness, stability,
                                    fluctanal_prop_r1,x_acf1,x_acf10, e_acf1, e_acf10,
                                    max_level_shift, max_var_shift, max_kl_shift,
                                    entropy, arch_acf, garch_acf, arch_r2, garch_r2,
                                    diff1_acf1,diff1_acf10,diff2_acf1,diff2_acf10)))
      
      # repML=counter+1
      # originML=length(train_sample[,1])
      
      chr_test_stack <- rbind(chr_test_stack, cbind(chr.test))#,repML,originML,tsid, hs))
      
     
        #Combine with actual values
        # all <- cbind(allf, test_sample_c) 
        # all$model <- model_id
        # all$tsid <- tsid
        # all$rep <- counter +1
        # all$fh <- c(1:fh)
        # all$period <- c((length(train_sample[,1])+1) : (length(train_sample[,1])+fh))
        # all$origin <- length(train_sample[,1])
        
        #Forecast_file_train <- rbind(Forecast_file_train, all)
        #
        
        ### ML model to classify the appropriate reconciliation method ####
      Reconcile_data <-  ML_data_train.dataframe[ML_data_train.dataframe$tsid==tsid,]
      Reconcile_data <- Reconcile_data[Reconcile_data$hs==hs,]
      
      
      y_value <- t(as.factor(Reconcile_data[,28]))
      ML_model <- randomForest(x = Reconcile_data[,1:27] , y= y_value)
      ML_output <- predict(ML_model, newdata = t(chr.train))
      
      
         if (ML_output==1){
              recon_f <- data.frame(allts(fore.wls.ro))[,hs+3]
          }else if (ML_output==2){
              recon_f <- data.frame(allts(fore.mint.shr.ro))[,hs+3]
          }else if (ML_output==3){
              recon_f <- data.frame(allts(fore.mint.shr.ro))[,hs+3]
                     #recon_f <- data.frame(allts(fore.mint.sam.ro))[,hs+3]
          }else if (ML_output==4){
                      recon_f <- data.frame(allts(fore.bu.ro))[,hs+3]
          }else if (ML_output==5){
            recon_f <- data.frame(allts(fore.mo.ro))[,hs+3]
          }else if (ML_output==6){
                      recon_f <- as.data.frame(allts(fore.tdfp.ro))[,hs+3]
          }else if (ML_output==7){
                recon_f <- data.frame(allts(fore.tdgsa.ro))[,hs+3]
          }else{
              recon_f <- data.frame(allts(fore.tdgsf.ro))[,hs+3]
          }
      recon_bottom_f <- cbind(recon_bottom_f,recon_f)
      }
    colnames(recon_bottom_f) <- paste0("F_", as.character(colnames(hts.data.scan.test$bts)))
    
    hts.bu.forecast <- hts((recon_bottom_f), nodes=list(2, c(6,6)), bnames = bnames)
    
    recon_bottom_f_all <- rbind(recon_bottom_f_all, cbind(recon_bottom_f, tsid, counter)) 
    ##Time series features ####
    
    ##### Forecasts accuracy #####
    
    #models <- c("WLS", "Shr", "Sam", "BU", "MO", "TDFP", "TDGSA", "TDGSF")
    test_sample <- as.data.frame(allts(hts.data.scan.test))
    colnames(test_sample) <- paste0("A_", as.character(colnames(test_sample)))
    ### All time series forecasts #### 
        error_list_test <- NULL
        #for (j in 1:ncol(allts(hts.data.scan.test))){
        #error_list_test <- error_cal(allts(hts.data.scan.train)[,hs],recon_bottom_f[,hs],allts(hts.data.scan.test)[,hs],ppy=1)
        #}
        for (j in 1:ncol(allts(hts.data.scan.test))){
          error_list_test <- rbind(error_list_test,error_cal(allts(hts.data.scan.train)[,j],allts(hts.bu.forecast)[,j],allts(hts.data.scan.test)[,j],ppy=1))
          #error_list_train <- rbind(error_list_train,error_cal(allts(hts.data.scan.train)[,hs],all[,hs],allts(hts.data.scan.test)[,hs],ppy=1))
        }
        Errors_test <- data.frame(error_list_test) ; colnames(Errors_test) <- ("CR")#,"AMSE","RMSSE")
        # Errors_test$Tags <- as.character(unlist(colnames(allts(hts.data.scan.test))))
        # for (z in 1:nrow(Errors_test)){
        #   Errors_test$Tags[z] <- substr(Errors_test$Tags[z],1, nchar(Errors_test$Tags[z]))
        # }
        # Errors_test$Level <- NA
        # for (idd in 1:nrow(Errors_test)){
        #   for (ldd in 1:length(structure$labels)){
        #     if (Errors_test$Tags[idd] %in% structure$labels[[ldd]]){
        #       Errors_test$Level[idd] <- ldd
        #     }
        #   }
        # }
        # Errors_test$model <- model_id
        # Errors_test$tsid <- tsid
        # Errors_test$origin_id <- origin_id
        # Errors_test$hs <- hs
        Summary_error_test_f <- rbind(Summary_error_test_f, cbind(Errors_test,Errors_test1))
        counter <- counter + 1
    }
   
  }
levels=c(1,2,2,3,3,3,3,3,3,3,3,3,3,3,3)
Allerrors <- cbind(Summary_error_test_f,rep(levels,36))
colMeans(Summary_error_test_f[Summary_error_test_f$level==1,])
#######################################################################
#Reconcile 
########################################################################
#setwd("C:/Users/FSU_Team/Desktop/Tourism/Reestimate")
#load("WD_ets_arima_theta.Rdata")
DoTheJob_ML <- function(x){
  dftest <- cbind(act_ml_array[,x], frc_ml_array) ; colnames(dftest)[1] <- "Y"
  model_ML <- randomForest(formula = Y ~ .,  data = dftest, ntree = 100)
  return(model_ML)
}

Summary_error<- c()
Forecasts_ml_saved <- NULL

originid= 52
for (sid in 1:55) {
for (originid in seq(52,max(Forecast_file$origin),1)){
  
  ## Mahdi: I replaced series type ID with model type ID, cause we have only one model
  model_type <- "ARIMA" ; series_type_id <- sid #Consider a base model of type ETS (1), ARIMA (2) or Theta (3)
  
  tempf <- Forecast_file[Forecast_file$origin == originid,]
  tempf <- tempf[tempf$series==series_type_id,]
  Outsample <- tempf[,(1:m_series)+m_series]
  Forecasts <- tempf[,1:m_series]
  Insample <- aggts(data[[sid]])[(1:originid),]
  
  Forecasts_st <- t(reconcile(Forecasts)) #WLS-Structural
  Forecasts_bu <- t(bu(Forecasts)) #BU
  Forecasts_td <- t(td(Forecasts)) #TD
  
  res = Residuals_file[(Residuals_file$series==series_type_id)&(Residuals_file$origin==originid),]
  res$series = res$origin = res$period = res$rep <- NULL
  row.names(res) <- NULL 
  res <- as.matrix(res) ; forc <- as.matrix(Forecasts)
  Forecasts_mint <- MinT(forc, get_nodes(structure.n), residual = res, 
              covariance = "shr", keep = "all", algorithms = "lu")
  
  #Train ML
  data_ml <- Forecast_file[Forecast_file$period <= originid,]
  data_ml <- data_ml[data_ml$fh==1,] #Consider just one step ahead forecasts
  data_ml <- data_ml[data_ml$series==series_type_id,] 
  data_ml$origin = data_ml$period = data_ml$rep = data_ml$fh = data_ml$model = data_ml$series <- NULL
  frc_ml_array <- data_ml[,c(1:m_series)] #Forecasts for all levels
  act_ml_array <- data_ml[,c(1:m_series)+m_series] #Actual values of bottom level
  act_ml_array <- act_ml_array[,c((m_series-b_series+1):m_series)]
  
  Models_ml <- NULL
  for (x in 1:ncol(act_ml_array)){
    Models_ml[length(Models_ml)+1] <- list(DoTheJob_ML(x))
  }
  
  Forecasts_ml_b <- NULL 
  for (iid in 1:b_series){
    Forecasts_ml_b <- cbind(Forecasts_ml_b, as.matrix(predict(Models_ml[[iid]] , newdata = Forecasts)))
  }
  Forecasts_ml <- Forecasts
  Forecasts_ml[,1:(ncol(Forecasts_ml)-b_series)] <- 0
  Forecasts_ml[,(ncol(Forecasts_ml)-b_series+1):ncol(Forecasts_ml)] <- Forecasts_ml_b
  Forecasts_ml_b <- NULL
  Forecasts_ml <- t(bu(Forecasts_ml))
  
  # just saving ML forecasts
  Forecasts_ml_s <- cbind(Forecasts_ml, series_type_id,originid )
  Forecasts_ml_saved <- (rbind(Forecasts_ml_s, Forecasts_ml_saved))

    ##########################################################################
  #End ML reconciliation
  ##########################################################################
  for (mid in 1:5){
    if (mid==1){
      Forecasts_temp <- Forecasts_bu
    }else if (mid==2){
      Forecasts_temp <- Forecasts_td
    }else if (mid==3){
      Forecasts_temp <- Forecasts_st
    }else if (mid==4){
      Forecasts_temp <- Forecasts_mint
    }else{
      Forecasts_temp <- Forecasts_ml
    }
    
    # forecast_rf <- NULL
    # forecast_rf1 <- rbind(forecast_rf1, Forecasts_ml)
    error_list <- NULL
    for (i in 1:ncol(Outsample)){
      error_list <- rbind(error_list, error_cal(Insample[,i],Forecasts_temp[,i],Outsample[,i],ppy=12))
    }
    Errors <- data.frame(error_list) ; colnames(Errors) <- c("MASE","AMSE","RMSSE")
    Errors$Tags <- as.character(unlist(colnames(Outsample)))
    for (i in 1:nrow(Errors)){
      Errors$Tags[i] <- substr(Errors$Tags[i],3, nchar(Errors$Tags[i]))
    }
    Errors$Level <- NA
    for (idd in 1:nrow(Errors)){
      for (ldd in 1:length(structure.n$labels)){
        if (Errors$Tags[idd] %in% structure.n$labels[[ldd]]){
          Errors$Level[idd] <- ldd
        }
      }
    }
    Errors$Tags <- NULL
    if (mid==1){
      Errors$mid <- "BU"
    }else if (mid==2){
      Errors$mid <- "TD"
    }else if (mid==3){
      Errors$mid <- "SS"
    }else if (mid==4){
      Errors$mid <- "MinT"
    }else if (mid==5){
      Errors$mid <- "ML"
    }
    
    Errors$origin <- originid
    Errors$series <- series_type_id
    Summary_error <- rbind(Summary_error, Errors)
  }
}
}

write.csv(Forecasts_ml_saved, "Forecast_ML_RF.csv")

Errors_plot <- ddply(Summary_error, .(mid, Level,origin), colwise(mean))
Errors_plot$mid <- as.factor(Errors_plot$mid)
Errors_plot$Level <- as.factor(Errors_plot$Level)
ggplot(data=Errors_plot, aes(x=origin, y=MASE, shape=mid, colour=Level)) + geom_line() + geom_point()


Errors_agg <- ddply(Summary_error, .(mid, Level), colwise(mean))
Errors_agg$origin <- NULL
ddply(Errors_agg, .(mid), colwise(mean))

write.csv(Summary_error, paste0("Summary_error_ML_RF",model_type,".csv"), row.names = F)


SummaryErrorAll <- read_xlsx("SummaryErrorAll.xlsx")
## Plots for the paper

plot.gts(hts.scan[[30]], labels = FALSE)


SummaryErrorAll %>% 
  select(-series, -origin, - AMSE, -RMSSE) %>% 
  group_by(Level, Method) %>%
  select(MASE, Level, Method) %>%
  ggplot(aes(x=Method, y=MASE, fill = Method))+
  labs(title="", x="")+
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_boxplot() + 
  facet_grid(~Level) + 
  coord_flip() +
  facet_grid(reorder(Level, -desc(Level))~.) +
  theme_bw(base_size = 14)


SummaryErrorAll %>% 
  select(-series, -origin, - MASE, -RMSSE) %>% 
  group_by(Level, Method) %>%
  select(AMSE, Level, Method) %>%
  ggplot(aes(x=Method, y=AMSE, fill = Method))+
  labs(title="", x="")+
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_boxplot() + 
  facet_grid(~Level) + 
  coord_flip() +
  facet_grid(reorder(Level, -desc(Level))~.) +
  theme_bw(base_size = 14)


SummaryErrorAll %>% 
  select(-series, -origin, - AMSE, -MASE) %>% 
  group_by(Level, Method) %>%
  select(RMSSE, Level, Method) %>%
  ggplot(aes(x=Method, y=RMSSE, fill = Method))+
  #labs(title="", x="")+
  theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) + 
  scale_y_continuous(limits = c(0,3)) + 
  geom_boxplot() + 
  facet_grid(~Level) + 
  coord_flip() +
  facet_grid(reorder(Level, -desc(Level))~.) +
  theme_bw(base_size = 14)

m1 <- SummaryErrorAll %>% 
     select(-series, -origin, - AMSE, -RMSSE) %>% 
     group_by(Level, Method) %>%
     select(MASE, Level, Method) %>% rename(value = MASE) %>% mutate(acc = "MASE") 
m2 <- SummaryErrorAll %>% 
     select(-series, -origin, - MASE, -RMSSE) %>% 
     group_by(Level, Method) %>%
     select(AMSE, Level, Method) %>% rename(value = AMSE) %>% mutate(acc = "AMSE") 
m3 <- SummaryErrorAll %>% 
     select(-series, -origin, - MASE, -AMSE) %>% 
     group_by(Level, Method) %>%
     select(RMSSE, Level, Method) %>% rename(value = RMSSE) %>% mutate(acc = "RMSSE")


m <- bind_rows(m1,m2,m3)

m %>% 
  ggplot(aes(x = factor(Method, levels = c("BU", "TD", "COM-SS", "COM-SHR", "ML-RF", "ML-XGB")), y = value, fill = Method)) +
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0,1.0)) +
  geom_boxplot() + facet_grid(Level~ acc) +
  labs(title = "",subtitle = "",  y = "", x = "") + 
  theme(axis.text.x = element_text(angle=50, hjust=1))


SummaryErrorAll %>% 
  select(-series, -origin) %>% 
  group_by(Level, Method) %>%
  gather() %>%
  ungroup() %>%
  select(RMSSE, Level, Method) %>%
  ggplot(aes(x=Method, y=RMSSE, fill = Method))+
  labs(title="", x="")+
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0,3)) + 
  geom_boxplot() + 
  facet_grid(~Level) + 
  coord_flip() +
  facet_grid(reorder(Level, -desc(Level))~.) +
  facet_wrap(. ~ Method)
  theme_bw(base_size = 14)




