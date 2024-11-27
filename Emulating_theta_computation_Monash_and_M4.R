#----------------------------------------------------------------------------------
# This file emulates the computation of mean MASE for the Theta model in the Monash 
# Time Series Forecasting Repository (https://forecastingdata.org/) and in the M4
# competition for the M4 yearly dataset. Computations can be easily changed to use 
# the quarterly and monthly datasets of the M4 competition.
# It also uses the same steps to compute the mean MASE of three of the model used
# in the paper.
# We recommend, at least the first time, execute code chunks of the file starting
# from the beginning
#----------------------------------------------------------------------------------
library(forecast)
library(tsfknn)
library(tsibble)
library(lubridate)

# The implementation of the theta method used in the Monash Repository consists
# of the following two functions, obtained from 
# https://github.com/rakshitha123/TSForecasting/blob/master/models/local_univariate_models.R
# Calculate theta forecasts
get_theta_forecasts <-function(time_series, forecast_horizon){
  tryCatch(
    forecast:::thetaf(y = time_series, h = forecast_horizon)$mean
    , error = function(e) {  
      warning(e)
      get_snaive_forecasts(time_series, forecast_horizon)
    })
}
# Calculate snaive forecasts
get_snaive_forecasts <- function(time_series, forecast_horizon){
  forecast:::snaive(time_series, h = forecast_horizon)$mean
}

# The following function computes the MASE for the prediction "forecast"
# when trying to predict "test_set". The model used to generate the
# prediction "forecast" has been trained with the time series "training_set".
# "f" is the frequency of the series
compute_mase <- function(forecast, test_set, training_set, f) {
  f <- switch(f,"yearly" = 1, "quarterly" = 4, "monthly" = 12, "weekly" = 52,
              "daily" = 7)
  mean(abs(forecast-test_set)) / mean(abs(diff(training_set, lag = f)))
}

# Function to read a M competition dataset of time series from a file with
# tsf format downloaded from https://forecastingdata.org/
# Returns a list in which every component is a list storing:
#    * a training set(a time series) 
#    * its associated test set (another time series),
#    * the name of the series (character)
#    * the frequency of the series (character)
# Example of use
# > ds <- load_dataset("m4_yearly_dataset.tsf")
# > ds[[1]]$training  # training set of the first series
# > ds[[1]]$test      # test set of the first series
# > ds[[1]]$frequency # frequency of the first series
# > library(vctsfr)
# > plot_ts(ds[[1]]$training, future = ds[[1]]$test) # get a plot
load_dataset <- function(file) {
  f <- file(file, "r")
  line <- readLines(f, n = 1)
  while (line != "@data") {
    if (substring(line, 1, 1) == "@") {
      l <- unlist(strsplit(line, split = " "))
      if (l[1] == "@horizon") {
        horizon <- as.numeric(l[2])
      }
      if (l[1] == "@frequency") {
        freq <-  l[2]
        frequency <- switch(l[2],
                            "yearly" = 1,
                            "quarterly" = 4,
                            "monthly" = 12,
                            "weekly" = 1,
                            "daily" = 7
        )
      }
    }
    line <- readLines(f, n = 1)
  }
  line <- readLines(f, n = 1)
  S <- list()
  while (length(line) == 1) {
    l <- unlist(strsplit(line, split = ":"))
    name <- l[1]
    date <- as.Date(l[2], format = "%Y-%m-%d %H-%M-%S")
    if (freq == "yearly") {
      year <- lubridate::year(date)
      start <-  year
    } else if (freq == "quarterly") {
      yq <- tsibble::yearquarter(date)
      year <- lubridate::year(yq)
      quarter <- lubridate::quarter(yq)
      start <-  c(year, quarter)
    } else if (freq == "monthly") {
      year <- lubridate::year(date)
      month <- lubridate::month(date)
      start <-  c(year, month)
    } else if (freq == "weekly") {
      start <- 1
    } else if (freq == "daily") {
      start <- c(1, 1)
    }
    timeS <- as.numeric(unlist(strsplit(l[3], split = ",")))
    timeS <- stats::ts(timeS, start = start, frequency = frequency)
    training <- stats::window(timeS, end = stats::time(timeS)[length(timeS)-horizon])
    test <- stats::window(timeS, start = stats::time(timeS)[length(timeS)-horizon+1])
    line <- readLines(f, n = 1)
    S[[length(S)+1]] <- list(name = name, frequency = freq, training = training, test = test)
  }
  close(f)
  S
}

# Code to compute the mean MASE of the Theta method used in the Monash Time
# Series forecasting Repository for the 23,000 series in the yearly dataset 
# of the M4 competition
# The expected mean MASE is 3.375
# The computation is quite fast, it took about 1 minute and 36 seconds in my computer
file <- "G:/Mi unidad/data/m4_yearly_dataset.tsf" # use your tsf file
# Download the M4 yearly dataset in tsf format from https://forecastingdata.org/
ds <- load_dataset(file) # load the dataset

mase_theta <- vector("numeric", length = length(ds)) # MASE for every series
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  f <- get_theta_forecasts(training, length(test))
  mase_theta[ind] <- compute_mase(f, test, training, ds[[ind]]$frequency)
}
# The mean MASE
cat("Mean MASE of the Theta approach implemented in Monash Repository:", mean(mase_theta), "\n")

# Now we use the same steps to compute the mean MASE of the knn model with 
# k = 3, additive transformation and targets averaged with the mean function
# The computation is quite fast, it took 48 seconds in my computer
# The mean MASE obtained is 3.290

mase_knn <- vector("numeric", length = length(ds)) # MASE of every series
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  f <- tsfknn::knn_forecasting(training, h = length(test), k = 3, 
                               transform = "additive", cf = "mean")$prediction
  mase_knn[ind] <- compute_mase(f, test, training, ds[[ind]]$frequency)
}
# The mean MASE
cat("Mean MASE of KNN with k = 3, additive transformation and mean of targets:", mean(mase_knn), "\n")

# Again, we follow the same methodology to compute the mean MASE of the ensemble
# of 27 KNN models with different configurations for the M4 yearly dataset
# In my computer the execution took about 18 minutes
# The mean MASE obtained is 3.163
st <- Sys.time()
mase_ensemble <- vector("numeric", length = length(ds)) # MASE of every series
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  f <- vector("numeric", length = length(test))
  for (k in c(3, 5, 7)) {
    for (transformation in c("additive", "multiplicative", "none")) {
      for (combination in c("mean", "median", "weighted")) {
        f <- f + tsfknn::knn_forecasting(training, h = length(test), k = k, 
                                         transform = transformation, 
                                         cf = combination)$prediction
      }
    }
  }
  f <- f / 27
  mase_ensemble[ind] <- compute_mase(f, test, training, ds[[ind]]$frequency)
}
print(Sys.time()-st)
cat("Mean MASE of the ensemble:", mean(mase_ensemble), "\n")

# This time we follow the same methodology to compute the mean MASE of the model
# that best performs on the validation set on the M4 yearly dataset
# In my computer took about 1 hour and 9 minutes
# The mean MASE obtained is 3.805
st <- Sys.time()
mase_bmvs <- vector("numeric", length = length(ds)) # MASE for every series
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  minimum <- .Machine$double.xmax        # RMSE of the best model on the validation set
  minimum_f <- rep(NaN, length(test))    # forecast of the best model
  for (k in c(3, 5, 7)) {
    for (transformation in c("additive", "multiplicative", "none")) {
      for (combination in c("mean", "median", "weighted")) {
        # train the model
        m <- tsfknn::knn_forecasting(training, h = length(test), k = k, 
                                         transform = transformation, 
                                         cf = combination)
        # make the forecast
        f <- predict(m, h = length(test))$prediction
        # test its accuracy on the validation test
        ro <- tryCatch(tsfknn::rolling_origin(m, h = length(test)), condition = function(c) NULL)
        if (!is.null(ro) && ro$global_accu["RMSE"] < minimum) {
          minimum <- ro$global_accu["RMSE"]
          minimum_f <- f # update the forecast of the best model on the validation set
        }
      }
    }
  }
  mase_bmvs[ind] <- compute_mase(minimum_f, test, training, ds[[ind]]$frequency)
}
print(Sys.time()-st)
# For some short time series is impossible to assess the forecast accuracy of
# the model on the validation set, because the training set is too small to
# hold k training examples. For this time series the MASE will be NaN
cat("Number of series too short to use a validation set:", sum(is.nan(mase_bmvs)), "of", length(ds), "series\n")
cat("Mean MASE of the best model on the validation set:", mean(mase_bmvs, na.rm = TRUE), "\n")


# Now we reproduce the calculation of the mean MASE for the 4Theta method used
# in the M4 competition on the 23,000 series of the M4 yearly dataset
# The implementation is downloaded from: https://github.com/Mcompetitions/M4-methods/blob/master/4Theta%20method.R
# In my computer took about 4 hours and 9 minutes
# The mean MASE obtained is 3.185

#Authors: E. Spiliotis and V. Assimakopoulos (2017) / Forecasting & Strategy Unit - NTUA
#Method Description: Generalizing the Theta model for automatic forecasting
#Method Type: Statistical - Decomposition

SeasonalityTest <- function(input, ppy){
  #Used for determining whether the time series is seasonal
  tcrit <- 1.645
  if (length(input)<3*ppy){
    test_seasonal <- FALSE
  }else{
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- ( abs(xacf[ppy]) > clim[ppy] )
    
    if (is.na(test_seasonal)==TRUE){ test_seasonal <- FALSE }
  }
  
  return(test_seasonal)
}

Theta.fit <- function(input, fh, theta, curve, model, seasonality , plot=FALSE){
  #Used to fit a Theta model
  
  #Check if the inputs are valid
  if (theta<0){ theta <- 2  }
  if (fh<1){ fh <- 1  }
  #Estimate theta line weights
  outtest <- naive(input, h=fh)$mean
  if (theta==0){
    wses <- 0
  }else{
    wses <- (1/theta)
  }
  wlrl <- (1-wses)
  #Estimate seasonaly adjusted time series
  ppy <- frequency(input)
  if (seasonality=="N"){
    des_input <- input ; SIout <- rep(1, fh) ; SIin <- rep(1, length(input))
  }else if (seasonality=="A"){
    Dec <- decompose(input, type="additive")
    des_input <- input-Dec$seasonal
    SIin <- Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }else{
    Dec <- decompose(input, type="multiplicative")
    des_input <- input/Dec$seasonal
    SIin <- Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }
  
  #If negative values, force to linear model
  if (min(des_input)<=0){ curve <- "Lrl" ; model <- "A"  }
  #Estimate theta line zero
  observations <- length(des_input)
  xs <- c(1:observations)
  xf = xff <- c((observations+1):(observations+fh))
  dat=data.frame(des_input=des_input, xs=xs)
  newdf <- data.frame(xs = xff)
  
  if (curve=="Exp"){
    estimate <- lm(log(des_input)~xs)
    thetaline0In <- exp(predict(estimate))+input-input
    thetaline0Out <- exp(predict(estimate, newdf))+outtest-outtest
  }else{
    estimate <- lm(des_input ~ poly(xs, 1, raw=TRUE))
    thetaline0In <- predict(estimate)+des_input-des_input
    thetaline0Out <- predict(estimate, newdf)+outtest-outtest
  }
  
  #Estimete Theta line (theta)
  if (model=="A"){
    thetalineT <- theta*des_input+(1-theta)*thetaline0In
  }else if ((model=="M")&(all(thetaline0In>0)==T)&(all(thetaline0Out>0)==T)){
    thetalineT <- (des_input^theta)*(thetaline0In^(1-theta))
  }else{
    model<-"A"
    thetalineT <- theta*des_input+(1-theta)*thetaline0In
  }
  
  #forecasting TL2
  sesmodel <- ses(thetalineT, h=fh)
  thetaline2In <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean
  
  #Theta forecasts
  if (model=="A"){
    forecastsIn <- as.numeric(thetaline2In*wses)+as.numeric(thetaline0In*wlrl)+des_input-des_input
    forecastsOut <- as.numeric(thetaline2Out*wses)+as.numeric(thetaline0Out*wlrl)+outtest-outtest
  }else if ((model=="M")&
            (all(thetaline2In>0)==T)&(all(thetaline2Out>0)==T)&
            (all(thetaline0In>0)==T)&(all(thetaline0Out>0)==T)){
    forecastsIn <- ((as.numeric(thetaline2In)^(1/theta))*(as.numeric(thetaline0In)^(1-(1/theta))))+des_input-des_input
    forecastsOut <- ((as.numeric(thetaline2Out)^(1/theta))*(as.numeric(thetaline0Out)^(1-(1/theta))))+outtest-outtest
  }else{
    model<-"A"
    thetalineT <- theta*des_input+(1-theta)*thetaline0In
    sesmodel <- ses(thetalineT,h=fh)
    thetaline2In <- sesmodel$fitted
    thetaline2Out <- sesmodel$mean
    forecastsIn <- as.numeric(thetaline2In*wses)+as.numeric(thetaline0In*wlrl)+des_input-des_input
    forecastsOut <- as.numeric(thetaline2Out*wses)+as.numeric(thetaline0Out*wlrl)+outtest-outtest
  }
  
  #Seasonal adjustments
  if (seasonality=="A"){
    forecastsIn <- forecastsIn+SIin
    forecastsOut <- forecastsOut+SIout
  }else{
    forecastsIn <- forecastsIn*SIin
    forecastsOut <- forecastsOut*SIout
  }
  
  #Zero forecasts become positive
  for (i in 1:length(forecastsOut)){
    if (forecastsOut[i]<0){ forecastsOut[i] <- 0 }
  }
  
  if (plot==TRUE){
    united <- cbind(input,forecastsOut)
    for (ik in 1:(observations+fh)){ united[ik,1] = sum(united[ik,2],united[ik,1], na.rm = TRUE) }
    plot(united[,1],col="black",type="l",main=paste("Model:",model,",Curve:",curve,",Theta:",theta),xlab="Time",ylab="Values",
         ylim=c(min(united[,1])*0.85,max(united[,1])*1.15))
    lines(forecastsIn, col="green") ; lines(forecastsOut, col="green")
    lines(thetaline2In, col="blue") ; lines(thetaline2Out, col="blue")
    lines(thetaline0In, col="red") ; lines(thetaline0Out, col="red")
  }
  
  output=list(fitted=forecastsIn,mean=forecastsOut,
              fitted0=thetaline0In,mean0=thetaline0Out,
              fitted2=thetaline2In,mean2=thetaline2Out,
              model=paste(seasonality,model,curve,c(round(theta,2))))
  
  return(output)
}

FourTheta<- function(input, fh){
  #Used to automatically select the best Theta model
  
  #Scale
  base <- mean(input) ; input <- input/base
  
  molist <- c("M","A") ; trlist <- c("Lrl","Exp")
  
  #Check seasonality & Create list of models
  ppy <- frequency(input) ; ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input, ppy) }
  if (ST==T){
    
    selist <- c("M","A")
    listnames <- c()
    for (i in 1:length(selist)){
      for (ii in 1:length(molist)){
        for (iii in 1:length(trlist)){
          listnames <- c(listnames,paste(selist[i], molist[ii], trlist[iii]))
        }
      }
    }
    
  }else{
    
    listnames <- c()
    for (ii in 1:length(molist)){
      for (iii in 1:length(trlist)){
        listnames <- c(listnames, paste("N", molist[ii], trlist[iii]))
      }
    }
    
  }
  
  modellist <- NULL
  for (i in 1:length(listnames)){
    modellist[length(modellist)+1] <- list(c(substr(listnames,1,1)[i], substr(listnames,3,3)[i],
                                             substr(listnames,5,7)[i]))
  }
  
  #Start validation
  errorsin <- c() ; models <- NULL
  
  #With this function determine opt theta per case
  optfun <- function(x, input, fh, curve, model, seasonality){
    mean(abs(Theta.fit(input=input, fh, theta=x, curve, model, seasonality , plot=FALSE)$fitted-input))
  }
  
  for (j in 1:length(listnames)){
    optTheta <- optimize(optfun, c(1:3),
                         input=input, fh=fh, curve=modellist[[j]][3], model=modellist[[j]][2],
                         seasonality=modellist[[j]][1])$minimum
    
    fortheta <- Theta.fit(input=input, fh=fh, theta=optTheta, curve=modellist[[j]][3], model=modellist[[j]][2],
                          seasonality=modellist[[j]][1], plot=F)
    models[length(models)+1] <- list(fortheta)
    errorsin <- c(errorsin, mean(abs(input-fortheta$fitted)))
  }
  
  #Select model and export
  selected.model <- models[[which.min(errorsin)]]
  description <- selected.model$model
  output <- list(fitted=selected.model$fitted*base,mean=selected.model$mean*base,
                 description=description)
  #Returns the fitted and forecasted values, as well as the model used (Type of seasonality, Type of Model, Type of Trend, Theta coef.)
  
  return(output)
  
}

st <- Sys.time()
mase_4theta <- vector("numeric", length = length(ds)) # mase of every forecast
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  f <- FourTheta(training, fh = length(test))$mean
  mase_4theta[ind] <- compute_mase(f, test, training, ds[[ind]]$frequency)
}
print(Sys.time()-st)
# The mean MASE
cat("Mean MASE of the 4Theta method of the M4 competition:", mean(mase_4theta), "\n")
