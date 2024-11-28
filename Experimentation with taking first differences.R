#----------------------------------------------------------------------------------
# This file contains the last experimentation in the paper with
# taking first differences
#----------------------------------------------------------------------------------

library(tsfknn)
library(forecast)
library(tsibble)
library(lubridate)

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

preprocess_ts <- function(timeS, differences = 1) {
  prepro <- list()
  prepro$original <- timeS
  
  last_values <- utils::tail(timeS, differences)
  timeS <- diff(timeS, differences = differences)
  prepro$last_values <- last_values
  prepro$differences <- differences
  
  prepro$preprocessed <- timeS
  return(prepro)
}

unpreprocess_ts <- function(forecast, prepro) {
  temp <- stats::diffinv(forecast, differences = prepro$differences, xi = prepro$last_values)
  forecast <- subset(temp, start = 1 + prepro$differences)
  return(forecast)
}

# Load the yearly dataset
# You have to download the M4 competition yearly dataset from: https://forecastingdata.org/
file <- "G:/Mi unidad/data/m4_yearly_dataset.tsf" # use your tsf file
ds <- load_dataset(file) # load the dataset
differences <- sapply(ds, function(d) ndiffs(d$training))
ds <- ds[differences == 1]
mase <- matrix(0, nrow = length(ds), ncol = 4)
colnames(mase) <- c("additive", "multiplicative", "none", "fd")
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  h <- length(test)
  f1 <- knn_forecasting(training, h = h, k = 5, transform = "additive")$prediction
  mase[ind, "additive"] <- compute_mase(f1, test, training, ds[[ind]]$frequency)
  f2 <- knn_forecasting(training, h = h, k = 5, transform = "multiplicative")$prediction
  mase[ind, "multiplicative"] <- compute_mase(f2, test, training, ds[[ind]]$frequency)
  f3 <- knn_forecasting(training, h = h, k = 5, transform = "none")$prediction
  mase[ind, "none"] <- compute_mase(f3, test, training, ds[[ind]]$frequency)
  p <- preprocess_ts(training, differences = 1)
  f4 <- knn_forecasting(p$preprocessed, h = h, k = 5, transform = "none")$prediction
  mase[ind, "fd"] <- compute_mase(unpreprocess_ts(f4, p), test, training, ds[[ind]]$frequency)
}
round(apply(mase, 2, mean), 3)


# Load the quarterly dataset
# You have to download the M4 competition quarterly dataset from: https://forecastingdata.org/
file <- "G:/Mi unidad/data/m4_quarterly_dataset.tsf" # use your tsf file
ds <- load_dataset(file) # load the dataset
differences <- sapply(ds, function(d) ndiffs(d$training))
ds <- ds[differences == 1]
mase <- matrix(0, nrow = length(ds), ncol = 4)
colnames(mase) <- c("additive", "multiplicative", "none", "fd")
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  h <- length(test)
  f1 <- knn_forecasting(training, h = h, k = 5, transform = "additive")$prediction
  mase[ind, "additive"] <- compute_mase(f1, test, training, ds[[ind]]$frequency)
  f2 <- knn_forecasting(training, h = h, k = 5, transform = "multiplicative")$prediction
  mase[ind, "multiplicative"] <- compute_mase(f2, test, training, ds[[ind]]$frequency)
  f3 <- knn_forecasting(training, h = h, k = 5, transform = "none")$prediction
  mase[ind, "none"] <- compute_mase(f3, test, training, ds[[ind]]$frequency)
  p <- preprocess_ts(training, differences = 1)
  f4 <- knn_forecasting(p$preprocessed, h = h, k = 5, transform = "none")$prediction
  mase[ind, "fd"] <- compute_mase(unpreprocess_ts(f4, p), test, training, ds[[ind]]$frequency)
}
round(apply(mase, 2, mean), 3)

# Load the monthly dataset
# You have to download the M4 competition monthly dataset from: https://forecastingdata.org/
file <- "G:/Mi unidad/data/m4_monthly_dataset.tsf" # use your tsf file
ds <- load_dataset(file) # load the dataset
differences <- sapply(ds, function(d) ndiffs(d$training))
ds <- ds[differences == 1]
mase <- matrix(0, nrow = length(ds), ncol = 4)
colnames(mase) <- c("additive", "multiplicative", "none", "fd")
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training
  test <- ds[[ind]]$test
  h <- length(test)
  f1 <- knn_forecasting(training, h = h, k = 5, transform = "additive")$prediction
  mase[ind, "additive"] <- compute_mase(f1, test, training, ds[[ind]]$frequency)
  f2 <- knn_forecasting(training, h = h, k = 5, transform = "multiplicative")$prediction
  mase[ind, "multiplicative"] <- compute_mase(f2, test, training, ds[[ind]]$frequency)
  f3 <- knn_forecasting(training, h = h, k = 5, transform = "none")$prediction
  mase[ind, "none"] <- compute_mase(f3, test, training, ds[[ind]]$frequency)
  p <- preprocess_ts(training, differences = 1)
  f4 <- knn_forecasting(p$preprocessed, h = h, k = 5, transform = "none")$prediction
  mase[ind, "fd"] <- compute_mase(unpreprocess_ts(f4, p), test, training, ds[[ind]]$frequency)
}
round(apply(mase, 2, mean), 3)
