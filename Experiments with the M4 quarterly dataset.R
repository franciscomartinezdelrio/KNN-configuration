#----------------------------------------------------------------------------------
# This file contains the experimentation with the M4 quarterly dataset
# In my computer this experimentation took about 1 hour and 47 minutes
#----------------------------------------------------------------------------------

library(tsfknn)
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

# Load the yearly dataset
# You have to download the M4 competition quarterly dataset from: https://forecastingdata.org/
file <- "G:/Mi unidad/data/m4_quarterly_dataset.tsf" # use your tsf file
ds <- load_dataset(file) # load the dataset

k <- c(3, 5, 7)                              # number of neighbors
t <- c("additive", "multiplicative", "none") # transformation
cf <- c("mean", "median", "weighted")        # combination function

# Compute the method's abbreviated names
mn <- character() # method's names
for (k_ in k) {
  for (cf_ in cf) {
    cf_name <- switch(cf_, "mean" = "M", "median" = "m", "weighted" = "D")
    for (t_ in t) {
      t_name <- switch(t_, "additive" = "A", "multiplicative" = "M", "none" = "N")
      mn[length(mn)+1] <- sprintf("%d%s%s", k_, t_name, cf_name)
    }
  }
}
mn <- c(mn, c("ensemble", "bmvs")) # bmvs - best model on the validation set

st <- Sys.time()
mase <- matrix(0, nrow = length(ds), ncol = length(mn))
colnames(mase) <- mn
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  training <- ds[[ind]]$training  # training set
  test <-  ds[[ind]]$test         # test set
  h <- length(test)               # forecasting horizon 
  p <- numeric(h)                 # forecast of the ensemble
  minimum <- .Machine$double.xmax # RMSE of the best model on the validation set
  minimum_f <- rep(NaN, h)        # forecast of the best model
  for (ind_t in seq_along(t)) {
    m <- knn_forecasting(training, h = h, k = 3, transform = t[ind_t], cf = "mean")
    for (ind_k in seq_along(k)) {
      m$model$k <- k[ind_k]
      for (ind_cf in seq_along(cf)) {
        m$model$cf <- cf[ind_cf]
        f <- predict(m, h = h)$prediction
        ro <- tryCatch(tsfknn::rolling_origin(m, h = h), condition = function(c) NULL)
        if (!is.null(ro) && ro$global_accu["RMSE"] < minimum) {
          minimum <- ro$global_accu["RMSE"]
          minimum_f <- f # update the forecast of the best model on the validation set
        }
        p <- p + f
        t_name <- switch(t[ind_t], "additive" = "A", "multiplicative" = "M", "none" = "N")
        cf_name <- switch(cf[ind_cf], "mean" = "M", "median" = "m", "weighted" = "D")
        name <- sprintf("%d%s%s", k[ind_k], t_name, cf_name)
        mase[ind, name] <- compute_mase(f, test, training, ds[[ind]]$frequency)
      }
    }
  }
  p <- p / (length(mn)-2)
  mase[ind, "ensemble"]  <- compute_mase(p, test, training, ds[[ind]]$frequency)
  mase[ind, "bmvs"]      <- compute_mase(minimum_f, test, training, ds[[ind]]$frequency)
}
print(Sys.time()-st)
setwd("G:/Mi unidad/code/R/time series forecasting/Experiments ensembles of KNN")
save(mase, file = "results_M4_quarterly.RData")
# To compute the mean MASE across the dataset (a small number of series with bmvs can produce some errors)
round(apply(mase, 2, mean, na.rm = TRUE), 3)

# To make the statistical differences analysis
# library(scmamp)                              # devtools::install_github("b0rxa/scmamp") to install it
# load("results_M4_quarterly.RData")           # to execute if the data is saved in a file
# scmamp::imanDavenportTest(-mase)             # Iman Davenport's test
# plotCD(results.matrix = -mase, alpha = 0.05) # critical differences diagram
