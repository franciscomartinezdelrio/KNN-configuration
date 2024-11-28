#----------------------------------------------------------------------------------
# This file allows you to visually compare the forecast with multiplicative and
# additive transformations for the M4 yearly dataset.
# Instructions: You can source the file and the, with CTRL-L execute the last
# command in the file: GUI_collection(collection)
# A shiny app will allow you to visually compare the forecasts for the 23,000
# series
#----------------------------------------------------------------------------------

library(tsfknn)
library(tsibble)
library(lubridate)
library(vctsfr)

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

file <- "G:/Mi unidad/data/m4_yearly_dataset.tsf" # use your tsf file
ds <- load_dataset(file) # load the dataset

collection <- vector("list", length = length(ds))
for (ind in seq_along(ds)) {
  if (ind %% 500 == 0) cat(ind, "/", length(ds), "\n", sep = "")
  t <- ds[[ind]]$training                            # time series
  name <- ds[[ind]]$name                        # time series's name
  test <- ds[[ind]]$test                           # "future" values
  af <- tsfknn::knn_forecasting(t, h = length(test), k = 5, transform = "additive")$prediction
  mf <- tsfknn::knn_forecasting(t, h = length(test), k = 5, transform = "multiplicative")$prediction
  collection[[ind]] <- ts_info(t, future = test,
                               prediction_info("add", af),
                               prediction_info("mul", mf),
                               name = name
  )
}
GUI_collection(collection)
