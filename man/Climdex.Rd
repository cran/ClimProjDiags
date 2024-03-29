% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Climdex.R
\name{Climdex}
\alias{Climdex}
\title{Wrapper for applying the climdex routine  ETCCDI climate change indices to 
n-dimensional arrays.}
\usage{
Climdex(
  data,
  metric,
  threshold = NULL,
  base.range = NULL,
  dates = NULL,
  timedim = NULL,
  calendar = NULL,
  ncores = NULL
)
}
\arguments{
\item{data}{A numeric n-dimensional array containing daily maximum or minimum 
temperature, wind speed or precipitation amount.}

\item{metric}{The metric to be computed, either 't90p', 't10p', 'Wx', 'cdd' or 
'rx5day'.}

\item{threshold}{For the 't90p' and 't10p' metrics, an array of the 90th/10th 
percentiles must be included. This parameter can be computed with the 
\code{Threshold} function.}

\item{base.range}{The years used for the reference period. If NULL 
(by default), all years are used.}

\item{dates}{A vector of dates with a calendar attributes. If NULL 
(by default), the 'time' attributes of parameter 'data' are considered.}

\item{timedim}{An integer number indicating the position of the time dimension 
in the parameter \code{data}. If NULL (by default), the dimension called 
'time' in parameter \code{data} is considered as temporal dimension.}

\item{calendar}{A character indicating the calendar type.}

\item{ncores}{The number of cores to be used when computing the index.}
}
\value{
A list of length 2:
\itemize{
 \item{\code{$result}, an array with the same dimensions as the input array, 
       except for the temporal dimension which is renamed to 'year', moved 
       to the first dimension position and reduce to annual resolution.}
 \item{\code{$years}, a vector of the corresponding years.}
}
}
\description{
This function computes the t90p, t10p, cdd or rx5day indices from 
n-dimensional arrays.
}
\examples{
##Example synthetic data:
data <- 1:(2 * 3 * 372 * 1)
dim(data) <- c(lon = 2, lat = 3, time = 372, model = 1)
time <- c(seq(ISOdate(1900, 1, 1), ISOdate(1900, 1, 31), "day"), 
         seq(ISOdate(1901, 1, 1), ISOdate(1901, 1, 31), "day"),
         seq(ISOdate(1902, 1, 1), ISOdate(1902, 1, 31), "day"),
         seq(ISOdate(1903, 1, 1), ISOdate(1903, 1, 31), "day"),
         seq(ISOdate(1904, 1, 1), ISOdate(1904, 1, 31), "day"),
         seq(ISOdate(1905, 1, 1), ISOdate(1905, 1, 31), "day"),
         seq(ISOdate(1906, 1, 1), ISOdate(1906, 1, 31), "day"),
         seq(ISOdate(1907, 1, 1), ISOdate(1907, 1, 31), "day"),
         seq(ISOdate(1908, 1, 1), ISOdate(1908, 1, 31), "day"),
         seq(ISOdate(1909, 1, 1), ISOdate(1909, 1, 31), "day"),
         seq(ISOdate(1910, 1, 1), ISOdate(1910, 1, 31), "day"),
         seq(ISOdate(1911, 1, 1), ISOdate(1911, 1, 31), "day"))
metadata <- list(time = list(standard_name = 'time', long_name = 'time', 
                            calendar = 'gregorian', 
                            units = 'days since 1970-01-01 00:00:00', 
                            prec = 'double', 
                            dim = list(list(name = 'time', unlim = FALSE))))
attr(time, "variables") <- metadata
attr(data, 'Variables')$dat1$time <- time

thres <- rep(10, 31 * 2 * 3)
dim(thres) <- c(jdays = 31, lon = 2, lat = 3,  model = 1)
str(thres)

clim <- Climdex(data, metric = "t90p", threshold = thres)
str(clim)
}
\references{
David Bronaugh for the Pacific Climate Impacts Consortium (2015).
 climdex.pcic: PCIC Implementation of Climdex Routines. R package
 version 1.1-6. http://CRAN.R-project.org/package=climdex.pcic
}
