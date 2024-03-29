% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WaveDuration.R
\name{WaveDuration}
\alias{WaveDuration}
\title{Heat and cold waves duration for n-dimensional arrays}
\usage{
WaveDuration(
  data,
  threshold,
  op = ">",
  spell.length = 6,
  by.seasons = TRUE,
  dates = NULL,
  calendar = NULL,
  ncores = NULL
)
}
\arguments{
\item{data}{A numeric n-dimensional array containing daily maximum or minimum 
temperature}

\item{threshold}{An array with the threshold to be/not to be reach, usually 
given by the 90th/10th percentiles for heat/cold waves computed with the 
\code{Threshold} function.}

\item{op}{A character ">" (by default) or ">=" for heat waves and "<" or "<=" 
for cold waves indicating the operator must be used  to compare data to 
threshold.}

\item{spell.length}{A number indicating the number of consecutive days with 
extreme temperature to be considered heat or cold wave.}

\item{by.seasons}{If TRUE (by default), the wave duration is computed for each 
season (DJF/MAM/JJA/SON) separately. If FALSE is specified, the monthly wave 
duration is computed.}

\item{dates}{A vector of dates including calendar attributes. If NULL (by 
default), the 'time' attributes of parameter 'data' is used.}

\item{calendar}{A character indicating the calendar type.}

\item{ncores}{The number of cores to be used when computing the wave duration.}
}
\value{
A list of length 2:
\itemize{
 \item{\code{$result}, an array with the same dimensions as the input 
       \code{data}, but with the time dimension reduce from daily to monthly 
       or seasonal resolution depending on the selected resolution in 
       \code{by.season}.}
 \item{\code{$years}, a vector of the years and season/months corresponding 
       to the resolution selected in \code{by.season} and temporal length of 
       the input \code{data}.}
}
}
\description{
This function computes the duration of a heat/cold wave as the number of 
consecutive days for which the maximum/minimum temperature is exceeding/below 
a threshold over a minimum number of days in month or seasonal resolution.
}
\examples{
##Example synthetic data:
data <- 1:(2 * 3 * 31 * 5)
dim(data) <- c(lon = 2, lat = 3, time = 31, model = 5)
time <- as.POSIXct(paste(paste(1900, 1, 1:31, sep = "-"), paste(12, 0, 0.0, 
                  sep = ":")), tz = "CET")
metadata <- list(time = list(standard_name = 'time', long_name = 'time', 
                calendar = 'standard', 
                units = 'days since 1970-01-01 00:00:00', prec = 'double', 
                dim = list(list(name ='time', unlim = FALSE))))
attr(time, "variables") <- metadata                 
attr(data, 'Variables')$dat1$time <- time
threshold <- rep(40, 31)

a <- WaveDuration(data, threshold, op = ">", spell.length = 6, 
                 by.seasons = TRUE, ncores = NULL)
str(a)
}
