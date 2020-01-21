#'Wrapper for applying the climdex routine  ETCCDI climate change indices to n-dimensional arrays.
#'
#'@description This function computes the t90p, t10p, cdd or rx5day indices from n-dimensional arrays.
#'
#'@param data A numeric n-dimensional array containing daily maximum or minimum temperature, wind speed or precipitation amount.
#'@param metric The metric to be computed, either 't90p', 't10p', 'Wx', 'cdd' or 'rx5day'.
#'@param threshold For the 't90p' and 't10p' metrics, an array of the 90th/10th percentiles must be included. This parameter can be computed with the \code{Threshold} function.
#'@param base.range The years used for the reference period. If NULL (by default), all years are used.
#'@param dates A vector of dates with a calendar attributes. If NULL (by default), the 'time' attributes of parameter 'data' are considered.
#'@param timedim An integer number indicating the position of the time dimension in the parameter \code{data}. If NULL (by default), the dimension called 'time' in parameter \code{data} is considered as temporal dimension.
#'@param calendar A character indicating the calendar type.
#'@param ncores The number of cores to be used when computing the index.
#'
#'@return A list of length 2:
#'\itemize{
#'  \item\code{$result} {An array with the same dimensions as the input array, except for the temporal dimension which is renamed to 'year', moved to the first dimension position and reduce to annual resolution.}
#'  \item\code{$years}  {A vector of the corresponding years.}}
#'
#'@import multiApply
#'@import PCICt
#'@examples 
#'##Example synthetic data:
#'data <- 1:(2 * 3 * 372 * 1)
#'dim(data) <- c(lon = 2, lat = 3, time = 372, model = 1)
#'time <- c(seq(ISOdate(1900, 1, 1), ISOdate(1900, 1, 31), "day"), 
#'          seq(ISOdate(1901, 1, 1), ISOdate(1901, 1, 31), "day"),
#'          seq(ISOdate(1902, 1, 1), ISOdate(1902, 1, 31), "day"),
#'          seq(ISOdate(1903, 1, 1), ISOdate(1903, 1, 31), "day"),
#'          seq(ISOdate(1904, 1, 1), ISOdate(1904, 1, 31), "day"),
#'          seq(ISOdate(1905, 1, 1), ISOdate(1905, 1, 31), "day"),
#'          seq(ISOdate(1906, 1, 1), ISOdate(1906, 1, 31), "day"),
#'          seq(ISOdate(1907, 1, 1), ISOdate(1907, 1, 31), "day"),
#'          seq(ISOdate(1908, 1, 1), ISOdate(1908, 1, 31), "day"),
#'          seq(ISOdate(1909, 1, 1), ISOdate(1909, 1, 31), "day"),
#'          seq(ISOdate(1910, 1, 1), ISOdate(1910, 1, 31), "day"),
#'          seq(ISOdate(1911, 1, 1), ISOdate(1911, 1, 31), "day"))
#'metadata <- list(time = list(standard_name = 'time', long_name = 'time',  calendar = 'gregorian', 
#'                             units = 'days since 1970-01-01 00:00:00', prec = 'double', 
#'                             dim = list(list(name = 'time', unlim = FALSE))))
#'attr(time, "variables") <- metadata
#'attr(data, 'Variables')$dat1$time <- time
#'
#'thres <- rep(10, 31 * 2 * 3)
#'dim(thres) <- c(jdays = 31, lon = 2, lat = 3,  model = 1)
#'str(thres)
#'
#'
#'clim <- Climdex(data, metric = "t90p", threshold = thres)
#'str(clim)
#'@references David Bronaugh for the Pacific Climate Impacts Consortium (2015).
#'  climdex.pcic: PCIC Implementation of Climdex Routines. R package
#'  version 1.1-6. http://CRAN.R-project.org/package=climdex.pcic
#'@export
Climdex <- function(data, metric, threshold = NULL, base.range = NULL, dates = NULL, timedim = NULL, 
                    calendar = NULL, ncores = NULL) {
  if (is.null(data) | is.null(metric)) {
    stop("Parameters 'data' and 'metric' cannot be NULL.")
  }
  if (!is.numeric(data)) {
    stop("Parameters 'data' must be numeric.")
  }
  if (is.null(dim(data))) {
    dim(data) = c(time = length(data))
    timedim = 1
  }
  if (is.null(timedim) & !is.null(names(dim(data)))) {
    timedim <- which(names(dim(data)) == "time")
  }
  if (is.null(timedim)) {
    stop("No time dimension provided in parameter 'timedim' nor as dimension names of parameter 'data'.")
  }
  if (is.null(dates)) {
    dates <- attr(data, 'Variables')$common$time
    if (is.null(dates)) {
      dates <- attr(data, 'Variables')$dat1$time
    }
  }
  if (is.null(dates)) {
    stop("No dates provided in parameter 'dates' nor as attribute of parameter 'data' or 'dates'.")
  }
  if (is.null(calendar)) {
    calendar <- attributes(dates)$calendar
    if (is.null(calendar)) {
      calendar <- attributes(dates)$variables$time$calendar
    }
    if (is.null(calendar)) {
      stop("The attribute 'calendar' must be present in the parameter 'dates' or specified in parameter 'calendar'.")
    }
  }  
  if (!any(class(dates) %in% 'POSIXct')) {
    dates <- try( {
      if (is.character(dates)) {
        as.POSIXct(dates, format = "%Y%m%d")
      } else {
        as.POSIXct(dates)
      }
    })
    if ('try-error' %in% class(dates) | sum(is.na(dates)) == length(dates)) {
      stop("Dates provided in parameter 'dates' or as attribute of parameter 'data' must be of class 'POSIXct' or convertable to 'POSIXct'.")
    }
  }
  stop_error <- FALSE
  if (length(dim(data)) == 1) {
    if (length(dates) != length(data)) {
      stop_error <- TRUE
    }
  } else {
    if (length(dates) != dim(data)[timedim]) {
      stop_error <- TRUE
    }
  }
  if (stop_error) {
    stop("Parameter 'dates' must be of the same length as the 'time' dimension of the parameter 'data'.")
  }
  dates <- as.PCICt(dates, cal = calendar)
  dates = as.character(dates)
  jdays <- as.numeric(strftime(dates, format = "%j"))
  if (calendar == "gregorian" | calendar == "standard" | calendar == "proleptic_gregorian") {
    year <- as.numeric(strftime(dates, format = "%Y"))
    if (length(unique(year)) > 1) {
      pos <- ((year / 100) %% 1 == 0) + ((year / 4) %% 1 == 0) + ((year / 400) %% 1 == 0)
      pos <- which(pos == 0 | pos == 2 | pos == 4)
      if (length(pos) > 0) {
        pos <- pos[which(jdays[pos] > 59)] 
        jdays[pos] <- jdays[pos] + 1
      }
    }
  }
  if (!is.character(metric) | (metric != "cdd" & metric != "t90p" & metric != "t10p" & metric != "rx5day" & 
                               metric != "Wx")) {
    stop("Parameter 'metric' must be a character indicating the metric to be computed, either 't90p', 't10p', 'Wx', 'cdd' or 'rx5day'.")
  }
  if (length(metric) > 1) {
    metric = metric[1]
    warning("Parameter 'metric' has length > 1 and only the first element will be used.")
  }
  if ((metric == "t90p" | metric == "t10p" | metric == "Wx") & is.null(threshold)) {
    stop("Parameter 'threshold' cannot be NULL for metric 't90p', 't10p' or 'Wx'.")
  }
  if ((metric == "cdd" | metric == "rx5day") & !is.null(threshold)) {
    threshold = NULL
    warning("Parameter 'threshold' haven't be used when computing metric 'cdd' or 'rx5day'.")
  }
  if ((metric == "t90p" | metric == "t10p" | metric == "Wx") & !is.null(threshold)) { 
    if (!is.null(names(dim(threshold)))) {
      time_dim_threshold <- which(names(dim(threshold)) == "jdays" | names(dim(threshold)) == "time") 
    } else {
      if (is.null(dim(threshold))) {
        time_dim_threshold <- 1
        dim(threshold) <- c(jdays = length(threshold))
      } else {
        stop("Parameter 'threshold' must have a dimension called 'jdays'.")
      }
    }
  }
  if (!is.null(threshold)) {
    if (!is.numeric(threshold)) {
      stop("Parameter 'threshold' must be numeric.")
    }
  }
  if (is.null(names(dim(data)))) {
    names(dim(data)) <- paste0("dims", 1 : length(dim(data)))
    names(dim(data))[timedim] <- "time"
  } 
  dims <- dim(data)
  if (metric == "t90p" | metric == "t10p" | metric == "Wx") {
    if (length(unique(jdays)) != dim(threshold)[time_dim_threshold]) {
      if (dim(threshold)[time_dim_threshold] != 1) {
        stop("Length of 'jdays' dimension provided in parameter 'threshold' must be the consistent with 'jdays' provided in parameter 'dates.")
      } else {
        threshold = rep(threshold, length(unique(jdays)))
        dim(threshold) = c(jdays = length(unique(jdays)))
        warning("Parameter 'threshold' has been recycled.")
      }
    }
    if (length(dim(threshold)) == 1) {
      if (length(dim(data)) > 1) {
         threshold <- rep(threshold, prod(dim(data)[-timedim]))
         dim(threshold) <- c(jdays = length(unique(jdays)), dim(data)[-timedim])
         names(dim(threshold)[-time_dim_threshold]) <- names(dim(data)[-timedim])
         warning("Parameter 'threshold' has been recycled.")
      }
    }
  }
  if (!is.numeric(ncores) & !is.null(ncores)) {
    stop("Parameter 'ncores' must be numeric.")
  }
  if (length(ncores) == 0 & !is.null(ncores)) {
    stop("Parameter 'ncores' must be of length 1.")
  }
  if (length(ncores) > 1) {
    ncores = ncores[1]
    warning("Parameter 'ncores' has length > 1 and only the first element will be used.")
  }
  if (!is.null(ncores)) {
    ncores <- round(ncores)
    if (ncores == 0) {
      ncores = NULL
    }
  }
  if (metric == "cdd" || metric == "rx5day") {
    dates <- as.factor(substr(dates, 1, 4))
    date.factor <- NULL
    jdays <- NULL
    years <- levels(dates)
    data <- list(data)
    target_dims <- list(timedim)
  } else if (metric == "t90p" || metric == "t10p" || metric == "Wx") {
    date.factor <- as.factor(substr(dates, 1, 4))
    years <- levels(date.factor)
    data <- list(data, threshold)
    target_dims <- list('time', 'jdays')
  }
  if (length(dims) > 1) {
    result <- Apply(data = data, fun = .Climdex,
                    # margins = list(2, 1),
                    target_dims = target_dims, 
                    output_dims = 'year',
                    dates = dates, date.factor = date.factor, metric = metric,
                    jdays = jdays, base.range = base.range, ncores = ncores)
  } else {
    result <- list()
    result$output1 <- .Climdex(data = data[[1]], threshold = data[[2]], dates = dates, metric = metric, 
                               jdays = jdays, date.factor = date.factor, base.range = base.range)
  }
  return(list(result = result$output1, years = as.numeric(years)))
}
.Climdex <- function(data, threshold, dates = dates, metric = metric,
                     jdays = jdays, date.factor = NULL, base.range = NULL) {
   if (metric == "cdd") {
    result <- .spell.length.max(daily.prec = as.vector(data), date.factor = dates, 1, "<", FALSE)
  } else if (metric == "rx5day") {
    result <- .nday.consec.prec.max(daily.prec = as.vector(data), date.factor = dates, 
                                   ndays = 5, center.mean.on.last.day = FALSE)
  } else if (metric == "t90p" || metric == "Wx") {
    result <- .percent.days.op.threshold(temp = data, dates = dates, jdays = jdays, 
                                        date.factor = date.factor, threshold.outside.base = threshold, 
                                        base.thresholds = threshold, base.range = base.range, 
                                        op = ">", 20)
  } else if (metric == "t10p") {
    result <- .percent.days.op.threshold(temp = data, dates = dates, jdays = jdays, 
                                        date.factor = date.factor, threshold.outside.base = threshold, 
                                        base.thresholds = threshold, base.range = base.range, 
                                        op = "<", 20)
  }
  return(result)
}

#' @title Maximum spell length
#' 
#' @description
#' This function returns the longest string of days which exceed or are below
#' the given threshold.
#' 
#' @details
#' This routine compares data to the threshold using the given operator,
#' generating a series of TRUE or FALSE values. It then computes the lengths of
#' sequences of TRUE values (spells) and chooses the longest spell in each
#' period (as defined by date.factor).
#' 
#' The \code{spells.can.span.years} option controls whether spells must always
#' terminate at the end of a period, or whether they may continue until the
#' criteria ceases to be met or the end of the data is reached. The default for
#' fclimdex is TRUE.
#' 
#' @param daily.prec Data to compute index on.
#' @param date.factor Date factor to split by.
#' @param threshold The threshold to compare to.
#' @param op The operator to use to compare data to threshold.
#' @param spells.can.span.years Whether spells can span years.
#' @return A timeseries of maximum spell lengths for each period.
#' @seealso \code{\link{climdex.cdd}}.
#' @keywords ts climate
#' @examples
#' 
#' prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
#' phony.date.factor <- factor(rep(1:2, each=5))
#' 
#' ## With spells spanning years...
#' cwd <- spell.length.max(prec.dat, phony.date.factor, 1, ">=", TRUE)
#' 
#' ## Without spells spanning years...
#' altcwd <- spell.length.max(prec.dat, phony.date.factor, 1, ">=", FALSE)
#' 
#' @noRd
.spell.length.max <- function(daily.prec, date.factor, threshold, op, spells.can.span.years) {
  bools <- match.fun(op)(daily.prec, threshold)

  if(spells.can.span.years) {
    all.true <- .tapply.fast(bools, date.factor, all)
    max.spell <- .tapply.fast(.get.series.lengths.at.ends(bools), date.factor, max)
    
    ## Mask out values which are in the middle of a spell with NA
    na.mask <- c(1, NA)[as.integer((max.spell == 0) & all.true) + 1]
    return(max.spell * na.mask)
  } else {
    return(.tapply.fast(bools, date.factor, function(x) { max(.get.series.lengths.at.ends(x)) }))
  }
}
#' Get series length at ends
#' 
#' This function takes a series of boolean values and returns a list of
#' integers of the same length corresponding to the lengths at the ends of
#' sequences of TRUE values.
#' 
#' It can often be useful to know how long a series of boolean values is. This
#' function provides a method of knowing where and how long such sequences are.
#' 
#' @param x Sequence of booleans.
#' @param na.value Value to replace NAs with.
#' @return A vector consisting of the lengths of sequences of TRUE values at
#' the location of the last TRUE value in the sequence, and zeroes elsewhere.
#' @keywords ts climate
#' @examples
#' 
#' ## Get lengths of sequences of TRUE values in a sequence
#' series.lengths <- get.series.lengths.at.ends(c(TRUE, TRUE, TRUE, FALSE,
#' TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))
#' 
#' 
#' @noRd
.get.series.lengths.at.ends <- function(x, na.value=FALSE) {
  stopifnot(is.logical(x) && is.logical(na.value))
  n <- length(x)
  if(n == 1)
    return(as.numeric(x))

  res <- rep(0, n)
  x[is.na(x)] <- na.value

  ## Compare series to lag-1 and lag+1 series; false added to trigger state transition from TRUE at ends of series
  start <- which(x & !(c(FALSE, x[1:(n - 1)])))
  end <- which(x & !(c(x[2:n], FALSE)))
  res[end] <- end - start + 1
  return(res)
}
#' @title Number of days (less than, greater than, etc) a threshold
#' 
#' @description
#' Produces sums of values that exceed (or are below) the specified threshold.
#' 
#' @details
#' This function takes a data series, the number of days in the running window,
#' a date factor to aggregate by, and an optional modifier parameter
#' (center.mean.on.last.day). It computes the n-day running sum of
#' precipitation and returns the maximum n-day total precipitation per unit
#' time, as defined by \code{date.factor}.
#' 
#' @param daily.prec Daily timeseries of precipitation.
#' @param date.factor Factor to aggregate by.
#' @param ndays Number of days in the running window.
#' @param center.mean.on.last.day Whether to center the n-day running mean on
#' the last day of the series, instead of the middle day.
#' @return A vector consisting of the maximum n-day sum of precipitation per
#' time interval.
#' @keywords ts climate
#' @examples
#' library(PCICt)
#' 
#' ## Parse the dates into PCICt.
#' tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' 
#' ## Load the data in.
#' ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP,
#' ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
#' tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))
#' 
#' ## Compute rx5day on a monthly basis.
#' rx5day <- nday.consec.prec.max(ci@@data$prec, ci@@date.factors$monthly, 5)
#' 
#' @noRd
.nday.consec.prec.max <- function(daily.prec, date.factor, ndays, center.mean.on.last.day=FALSE) {
  if(ndays == 1) {
    return(suppressWarnings(.tapply.fast(daily.prec, date.factor, max, na.rm=TRUE)))
  }
  ## Ends of the data will be de-emphasized (padded with zero precip data); NAs replaced with 0
  daily.prec[is.na(daily.prec)] <- 0
  prec.runsum <- .running.mean(daily.prec, ndays)
  prec.runsum[is.na(prec.runsum)] <- 0
  if(center.mean.on.last.day) {
      k2 = ndays %/% 2
      prec.runsum <- c(rep(0, k2), prec.runsum[1:(length(prec.runsum) - k2)])
  }
  return(.tapply.fast(prec.runsum, date.factor, max) * ndays)
}
#' @title Running Mean of a Vector
#'
#' @description Calculates the running means of a vector with a shifting window
#'
#' @details Returns a new vector the same length as vec, where the ith
#' element is the mean of the bin of elements centered at the ith element
#' of the original vector. Means cannot be calculated for elements less
#' than half the width of the bin from the beginning or end of the vector;
#' the result vector has NA in those positions.
#'
#' @param vec A vector
#' @param bin The number of entries to average over for each mean
#'
#' @return a vector containing the running mean of bin elements of vec
#'
#' @example
#' \dontrun { 
#' running.mean(c(1, 2, 3, 4, 5, 6), 2) 
#' }
#' \dontrun { 
#' running.mean(c(5, 5, 5, 5, 5), 4) 
#' }
#'@noRd
.running.mean <- function(vec, bin){
  vec = as.vector(vec)
  len = length(vec)
  if (bin<=1) {
    return (vec)
  }
  if (bin > len) {
    bin = len
  }
  left.bin = bin%/%2

  means = double(len)

  right.bin = bin - left.bin - 1
  means = c( sum(vec[1:bin]), diff(vec,bin) ) # find the first sum and the differences from it
  means = cumsum(means)/bin                  # apply precomputed differences
  means = c(rep(NA,left.bin), means, rep(NA,right.bin))   # extend to original vector length
  return(means)
}
#' Lengths of strings of TRUE values
#' 
#' Computes fraction of days above or below the baseline threshold for each
#' day, and averages them using the date factor passed in.
#' 
#' This function computes fractions of days above or below baseline thresholds
#' for each day, then aggregates them using \code{date.factor}. It is used to
#' implement TN/TX 10/90p.
#' 
#' @param temp Sequence of temperature values.
#' @param dates Sequence of associated dates.
#' @param jdays Sequence of associated days of year.
#' @param date.factor Factor to aggregate data using.
#' @param threshold.outside.base Sequence of thresholds to be used for data
#' outside the base period.
#' @param base.thresholds Data structure containing sets of thresholds to be
#' used inside the base period; see \link{climdexInput-class}.
#' @param base.range Date range (type PCICt) of the baseline period.
#' @param op Comparison operator to use.
#' @param max.missing.days Maximum number of NA values per time period.
#' @return A vector consisting of the mean fraction of days above or below the
#' supplied set of thresholds.
#' @note If date.factor is omitted, daily series will be returned.
#' @seealso \link{climdexInput-class}.
#' @keywords ts climate
#' @examples
#' library(PCICt)
#' 
#' ## Parse the dates into PCICt.
#' tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' 
#' ## Load the data in.
#' ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP,
#' ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
#' tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))
#' 
#' ## Compute monthly tx90p.
#' tx90p <- percent.days.op.threshold(ci@@data$tmax, ci@@dates, ci@@jdays,
#'                                    ci@@date.factors$monthly,
#'                                    ci@@quantiles$tmax$outbase$q90,
#'                                    ci@@quantiles$tmax$inbase$q90,
#'                                    ci@@base.range, ">",
#'                                    ci@@max.missing.days['monthly']) *
#'          ci@@namasks$monthly$tmax
#' 
#' @noRd
.percent.days.op.threshold <- function(temp, dates, jdays, date.factor, threshold.outside.base, base.thresholds, base.range, op='<', max.missing.days) {
  f <- match.fun(op)
  dat <- f(temp, threshold.outside.base[jdays])
  
  inset <- dates >= base.range[1] & dates <= base.range[2]
  ## Don't use in-base thresholds with data shorter than two years; no years to replace with.
  if(sum(inset) > 0 && length(dates) >= 360 * 2) {
    jdays.base <- jdays[inset]
    years.base <- .get.years(dates[inset])

    ## Get number of base years, subset temp data to base period only.
    temp.base <- temp[inset]
    years.base.range <- range(years.base)
    byrs <- (years.base.range[2] - years.base.range[1] + 1)

    ## Linearize thresholds, then compare them to the temperatures
    bdim <- dim(base.thresholds)
    dim(base.thresholds) <- c(bdim[1] * bdim[2], bdim[3])
    yday.byr.indices <- jdays.base + (years.base - .get.years(base.range)[1]) * bdim[1]
    f.result <- f(rep(temp.base, byrs - 1), base.thresholds[yday.byr.indices,])
    dim(f.result) <- c(length(yday.byr.indices), bdim[3])

    ## Chop up data along the 2nd dim into a list; sum elements of the list
    dat[inset] <- rowSums(f.result, na.rm=TRUE) / (byrs - 1)
  }
  dat[is.nan(dat)] <- NA
  if(missing(date.factor))
    return(dat)
  na.mask <- .get.na.mask(dat, date.factor, max.missing.days)
  ## FIXME: Need to monthly-ize the NA mask calculation, which will be ugly.
  ret <- .tapply.fast(dat, date.factor, mean, na.rm=TRUE) * 100 * na.mask
  ret[is.nan(ret)] <- NA
  return(ret)
}
## Get year
.get.years <- function(dates) {
  return(as.POSIXlt(dates)$year + 1900)
}

