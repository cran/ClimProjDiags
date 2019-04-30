#'Calculate spatial area-weighted average of multidimensional arrays
#'
#'
#'This function computes a spatial area-weighted average of n-dimensional arrays being possible to select a region and to add a mask to be applied when computing the average.
#'
#'@param data An array with minimum two dimensions of latitude and longitude.
#'@param lon Numeric vector of longitude locations of the cell centers of the grid of \code{data}. This vector must be the same length as the longitude dimension in the parameter \code{data}.
#'@param lat Numeric vector of latitude locations of the cell centers of the grid of \code{data}. This vector must be the same length as the latitude dimension in the parameter \code{data}. 
#'@param region A vector of length four indicating the minimum longitude, the maximum longitude, the minimum latitude and the maximum latitude of the region to be averaged.
#'@param mask A matrix with the same spatial dimensions of \code{data}. It can contain either a) TRUE where the value at that position is to be accounted for and FALSE where not, or b) numeric values, where those greater or equal to 0.5 are to be accounted for, and those smaller are not. Attention: if the longitude and latitude dimensions of the data and mask coincide in length, the user must ensure the dimensions of the mask are in the same order as the dimensions in the array provided in the parameter \code{data}.
#'@param londim An integer number indicating the position of the longitude dimension in the \code{data} object.
#'@param latdim An integer number indicating the position of the latitude dimension in the \code{data} object.
#'
#'@return An array, matrix or vector containig the area-weighted average with the same dimensions as \code{data}, except for the spatial longitude and latitude dimensions, which disappear.
#'
#'@importFrom plyr aaply
#'@examples
#'##Example synthetic data 1:
#'data <- 1:(2 * 3 * 4 * 5)
#'dim(data) <- c(lon = 2, lat = 3, time = 4, model = 5)
#'lat <- c(1, 10, 20)
#'lon <- c(1, 10)
#'
#'a <- WeightedMean(data = data, lon = lon, lat = lat, region = NULL, 
#'                  mask = NULL, londim = 1, latdim = 2)
#'str(a)
#'
#'mask <- c(0, 1, 0, 1, 0, 1)
#'dim(mask) <- c(lon = 2, lat = 3)
#'a <- WeightedMean(data = data, lon = lon, lat = lat, region = NULL, 
#'                  mask = mask, londim = 1, latdim = 2)
#'str(a)
#'
#'region <- c(1, 10, 1, 10)
#'a <- WeightedMean(data = data, lon = lon, lat = lat, region = region, 
#'                  mask = mask, londim = 1, latdim = 2)
#'str(a)
#'
#'##Example synthetic data:
#'data <- 1:(2 * 3 * 4)
#'dim(data) <- c(lon = 2, lat = 3, time=4)
#'lat <- c(1, 10, 20)
#'lon <- c(1, 10)
#'
#'a <- WeightedMean(data = data, lon = lon, lat = lat, region = NULL,
#'                  mask = NULL, londim = 1, latdim = 2)
#'str(a)
#'@export
WeightedMean <- function(data, lon, lat, region = NULL, mask = NULL, londim = NULL, latdim = NULL){
  if (is.null(data) | is.null(lon) | is.null(lat)) {
    stop("Parameter 'data', 'lon' and 'lat' cannot be NULL.")
  }
  if (!is.numeric(data) | !is.numeric(lon)| !is.numeric(lat)) {
    stop("Parameter 'data', 'lon' and 'lat' must be a numeric.")
  }
  if (length(dim(data)) < 2) {
    stop("Parameter 'data' needs to have dimensions lon and lat.")
  }
  if (!is.null(dim(lat)) | !is.null(dim(lon))) {
    stop("Parameter 'lon' and lat' need to be a vector.")
  }
  dim_names <- names(dim(data))
  dims <- 1 : length(dim(data))
  if (is.null(londim)) {
    if (!is.null(dim_names)) {
      londim <- which(dim_names == 'lon')
    }
    if(is.null(londim)) {
      londim <- which(dim(data) == length(lon))
    }
    if (length(londim) == 0)  {
      stop("No longitudinal dimension provided in parameter 'londim' nor as attribute of parameter 'data'.")
    }
  }
  if (is.null(latdim)) {
    if (!is.null(dim_names)) {
      latdim <- which(dim_names == 'lat')
    }
    if (is.null(latdim)) {
      latdim <- which(dim(data) == length(lat))
    }  
    if (length(latdim) == 0) {
        stop("No latitudinal dimension provided in parameter 'latdim' nor as attribute of parameter 'data'.")
    }
  }
  if (londim == latdim) {
    stop("Parameter 'londim' and 'latdim' cannot be equal.")
  }
  if (dim(data)[londim] != length(lon)){
    stop("The longitudinal dimension of parameter 'data' must be the same length of parameter 'lon'.")
  }
  if (dim(data)[latdim] != length(lat)){
    stop("The latitudinal dimension of parameter 'data' must be the same length of parameter 'lat'.")
  }
  nlon <- length(lon)
  nlat <- length(lat)
  if (!is.null(region)) {
    aux <- SelBox(data, lon = lon, lat = lat, region = region, londim = londim, latdim = latdim, mask = mask)
    data <- aux$data
    lon <- aux$lon
    lat <- aux$lat
    if (!is.null(mask)) {
      mask <- aux$mask
    }
  }
  wtmean <- aaply(data, .margins = dims[c(-londim, -latdim)], 
                        .fun = .WeightedMean, lon, lat, mask, .drop = FALSE)
  dim(wtmean) <- dim(data)[-c(latdim,londim)]
  if(length(dim(data)) > 3) {
    if (is.null(dim_names)) {
      dim_names <- paste0("dim", 1:length(dim(data)))
    }
    names(dim(wtmean)) = dim_names[-c(londim, latdim)]
  } else {
    attributes(wtmean) <- NULL
  }
  wtmean
}
.WeightedMean <- function(data, lon, lat, mask) {
  cosphi <- t(array(cos(lat * pi / 180), dim = c(length(lat), length(lon))))
  nblat <- length(lat)
  nblon <- length(lon)
  dlon <- abs(c(lon[2 : nblon] - lon[1 : nblon - 1])) * pi / 180
  dlon <- c(dlon, dlon[1])
  dlon <- array(dlon, dim = c(nblon, nblat))
  dlat <- abs(c(lat[2 : nblat] - lat[1 : nblat - 1])) * pi / 180
  dlat <- c(dlat, dlat[1])
  dlat <- t(array(dlat, dim = c(nblat, nblon)))
  weight <- (dlon * dlat * cosphi) 
  if ((is.null(mask) == FALSE) && all(dim(data) != dim(mask))) {
    stop("The provided parameter 'mask' must have the same size as the parameter 'data'.")
  }
  if ((nblat == dim(data)[1]) && (nblon == dim(data)[2])) {
    data <- t(data)
    if (!is.null(mask)) {
      mask <- t(mask)
    }
  } 
  if ((nblon != dim(data)[1]) || (nblat != dim(data)[2])) {
    stop("The parameter 'data' needs to have at least two dimensions of 'lon' and 'lat'.")
  }
  if(!is.null(mask)) {
    data[mask < 0.5] <- NA
  }
  weight[is.na(data)] <- NA
  coeff <- sum(weight, na.rm = TRUE)
  mean <- sum(weight * data, na.rm = TRUE) / coeff
  output <- mean
}
