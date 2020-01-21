#'Select apatial region from multidimensional arrays
#'
#'@description This function subsets an spatial region from spatial data giving a vector with the maximum and minimum of latitudes and longitudes of the selected region.
#'
#'@param data An array with minimum two dimensions of latitude and longitude.
#'@param lon Numeric vector of longitude locations of the cell centers of the grid of \code{data}'. 
#'@param lat Numeric vector of latitude locations of the cell centers of the grid of \code{data}'.
#'@param region A vector of length four indicating the minimum longitude, the maximum longitude, the minimum latitude and the maximum latitude.
#'@param londim An integer number indicating the position of the longitude dimension in the \code{data} object. If NULL (by deafault), the function search for a dimension call 'lon' in the \code{data} input.
#'@param latdim An integer number indicating the position of the latitude dimension in the \code{data} object.  If NULL (by deafault), the function search for a dimension call 'lat' in the \code{data} input.
#'@param mask A matrix with the same spatial dimensions of \code{data}.
#'
#'@return A list of length 4:
#'\itemize{
#'  \item\code{$data}{An array with the same dimensions as the input \code{data} array, but with spatial dimension reduced to the selected \code{region}}
#'  \item\code{$lat}{A vector with the new corresponding latitudes for the selected \code{region}}
#'  \item\code{$lon}{A vector with the new corresponding longitudes for the selected \code{region}}
#'  \item\code{$mask}{If parameter \code{mask} is supplied, an array  with reduced length of the dimensions to the selected \code{region}. Otherwise, a NULL element is returned.}}
#'
#'@examples 
#'## Example with synthetic data:
#'data <- 1:(20 * 3 * 2 * 4)
#'dim(data) <- c(lon = 20, lat = 3, time = 2, model = 4)
#'lon <- seq(2, 40, 2)
#'lat <- c(1, 5, 10)
#'
#'a <- SelBox(data = data, lon = lon, lat = lat, region = c(2, 20, 1, 5), 
#'            londim = 1, latdim = 2, mask = NULL)
#'str(a)
#'@export
SelBox <- function(data, lon, lat, region, londim = NULL, latdim = NULL, mask = NULL) {
  if (is.null(data) | is.null(lon) | is.null(lat) | is.null(region)){
    stop("Parameters 'data', 'lon', 'lat' or 'region' cannot be NULL.")
  }
  if (!is.numeric(data) | !is.numeric(lon) | !is.numeric(lat) | !is.numeric(region)){
    stop("Parameters 'data', 'lon', 'lat' or 'region' must be numeric.")
  }
  if (!is.array(data) && !is.matrix(data)) {
    stop("Parameter 'data' must be an array or matrix.")
  }
  if (!is.null(dim(lat)) | !is.null(dim(lon))) {
    stop("Parameter 'lon' and lat' need to be a vector.")
  }
  if (length(region) != 4){
    stop("The region argument has to be a vector of length four indicating the minimum longitude, the maximum longitude, the minimum latitude and the maximum latitude.")
  }
  dims <- 1:length(dim(data))
  if (is.null(londim)) {
    if ("lon" %in% names(dim(data))) {
      londim <- which(names(dim(data)) == "lon")
    } else if (length(lon) %in% dim(data)) {
      londim <- which(dim(data) == length(lon))
      if (length(londim) > 1) {
        stop("More than one dimension of the parameter 'data' has the same length as 'lon' parameter.")
      }
    } else {
      stop("Non of the dimensions of the parameter 'data' are of the same length as 'lon'.")
    }
  } 
  if (is.null(latdim)) {
    if ("lat" %in% names(dim(data))) {
      latdim <- which(names(dim(data)) == "lat")
    } else if (length(lat) %in% dim(data)) {
      latdim <- which(dim(data) == length(lat))
      if (length(latdim) > 1) {
        stop("More than one dimension of the parameter 'data' has the same length as 'lat' parameter.")
      }
    } else {
      stop("Non of the dimensions of the parameter 'data' are of the same length as 'lat' parameter.")
    }
  }
  if (londim == latdim) {
    stop("Parameter 'londim' and 'latdim' cannot be equal.")
  }
  if (dim(data)[londim] != length(lon)){
    stop("The longitudinal dimension of parameter 'data' must be of the same length of parameter  'lon'.")
  }
  if (dim(data)[latdim] != length(lat)){
    stop("The latitudinal dimension of parameter 'data' must be of the same length of parameter  'lat'.")
  }
  if (region[3] <= region[4]) {
    LatIdx <- which( lat >= region[3] & lat <= region[4])
  } else {
    LatIdx <- which(lat <= region[3] | lat >= region[4])
  }
  #if (region[1] <= region[2]) {
  #  LonIdx <- which(lon >= region[1] & lon <= region[2])
  #} else {
  #  LonIdx <- which(lon >= region[1] | lon <= region[2])
  #}
  LonIdx <- Lon2Index(lon, lonmin = region[1], lonmax = region[2])

  data <- Subset(data, along = londim, indices = LonIdx, drop = "none")
  data <- Subset(data, along = latdim, indices = LatIdx, drop = "none")
  if (!is.null(mask)) {
    mask <- Subset(mask, along = latdim - length(dim(data)) + 2, indices = LatIdx, drop = "none")
    mask <- Subset(mask, along = londim - length(dim(data)) + 2, indices = LonIdx, drop = "none")
  } else {
    mask <- NULL
  }
  list(data = data, lon = lon[LonIdx], lat = lat[LatIdx], mask = mask)
}
