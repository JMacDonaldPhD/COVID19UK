#' Data


importPopnData = function(filename, replaceNAwZero = F, ...){
  tmp = read.csv(filename, stringsAsFactors = FALSE, header = FALSE, ...)
  n = ncol(tmp) - 1
  Mat = data.matrix(tmp[,-(1:n)])
  rownames(Mat) = tmp[,2]
  if(isTRUE(replaceNAwZero)){
    Mat[is.na(Mat)] = 0
  }
  Mat
}


importUACaseData = function(filename, replaceNAwZero = F, ...){
  dat = read.csv(filename,
                 stringsAsFactors = F)
  cases = dat$TotalCases
  names(cases) = dat$GSS_NM
  cases
}



#' Population in English Unitary Authorities
#'
#' A Dataset containing the estimated population of Unitary Authority
#' Areas in England (mid-2017)
#'
#' @format A matrix with 149 named rows and 1 variable
#'
#' @source \url{https://www.ons.gov.uk/datasets/mid-year-pop-est/editions/time-series/versions/4}
#'
#'

#' Centroids of English Unitart Authorities
#'
#' A dataset containing the Eastling, Northling, Latitude and Longitude of the centroids
#' od Unitary Authority Areas in England.
#'
#' @format A matrix with 149 named rows and 4 variables
#'
#' @source
#'



