#' UA kernels


#' Adjacency Matrix (Bordering or not borderning)

#' Distance

#' Population (number of people travelling to and forth is dependent
#'             on 2 populations)

kernel.pop = function(N){
  scaledPop = N/1e5
  n = length(N)
  mat = matrix(1, nrow = n, ncol = n) - diag(1, n)
  for(i in 1:n){
    for(j in 1:n){
      mat[i, j] = (i != j)*scaledPop[i]*scaledPop[j]
    }
  }
  mat
}

#' @func
#' @param centroids Central point of each UA area (by average postcode)
kernel.dist = function(centroids){
  dist_matrix = sp::spDists(centroids)
  k = function(param){
    return(exp(-dist_matrix/param))
  }
}


