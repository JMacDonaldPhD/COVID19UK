
#' kernel COVID

kernelprelockdown = function(ageLimits){
  polymod = socialmixr::polymod
  K = socialmixr::contact_matrix(polymod, countries = "United Kingdom",
                                 age.limits = ageLimits)$matrix
  kernel = function(beta){
    K = beta*K
    K
  }
}

# Contact is reduced in general
# Need to think about social mixing for:

# key workers


# Household mixing. Children mix with brothers/sisters and a parents


kernelpostlockdown = function(ageLimits){
  polymod = socialmixr::polymod
  K = socialmixr::contact_matrix(polymod, countries = "United Kingdom",
                                 age.limits = ageLimits)$matrix
  n = length(ageLimits)
  kernel = function(beta){
    #l_matrix = matrix(l[1], nrow = n, ncol = n)
    #[, n] = l_matrix[n, ] = l[2]
    K = beta*K
    K
  }
}




kernelKW = function(ageLimits){
  polymod = socialmixr::polymod
  K = socialmixr::contact_matrix(polymod, countries = "United Kingdom",
                                 age.limits = ageLimits)$matrix
  noGroups <- 2*length(ageLimits)
  kernel = function(a, b){
    KMatrix <- matrix(nrow = noGroups, ncol = noGroups)
    NKWxNKW <- (a*a)*K
    KWxKW <- (b*b)*K
    KWxNKW <- (a*b)*K

    NKW_index <- 1:(noGroups/2)
    KW_index <- (noGroups/2 + 1):noGroups

    KMatrix[NKW_index, NKW_index] <- NKWxNKW
    KMatrix[NKW_index, KW_index] <- KWxNKW
    KMatrix[KW_index, NKW_index] <- KWxNKW
    KMatrix[KW_index, KW_index] <- KWxKW

    KMatrix
  }
}


