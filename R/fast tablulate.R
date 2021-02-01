#' Fast Tabulate

fast.tabulate <- function(draws, groups){
  a <- data.table::data.table(V1 = groups, V2 = draws)
  b <- a[, .N, by = list(V1, V2)]
  c <- tapply(b$N, list(b$V1, b$V2), function(x) sum(x))
  c[] <- sapply(c, function(x) replace(x, is.na(x), 0))
  return(c)
}
