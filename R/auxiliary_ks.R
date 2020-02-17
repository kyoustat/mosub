#' @keywords internal
mygibbs.step.c <- function(wi){
  n = nrow(wi)
  k = ncol(wi)
  
  iteration = 100000
  for (i in 1:iteration){
    labels = rep(0,n)
    for (j in 1:n){
      labels[j] = sample(1:k, 1, prob = wi[j,])
    }
    if (length(unique(labels))==k){
      break
    } 
  }
  return(labels)
}
#' @keywords internal
mygibbs.step.b <- function(eik, kappa=1){
  n = nrow(eik)
  k = ncol(eik)
  wi = array(0,c(n,k))
  for (i in 1:n){
    tgt = exp(-kappa*as.vector(eik[i,]))
    wi[i,] = tgt/sum(tgt)
  }
  wi[is.na(wi)] = 1
  wi = wi/rowSums(wi)
  return(wi)
}