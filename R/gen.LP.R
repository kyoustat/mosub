#' Generate Line and Plane Example
#' 
#' This function generates a toy example of 'line and plane' data in \eqn{R^3} that 
#' data are generated from a mixture of lines (one-dimensional) planes (two-dimensional).
#' When \code{K=2}, it gives one line and one plane. For \code{K>2}, it gives one line, 
#' one plane, and \code{K-2} componenets are randomly chosen to be either a line or a plane.
#' 
#' @param n the number of data points for each line and plane.
#' @param K the number of mixture components.
#' @param iso.var degree of isotropic variance.
#' 
#' @return a named list containing:\describe{
#' \item{data}{an \eqn{(2n\times 3)} data matrix.}
#' \item{class}{length-\eqn{2n} vector for class label.}
#' }
#' 
#' @examples 
#' ## test for visualization
#' set.seed(10)
#' tester = gen.LP(n=100, K=2, iso.var=0.1)
#' data   = tester$data
#' label  = tester$class
#' 
#' ## do PCA for data reduction
#' proj = base::eigen(stats::cov(data))$vectors[,1:2]
#' dat2 = data%*%proj
#' 
#' ## visualize
#' opar <- par(mfrow=c(2,2))
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.5,col=label,main="PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.5,col=label,main="Axis 1 vs 2")
#' plot(data[,1],data[,3],pch=19,cex=0.5,col=label,main="Axis 1 vs 3")
#' plot(data[,2],data[,3],pch=19,cex=0.5,col=label,main="Axis 2 vs 3")
#' par(opar)
#' 
#' \dontrun{
#' ## visualize in 3d
#' x11()
#' scatterplot3d::scatterplot3d(x=data, pch=19, cex.symbols=0.5, color=label)
#' }
#' 
#' @export
gen.LP <- function(n=100, K=2, iso.var=0.1){
  
  # parameters
  m=3
  K=round(K)
  n=round(n)
  isotropic.var=as.double(iso.var)
  
  k.true= c(1:2)
  if (K>2){
    k.true = c(k.true, sample(1:2, size=K-2, replace=TRUE))
  }
  k.true = sort(k.true, decreasing = F)
  
  
  
  #J.true=2
  #Generate the subspaces
  Utrue.list = list()
  for(k in 1:K){
    Utrue.list[[k]] = rstiefel::rustiefel(m=m,R=k.true[k])
  }
  NU.true = list()
  for(k in 1:K){
    NU.true[[k]]=MASS::Null(Utrue.list[[k]])
  }
  PNU.list = list()
  for(k in 1:K){
    PNU.list[[k]] = NU.true[[k]]%*%t(NU.true[[k]])
  }
  #generate means in subspace coordinates
  mutrue.list = list()
  for(k in 1:K){
    mutrue.list[[k]] = rnorm(k.true[k])
  }
  
  #generate the residual space noise level
  phitrue.list = rep(10,K)
  sigmatrue.list = rep(isotropic.var,K)
  #generate the subspace variances
  sigma0.true.list = list()
  for(k in 1:K){
    sigma0.true.list[[k]]=runif(k.true[k],isotropic.var,5.1)
  }
  Sigma0.true.list = list()
  for(k in 1:K){
    Sigma0.true.list[[k]]=diag(sigma0.true.list[[k]],k.true[k])
  }
  
  
  #generate the euclidean space coordinate mean vector theta
  theta.true.list = list()
  for(k in 1:K){
    theta.true.list[[k]] = PNU.list[[k]]%*%rnorm(m)
  }
  X = c()
  label = c()
  for(k in 1:K){
    X = rbind(X,MASS::mvrnorm(n,mu=Utrue.list[[k]]%*%mutrue.list[[k]]+theta.true.list[[k]],
                              Sigma = sigmatrue.list[k]^2*diag(m)/10+Utrue.list[[k]]%*%(Sigma0.true.list[[k]]-sigmatrue.list[k]*diag(k.true[k]))%*%t(Utrue.list[[k]])))
    label = c(label, rep(k,n))
  }
  # scatterplot3d(x=X[,1], y=X[,2], z=X[,3], pch = 19, color = label)
  

  # return
  output = list()
  output$data  = X
  output$class = label
  return(output)
}