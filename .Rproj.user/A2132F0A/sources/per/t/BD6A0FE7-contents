#' Generate Line and Plane Example
#' 
#' This function generates a toy example of 'line and plane' data in \eqn{R^3} that 
#' data are generated from a mixture of a line (one-dimensional) and a plane (two-dimensional).
#' 
#' @param n the number of data points for each line and plane.
#' @param shape shape parameter inverse-gamma for sampling variance parameter.
#' @param scale scale parameter of inverse-gamma for sampling variance parameter.
#' @param sig2 variance parameter for Gaussian errors in ambient space.
#' 
#' @return a named list containing:\describe{
#' \item{data}{an \eqn{(2n\times 3)} data matrix.}
#' \item{class}{length-\eqn{2n} vector for class label.}
#' }
#' 
#' @examples 
#' ## test for visualization
#' set.seed(10)
#' tester = gen.LP(shape=2, scale=5, sig2=0.05)
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
#' @export
gen.LP <- function(n=100, shape=2, scale=5, sig2=0.05){
  
  my.alpha = as.double(shape)
  my.beta  = as.double(scale)
  
  # draw line
  U1  = rstiefel::rustiefel(3,1)
  mu1 = stats::rnorm(1)
  sig1inv = 1/(stats::rgamma(1, shape=my.alpha, scale=my.beta))
  theta1  = (diag(3)-(U1%*%t(U1)))%*%(rnorm(3))
  
  # draw plane
  U2 = rstiefel::rustiefel(3,2)
  mu2 = stats::rnorm(2)
  sig2inv = diag(1/(stats::rgamma(2, shape=my.alpha, scale=my.beta)))
  theta2  = (diag(3)-(U2%*%t(U2)))%*%(rnorm(3))
  
  
  # parameters for generation
  vec1 = as.vector(U1*mu1)   + as.vector(theta1)
  vec2 = as.vector(U2%*%mu2) + as.vector(theta2)
  
  sigma1 = (1/sig1inv - sig2)*(U1%*%t(U1)) + (sig2*diag(3))
  sigma2 = U2%*%(base::solve(sig2inv)-(sig2*diag(2)))%*%t(U2) + (sig2*diag(3))
  
  # real generation
  dat1 = mvtnorm::rmvnorm(round(n), mean=vec1, sigma=sigma1)
  dat2 = mvtnorm::rmvnorm(n, mean=vec2, sigma=sigma2)

  # return
  output = list()
  output$data  = rbind(dat1, dat2)
  output$class = c(rep(1,round(n)),rep(2,n))
  return(output)
}