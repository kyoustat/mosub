#' Predict Class Labels for k-Plane Clustering
#' 
#' Once the run from main function \code{kplane} is done, \code{kplanepred} predicts 
#' the class label for the given data. In this function, there is no need for an user 
#' to modify any output. It simply takes a returned object from \code{kplane} function, 
#' which is a list of lists, and new data.
#' 
#' @param X an \eqn{(n\times p) data matrix.}
#' @param kpoutput output of \code{kplane} function. 
#' 
#' @return  a length-\eqn{n} vector of class labels.
#' 
#' @examples 
#' ## generate a toy example of two plane components.
#' set.seed(8128)
#' alldat = gen.LP2(n=500, nl=0, np=2, iso.var=0.25)
#' 
#' ## separate as train/test data
#' id.train = sample(1:1000, 800, replace=FALSE)
#' id.test  = setdiff(1:1000, id.train)
#' 
#' train.dat = alldat$data[id.train,]
#' train.lab = alldat$class[id.train]
#' 
#' test.dat = alldat$data[id.test,]
#' test.lab = alldat$class[id.test]
#' 
#' ## run k-plane algorithm with K=2
#' output2   = kplane(train.dat, K=2)
#' test.pred = kplanepred(test.dat, output2)
#' 
#' ## visualize with axis Y and Z for the test data
#' tX = test.dat[,2]
#' tY = test.dat[,3]
#' 
#' opar <- par(mfrow=c(1,2))
#' plot(tX,tY,col=3-test.lab,pch=19,cex=0.5,main="true label")
#' plot(tX,tY,col=test.pred, pch=19,cex=0.5,main="predicted label")
#' par(opar)
#' 
#' 
#' @export
kplanepred <- function(X, kpoutput){
  ## simply apply cluster assignment
  K = ncol(kpoutput$w)
  pred.assign = kplane.assign.naive(X, K, kpoutput$w, kpoutput$gamma)
  return(pred.assign$cluster)
}
