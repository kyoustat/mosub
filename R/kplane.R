#' k-Plane Clustering Algorithm
#' 
#' It is one of the first method proposed in 2000 for fitting \eqn{K} planes by EM-like algorithm 
#' that updates parametrization of a plane and assigning data points to the nearest plane in an alternating manner.
#' Empirically, its performance is not competitive against newer methods but we include this as a test benchmark 
#' as well as it is simple and fast enough. 
#' 
#' @param X an \eqn{(n\times p) data matrix.}
#' @param K the number of clusters.
#' @param iter the number of EM-type updating.
#' @param print.progress a logical; \code{TRUE} to show completion of iterations by 10, \code{FALSE} otherwise.
#' 
#' @return a list whose elements are also lists of following elements: \describe{
#' \item{cluster}{length-\code{n} vector of cluster label.}
#' \item{w}{an \eqn{(p\times K)} matrix of orthonormal columns.}
#' \item{gamma}{a length-\eqn{K} vector of projection values.}
#' }
#' 
#' @examples 
#' ## generate a toy example of two plane components
#' set.seed(18)
#' tester = gen.LP2(n=100, nl=0, np=2)
#' data   = tester$data
#' label  = tester$class
#' 
#' ## do PCA for data reduction
#' proj = base::eigen(stats::cov(data))$vectors[,1:2]
#' dat2 = data%*%proj
#' 
#' ## run k-plane algorithm with K=2, 3, and 4
#' kplane2 = kplane(data, K=2)
#' kplane3 = kplane(data, K=3)
#' kplane4 = kplane(data, K=4)
#' 
#' ## extract clustering
#' finc2 = kplane2$cluster
#' finc3 = kplane3$cluster
#' finc4 = kplane4$cluster
#' 
#' ## visualize
#' opar <- par(mfrow=c(3,4))
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc2+1,main="K=2:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc3+1,main="K=3:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc4+1,main="K=4:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(2,3)")
#' par(opar)
#' 
#' @references 
#' \insertRef{bradley_k-plane_2000}{mosub} 
#' 
#' @export
kplane <- function(X, K=2, iter=496, print.progress=TRUE){
  #########################################################
  # Initialization : use the notation from the paper
  m = nrow(X)
  n = ncol(X)
  K = round(K)
  
  kmeans.init   = stats::kmeans(X, center=K)$cluster
  kmeans.update = kplane.update(X, K, kmeans.init)
  
  w.old     = kmeans.update$w
  gamma.old = kmeans.update$gamma
  clust.old = kmeans.init
  
  
  
  #########################################################
  # Iterate
  for (it in 1:iter){
    # (a) Cluster Assignment
    clust.tmp = kplane.assign.naive(X, K, w.old, gamma.old)
    clust.new = clust.tmp$cluster
    if (length(unique(clust.new))<K){
      clust.new = kplane.assign.probabilistic(clust.tmp$distmat)
    }  
    if ((length(clust.new)==1)&&(clust.new==FALSE)){
      if (print.progress){
        print(paste('* kplane : iteration ', it,'/',iter,' terminated due to collapse at ', Sys.time(),sep=""))  
      }
      break
    }
    
    # (b) Cluster Update
    update.list = kplane.update(X, K, clust.new)
    w.new     = update.list$w
    gamma.new = update.list$gamma
    
    # (c) Update information
    inc.clust = as.double(mclustcomp::mclustcomp(clust.old, clust.new, type="jaccard")[2])
    
    clust.old = clust.new
    w.old     = w.new
    gamma.old = gamma.new
    
    if (it >= 5){
      if (inc.clust > 0.99){
        if (print.progress){
          print(paste('* kplane : iteration ', it,'/',iter,' terminated at ', Sys.time(),sep=""))
        }
        break
      }
    }
    if (print.progress){
      if(it%%10==0){
        print(paste('* kplane : iteration ', it,'/',iter,' complete at ', Sys.time(),sep=""))
      }
    }
  }
  
  #########################################################
  # Return output
  output = list()
  output$cluster = clust.old
  output$w = w.old
  output$gamma = gamma.old
  return(output)
}

#' @keywords internal
kplane.update <- function(X, K, cluster){
  # cluster information
  n = ncol(X)
  list.clust = list()
  for (k in 1:K){
    list.clust[[k]] = which(cluster==k)
  }
  
  # prepare
  w.mat = array(0,c(n,K))
  gamma.vec = rep(0,K)
  count.group = rep(0,K)
  
  # compute w.mat
  for (k in 1:K){
    idk = list.clust[[k]]
    ml = length(idk)
    count.group[k] = ml
    
    Al = X[idk,]
    if (length(idk)==1){
      Al = matrix(Al, nrow = 1)
    }
    Bl = t(Al)%*%(diag(ml)-outer(rep(1,ml),rep(1,ml))/ml)%*%Al
    wl = as.vector(eigen(Bl)$vectors[,nrow(Bl)])
    
    w.mat[,k] = wl
    gamma.vec[k] = sum(as.vector(array(1,c(1,ml))%*%Al)*wl)/ml
  }
  
  # return
  output = list()
  output$w = w.mat
  output$gamma = gamma.vec
  return(output)
}
#' @keywords internal
kplane.assign.naive <- function(X, K, w.mat, gamma.vec){
  m = nrow(X)
  distmat = array(0,c(m,K))
  for (k in 1:K){
    distmat[,k] = as.vector(X%*%w.mat[,k])-gamma.vec[k]
  }
  output = list()
  output$cluster = (apply(distmat, 1, which.min))
  output$distmat = distmat
  return(output)
}
#' @keywords internal
kplane.assign.probabilistic <- function(distmat){
  wi = exp(-distmat)
  wi = wi/rowSums(wi)
  wi[is.na(wi)] = 1
  return(kplane.mygibbs.step.c(wi))
}
#' @keywords internal
kplane.mygibbs.step.c <- function(wi){
  n = nrow(wi)
  k = ncol(wi)
  
  cflag = TRUE
  iter  = 1
  while (cflag){
    iter = iter + 1
    labels = rep(0,n)
    for (j in 1:n){
      labels[j] = sample(1:k, 1, prob = as.vector(wi[j,]))
    }
    tbl.lab = table(label)
    if (length(tbl.lab)==k){
      if (all(tbl.lab>1)){
        cflag = FALSE
      }
    }
    if (iter == 100){
      print("* kplane : updating the cluster assignment collapses. Return the last successful update.")
      return(FALSE)
    }
  }
  return(labels)
}
