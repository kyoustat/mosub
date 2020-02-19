#' Predict Class Labels for Mixtures of Subspaces
#' 
#' Once the run from main function \code{msm} is done, \code{msmpred} predicts 
#' the class label for the given data. In this function, there is no need for an user 
#' to modify any output. It simply takes a returned object from \code{msm} function, 
#' which is a list of lists, and new data.
#' 
#' @param X an \eqn{(n\times p) data matrix.}
#' @param msmoutput output of \code{msm} function. 
#' @param print.progress a logical; \code{TRUE} to show completion of iterations by 10, \code{FALSE} otherwise.
#' 
#' @return a length-\eqn{n} vector of class labels.
#' 
#' @examples 
#' ## generate a toy example
#' set.seed(10)
#' alldat = gen.LP(n=500, K=2, iso.var=0.1)
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
#' ## run MSM algorithm with K=2
#' maxiter   = 10000
#' output2   = msm(train.dat, K=2, iter=maxiter)
#' test.pred = msmpred(test.dat, output2)
#' 
#' ## visualize with axis Y and Z for the test data
#' tX = test.dat[,2]
#' tY = test.dat[,3]
#' 
#' opar <- par(mfrow=c(1,2))
#' plot(tX,tY,col=test.lab, pch=19,cex=0.5,main="true label")
#' plot(tX,tY,col=test.pred,pch=19,cex=0.5,main="predicted label")
#' par(opar)
#' 
#' @export
msmpred = function(X, msmoutput, print.progress=TRUE){
  #Given the list of lists of P[[iter]][[k]] as the projection matrices
  #and theta[[iter]][[k]] as the list of thetas
  #assign observations X[n,m] to K clusters. 
  #While we are here, give a confidence interval for the distance of subspaces by calculating the conway distance at each iter
  
  # parameters
  n = nrow(X)
  m = ncol(X)
  iterations = length(msmoutput);
  K = length(msmoutput[[1]]$P)
  
  #store cluster assignments
  cluster.mat = matrix(0,nrow=iterations, ncol=n);
  
  #store distances for each iteration
  dist.mat= matrix(0,nrow=K,ncol=n);
  
  #store each pairwaise distance
  n.pairs = K*(K-1)/2
  subdists = matrix(0,nrow = iterations, ncol = n.pairs);
  # print(paste("Iteration ", 0," at ", Sys.time()))
  for(i in 1:iterations){
    if(i%%10==0){
      if (print.progress){
        print(paste('* msmpred : iteration ', i,'/',iterations,' complete at ', Sys.time(),sep=""))  
      }
    }
    subdist.col = 1 #getting pairwise entry for subdist.col 1 first, incrementing from there
    for(j in 1:(K-1)){
      Pij = msmoutput[[i]]$P[[j]]
      for(h in (j+1):K){
        Pih = msmoutput[[i]]$P[[h]]
        subdists[i,subdist.col] = distance(Pij,Pih,subspace=F)
        # print("distance complete...")
        subdist.col=subdist.col+1;
      }
    }
    #Now get the distance for each observation from the subspace
    # print("42 : start")
    for(k in 1:K){
      P.ik = msmoutput[[i]]$P[[k]]
      theta.ik = (msmoutput[[i]]$theta[[k]])
      dist.mat[k,] = gibbs.loss.prj(x=X,P=P.ik,mu=theta.ik,subspace=F)
    }
    # print("42 : finished")
    #now determine which cluster was closest to the point
    
    cluster.mat[i,]=apply(dist.mat,2,which.min)
  }
  
  # returned object
  # return(list("cluster"=cluster.mat,"sub.dist"=subdists))
  predicted = rep(0,n)
  for (i in 1:n){
    tgt    = as.vector(cluster.mat[,i])
    tbtgt  = table(tgt)
    rnames = rownames(tbtgt)
    predicted[i] = round(as.numeric(rnames[as.numeric(which.max(tbtgt))]))
  }
    # apply(apply(cluster.mat, 2, table),2,which.max)
  
  return(predicted)
}
