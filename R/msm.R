#' Mixtures of Varying-Dimensional Subspaces Model
#' 
#' \code{msm} is a Bayesian model inferring mixtures of subspaces that are of possibly different dimensions. 
#' For simplicity, this function returns only a handful of information that are most important in 
#' representing the mixture model, including projection, location, and hard assignment parameters. 
#' 
#' @param X an \eqn{(n\times p) data matrix.}
#' @param K the number of mixtures.
#' @param iter the number of MCMC runs.
#' @param prop.var proposal variance parameter.
#' @param temperature temperature value for Gibbs posterior.
#' @param print.progress a logical; \code{TRUE} to show completion of iterations by 10, \code{FALSE} otherwise.
#' 
#' @return a length-\code{iter} list whose elements are also lists of following elements: \describe{
#' \item{P}{length-\code{K} list of projection matrices.}
#' \item{U}{length-\code{K} list of orthonormal basis.}
#' \item{theta}{length-\code{K} list of center locations of each mixture.}
#' \item{cluster}{length-\code{n} vector of cluster label.}
#' }
#' 
#' @examples 
#' ## generate a toy example
#' set.seed(10)
#' tester = gen.LP(n=100, K=2, iso.var=0.1)
#' data   = tester$data
#' label  = tester$class
#' 
#' ## do PCA for data reduction
#' proj = base::eigen(stats::cov(data))$vectors[,1:2]
#' dat2 = data%*%proj
#' 
#' ## run MSM algorithm with K=2, 3, and 4
#' maxiter = 5000
#' output2 = msm(data, K=2, iter=maxiter)
#' output3 = msm(data, K=3, iter=maxiter)
#' output4 = msm(data, K=4, iter=maxiter)
#' 
#' ## extract final clustering information
#' finc2 = output2[[maxiter]]$cluster
#' finc3 = output3[[maxiter]]$cluster
#' finc4 = output4[[maxiter]]$cluster
#' 
#' ## visualize
#' opar <- par(mfrow=c(3,4))
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.5,col=finc2,main="K=2:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.5,col=finc2,main="K=2:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.5,col=finc2,main="K=2:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.5,col=finc2,main="K=2:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.5,col=finc3,main="K=3:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.5,col=finc3,main="K=3:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.5,col=finc3,main="K=3:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.5,col=finc3,main="K=3:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.5,col=finc4,main="K=4:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.5,col=finc4,main="K=4:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.5,col=finc4,main="K=4:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.5,col=finc4,main="K=4:Axis(2,3)")
#' par(opar)
#' 
#' @references 
#' \insertRef{thomas_learning_2015}{mosub}
#' 
#' @export
msm <- function(X, K=2, iter=496, prop.var = 1.0, temperature=1e-6, print.progress=TRUE){
  # MY ADDITION
  m = ncol(X)
  n = nrow(X)
  temperature = as.double(temperature)
  talk = print.progress
  report = 10
  sub.met = as.double(abs(prop.var))
  sub.met.sd = rep(sub.met, K)
  R = base::sample(1:(m-1), K, replace = TRUE)
  Urandom  = TRUE # local pca
  my.kappa = 10
  
  ############# copied code
  #Intialize the list of K subspaces
  #Intialize the list of K subspaces and projection matrices
  kmeans.init = kmeans(X,K)
  U.list = list()
  P.list = list()
  if (Urandom){
    for(k in 1:K){
      unow        = rustiefel(m=m, R = R[k])
      U.list[[k]] = unow
      P.list[[k]] = (unow%*%t(unow))
    }  
  } else {
    for (k in 1:K){
      idnow = which(kmeans.init$cluster==k)
      U.list[[k]] = base::eigen(stats::cov(X[idnow,]))$vectors[,(1:R[k])]
      P.list[[k]] = kisung.outer(U.list[[k]])
    }
  }
  #Initialize theta list
  theta.list = list()
  for(k in 1:K){
    theta.list[[k]] =Null(U.list[[k]])%*%t(Null(U.list[[k]]))%*% kmeans.init$centers[k,]
  } 
  
  #Initialize theta storage
  theta.mat = vector("list",iter)
  
  #Intialize subspace storage. Each subspace iteration will store U.mat[[iter]] = U.list
  #Allowing each subspace to be accessed as U.mat[[iter]][[k]]
  U.mat = vector("list",iter)
  P.mat = vector("list",iter)
  
  
  #Initialize the distance matrix to be stored every iteration
  distance = matrix(0,nrow=K,ncol = n)
  
  #Initialize storage for the latent normal and sphere walks 
  z.store = array(0,dim = c(iter,K,m*(m+1)/2))
  s.store = array(0,dim = c(iter, K, m*(m+1)/2))
  
  #initialize a random normal location
  z = matrix(0,K,m*(m+1)/2)
  s = matrix(0,K,m*(m+1)/2)
  for(k in 1:K){
    temp = conway.sphere.step(z=rep(0,m*(m+1)/2),sd=1,m=m)
    z[k,]=temp$z
    s[k,]=temp$s
  }
  
  #Initialize acceptance counter
  accept = rep(0,K)
  
  #Get loss for initial estimates
  curr.lossclus = fast.log.loss(x=X, P.list = P.list,mu = theta.list,temperature=temperature)
  curr.loss = curr.lossclus$loss
  curr.clus = curr.lossclus$clus
  
  #initialize and store the clustering
  lat.clus = which(rmultinom(n=n, size = 1, prob=rep(1/K,K))==1,arr.ind=T)[,1]
  lat.clus.mat = matrix(0,nrow=iter,ncol=n)
  pi.mat = matrix(0,n,K)
  
  #Intialize the multinomial weights
  pi.mat = matrix(1/K,nrow=n,ncol=K)
  dir.prior = rep(1/K,K)
  n.vec = rep(0,K)
  r0 = rep(1,K)
  
  #set up tuning parameter
  tune = 0 
  tune.accept = rep(0,K)
  record.clus = list()
  for(i in 1:iter){
    # print(paste('iter is ',i))
    tune = tune+1
    if(talk){
      if(i%%report==0){
        print(paste('* msm : iteration ', i,'/',iter,' complete at ', Sys.time(),sep=""))
        
      }
    }
    
    
    #For each component, generate a proposal subspace, and then
    #accept or reject it based on the gibbs posterior
    for(k in 1:K){
      #Get proposal projection matrix
      #print('get proposal')
      proposal = conway.step(z=z[k,],sd=sub.met.sd[k],m=m)
      #print('get Subspace')
      prop.sub = con2sub(P=unembed(proposal$s),return.proj = F)
      #restrict samples to m/2 ([m+1]/2 if m is odd ) to stay on lower half of sphere
      #print('Restrict')
      if(is.matrix(prop.sub)){
        if(dim(prop.sub)[2]>ceiling(m/2)){
          if(dim(prop.sub)[2]==m){
            #print('Proposal full dimension, replace with 1 dimension')
            prop.sub= rustiefel(R=1,m=m)
            proposal$z=embed(prop.sub,subspace=T)
          }else{
            #print('Proposal not full dimension')
            prop.sub = Null(prop.sub)
            proposal$z = embed(prop.sub,subspace=T)
          }
        }
      }
      # print('Get projection')
      prop.proj = prop.sub%*%t(prop.sub)
      
      #Set the proposal list
      prop.P.list = P.list
      prop.P.list[[k]] = prop.proj
      
      #Choose an appropriate theta in the null space
      prop.null = Null(prop.sub)
      prop.nullproj = prop.null%*%t(prop.null)
      prop.theta.list = theta.list
      prop.theta.list[[k]] = prop.nullproj%*%theta.list[[k]]
      for(l in 1:K){
        if(is.null(dim(prop.theta.list[[l]]))){
          prop.theta.list[[l]]=matrix(prop.theta.list[[l]],nrow=m,ncol=1)
        }
      }
      
      #Get the loss of the proposal
      prop.lossclus = fast.log.loss(x=X,P.list = prop.P.list, mu = prop.theta.list,temperature=temperature)
      prop.loss = prop.lossclus$loss
      #The coin toss
      toss = log(runif(1,0,1))
      diff.loss = prop.loss - curr.loss
      if(toss<diff.loss){
        accept[k] = accept[k] +1
        tune.accept[k] = tune.accept[k]+1
        #int.acc = int.acc +1
        P.list[[k]] = prop.proj
        U.list[[k]] = (prop.sub) #--------------------------------
        theta.list[[k]] = prop.theta.list[[k]]
        z[k,] =proposal$z 
        curr.loss = prop.loss
        curr.clus = prop.lossclus$clus
      }
    }
    
    # ## ADD : I want no cluster to be empty
    if (length(unique(curr.clus))<K){
      add.lossclus = fast.log.loss(x=X,P.list=P.list,mu=theta.list,temperature=temperature)
      add.wi       = mygibbs.step.b(t(add.lossclus$dist), kappa=my.kappa)
      curr.clus    = mygibbs.step.c(add.wi)  
    }
    record.clus[[i]] = curr.clus

    #Set new means based on closest clustering
    for(k in 1:K){
      idnow = which(curr.clus==k)
      if (length(idnow)==1){
        theta.temp = X[idnow,]
      } else if (length(idnow)>1){
        theta.temp = apply(X[curr.clus==k,],2,mean);  
      } else {
        stop("* there is no cluster.")
      }
      theta.list[[k]]=theta.temp
      theta.list[[k]] = Null(U.list[[k]])%*%t(Null(U.list[[k]]))%*%theta.list[[k]]  
    }
    
    #increase or decrease variance to adjust acceptance rates
    if(tune == 100){
      for(k in 1:K){
        if(tune.accept[k]<10){
          sub.met.sd[k] = .1*sub.met.sd[k]
        }else{
          if(tune.accept[k]<30){
            sub.met.sd[k]=.5*sub.met.sd[k]
          }else{
            if(tune.accept[k]<60){
              sub.met.sd[k]=2*sub.met.sd[k]
            }else{
              if(tune.accept[k]<90){
                sub.met.sd[k] = 5*sub.met.sd[k]
              }else{
                sub.met.sd[k]=10*sub.met.sd[k]
              }
            }
          }
        }
        tune.accept[k]=0
      }
      tune=0
    }
    
    
    #For storage at the end
    P.mat[[i]] = P.list
    U.mat[[i]] = U.list
    theta.mat[[i]]=theta.list
  }
  
  # return output
  output = list()
  for (i in 1:iter){
    iterate = list()
    iterate$P = P.mat[[i]]
    iterate$U = U.mat[[i]]
    iterate$theta = theta.mat[[i]]
    iterate$cluster = record.clus[[i]]
    output[[i]] = iterate
  }
  return(output)
}