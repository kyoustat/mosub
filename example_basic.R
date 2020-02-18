# example_basic -----------------------------------------------------------
#
# In this script, we do the followings,
#   (a) show generated data for line-and-plane example,
#   (b) prediction using the fitted model, and
#   (c) impact of isotropic variance in clustering efficacy.
#
# Warning that (c) is quite time-consuming, taking about an hour in a normal desktop
# if run without parallel computing. If you want to run (c) as well, 
# change the following logical as 'run.part.c = TRUE'.

rm(list=ls())
run.part.c = FALSE

# preliminary : install and load required libraries -----------------------
# install
package.cran = c("devtools","scatterplot3d","mclustcomp")
for (pkg in package.cran){
  if (!is.element(pkg, installed.packages()[,1])){
    install.packages(pkg, dependencies = TRUE)
  }  
}
if (!is.element("mosub", installed.packages()[,1])){
  devtools::install_github("kyoustat/mosub")
}
# load libraries
library(mosub)
library(mclustcomp)
library(scatterplot3d)


# (a) show generated data for line-and-plane example ----------------------
set.seed(6)
mydata = mosub::gen.LP(n=500, K=2, iso.var=0.25) 
X      = mydata$data
label  = mydata$class
if (run.part.c){
  par(mfrow=c(2,2))  
} else {
  par(mfrow=c(1,3))
}
scatterplot3d(x=X, pch=19, cex.symbols=0.5, color=label+2,
              main="(a) 3d visualization",xlab = "",ylab="",zlab="")



# (b) prediction using the fitted model -----------------------------------
# (b-1) separate as train/test data
id.train = sample(1:1000, 700, replace=FALSE) # select 700 data points for training
id.test  = setdiff(1:1000, id.train)          # rest of 300 data are for prediction

train.dat = mydata$data[id.train,]
train.lab = mydata$class[id.train]

test.dat = mydata$data[id.test,]
test.lab = mydata$class[id.test]

# (b-2) run MSM algorithm with K=2
maxiter   = 10000
output2   = msm(train.dat, K=2, iter=maxiter)
test.pred = msmpred(test.dat, output2)

# (b-3) visualize the test data
scatterplot3d(x=test.dat, color=test.lab+1,  pch=19, cex.symbols=0.5,
              xlab = "",ylab="",zlab="",main="(b) true label")
scatterplot3d(x=test.dat, color=test.pred+1, pch=19, cex.symbols=0.5,
              xlab = "",ylab="",zlab="",main="(b) predicted label")

# (c) impact of isotropic variance in clustering efficacy -----------------
if (run.part.c){
  # (c-1) set up parameters
  vec.iso  = c(0.1,0.5,1,2,5)
  vec.mean = rep(0,5)
  vec.sdev = rep(0,5)
  
  # (c-2) let's do the iteration
  iter = 0
  for (i in 1:5){        # loop for levels of isotropic variance
    now.iso = vec.iso[i] # current variance value
    tmp.rec = rep(0,10)  # temporary recording for 10 iterations at each setting
    for (j in 1:10){
      # generate data
      mydata = mosub::gen.LP(n=500, K=2, iso.var=now.iso)
      # separate data as in (b)
      id.train = sample(1:1000, 700, replace=FALSE) # select 700 data points for training
      id.test  = setdiff(1:1000, id.train)          # rest of 300 data are for prediction
      
      train.dat = mydata$data[id.train,]
      test.dat  = mydata$data[id.test,]
      test.lab  = mydata$class[id.test]
      # run the algorithm and do the prediction
      output2   = msm(train.dat, K=2, iter=maxiter, print.progress = FALSE)
      test.pred = msmpred(test.dat, output2, print.progress = FALSE)
      # compute Jaccard index between two labels
      tmp.rec[j] = as.double(mclustcomp(test.lab, test.pred, types="jaccard"))[2]
      # show progress
      iter = iter + 1
      print(paste("(c) iteration ",iter,"/50 complete...",sep=""))
    }
    vec.mean[i] = base::mean(tmp.rec) # record mean
    vec.sdev[i] = stats::sd(tmp.rec)  #    and standard deviation of Jaccard indices
  }
  # (c-3) visualize
  plot(1:5, vec.mean,
       ylim=range(c(vec.mean-vec.sdev, vec.mean+vec.sdev)),
       pch=19, xlab="isotropic variances", ylab="Jaccard index",
       main="(c) effect of isotropic variances")
  arrows(1:5, vec.mean-vec.sdev, 1:5, vec.mean+vec.sdev, length=0.05, angle=90, code=3)
}
