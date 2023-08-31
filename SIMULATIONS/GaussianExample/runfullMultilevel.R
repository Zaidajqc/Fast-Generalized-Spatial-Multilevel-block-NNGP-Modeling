
############################################################################
## Author: Z. Quiroz, M. O. Prates, Zhiyong Hu and  D. Dey.
## Date: 31.08.2023
##
## Description:
##
##    The code performs a Bayesian analysis of full multilevel block-NNGP models using
##    Integrated Nested Laplace approximation (INLA). The BlockNNGP latent effect is implemented
##    using the rgeneric funcion. 


#############################################################################


setwd("/set_your_directory/")

## full multilevel coefficients models - gaussian distribution

source("find neighbor.r")

library(Matrix)
library(parallel)
library(sp)
library(pdist)
library(invgamma)
library(akima)
library(fields)
library(INLA)
library(inlabru)



inla.rgeneric.blocknngp.model <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {

  envir = parent.env(environment())
  #Internal function

  interpret.theta <- function() {
    return(
      list(var1 = exp(theta[1L]),
           var2 = exp(theta[2L]))
    )
  }

  graph = function() {
    require(Matrix)

    return (Q())
  }

  Q = function() {
    require(Matrix)
    require(parallel)

    param = interpret.theta()

    sigma2<- param$var1^2
    phis<- 2/param$var2

    #BLOCK LEVEL DISTRIBUTION
    build.BF<-function(i){
      if(i==1){
        B<-0
        F<-sigma2*exp(-phis*D_b[[i]])
      }else{
        C_bN<-sigma2*exp(-phis*D_bN[[i]])
        C_N<-sigma2*exp(-phis*D_NN[[i]])
        C_b<-sigma2*exp(-phis*D_b[[i]])
        B<-C_bN%*%solve(C_N)
        F<-C_b-B%*%t(C_bN)
      }
      return(list(B=B,F=F))
    }
    ####Construct Precision Matrix
    res.BF<-mclapply(1:nb^2,build.BF,mc.cores = 4)
    K<-nzip
    construct.Bs<-function(i){
      if(i>1){
        nb.i<-dim(res.BF[[i]]$F)[1]
        Bs.star<-matrix(rep(0,K*nb.i),nrow = nb.i,ncol = K)
        Bs.star[,train.coord$block.ind==i]<-diag(nb.i)
        Bs.star[,which(train.coord$block.ind %in% neighbor.block$NN_ind[i-1,])]<- -res.BF[[i]]$B
      }else{
        nb.i<-dim(res.BF[[i]]$F)[1]
        Bs.star<-matrix(rep(0,K*nb.i),nrow = nb.i,ncol = K)
        Bs.star[,train.coord$block.ind==i]<-diag(nb.i)
      }
      return(Matrix(t(Bs.star),sparse = T))
    }

    Bs<-mclapply(1:nb^2,construct.Bs,mc.cores = 4)
    Bs<-do.call(cbind,Bs)

    F.inv<-bdiag(mclapply(1:nb^2,function(i){ solve(res.BF[[i]]$F)}))

    Q<-Bs%*%F.inv%*%t(Bs)

    return (Q)
  }

  mu = function() {
    return(numeric(0))
  }

  log.norm.const = function() {
    return(numeric(0))

  }

  log.prior <- function() {
    param = interpret.theta()
    res <- log(lam1) + log(lam2) -2*(param$var2)- (lam1/(param$var2)) - lam2*((param$var1))
    return(res)
  }


  initial <- function() {
    return(c(initial.sigma, initial.range))
  }


  quit <- function() {
    return(invisible())
  }

  if (!length(theta)) theta = initial()
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}



multilevelsim1 <- function(indseed){

seed <- indseed
set.seed(seed)


# True parameters

beta1_true <- 2;beta2_true<-1; beta3_true <- -1;

phis_true <- 20; sig2 <- 0.2;  tau<-0.03; p.cov<-0.5;
phis_true2 <- 10; sig22 <- 0.1;
phis_true3 <- 15; sig23 <- 0.05;


nzip <- 2500

coord<-cbind(c(10001:(nzip+10000)),runif(nzip),runif(nzip))

nhouse<-sample(1:5,nzip,replace = T)
coord<-cbind(coord,nhouse)
colnames(coord)<-c('zip','horizontal','vertical','nhouse')
coord<-as.data.frame(coord)


distM <- as.matrix(dist(coord[,c('horizontal','vertical')])) #Euclidean

R <- sig2*exp(-distM*phis_true)
R2 <- sig22*exp(-distM*phis_true2)
R3 <- sig23*exp(-distM*phis_true3)


cov2<- unlist(sapply(nhouse, function(x) rnorm(x,100,20)))

cov3 <- unlist(sapply(coord$nhouse, function(x) rpois(x,1)+1))

set.seed(seed)
trueZ<-mvtnorm::rmvnorm(1, rep(0,nzip), sigma=R, method = "chol")
trueZ<-as.vector(trueZ)

set.seed(seed)
trueZ2<-mvtnorm::rmvnorm(1, rep(0,nzip), sigma=R2, method = "chol")
trueZ2<-as.vector(trueZ2)

set.seed(seed)
trueZ3<-mvtnorm::rmvnorm(1, rep(0,nzip), sigma=R3, method = "chol")
trueZ3<-as.vector(trueZ3)

names(trueZ) <- coord$zip
names(trueZ2) <- coord$zip
names(trueZ3) <- coord$zip

set.seed(seed)
nclaims <- sapply(beta1_true + beta2_true*cov2 + beta3_true*cov3 +
                    rep(trueZ,coord$nhouse) +
                    cov2*rep(trueZ2,coord$nhouse) +
                    cov3*rep(trueZ3,coord$nhouse),
                  function(x) rnorm(1,mean=x,sd=sqrt(tau)),simplify = T)

sim_data<- data.frame(zip = names(nclaims), nclaims = as.numeric(nclaims),  cov2=cov2, cov3=cov3, z1= rep(trueZ,coord$nhouse),  z2=rep(trueZ2,coord$nhouse), z3= rep(trueZ3,coord$nhouse) )

ind.train<-sort(sample(nzip,floor(nzip*0.8)))
ind.test<-c(1:nzip)[-ind.train]


train.coord<-coord[ind.train,]
test.coord<-coord[ind.test,]


train<-sim_data[which(sim_data$zip%in%train.coord$zip),]
test<-sim_data[which(sim_data$zip%in%test.coord$zip),]


original.train.coord<-train.coord
original.test.coord<-test.coord
original.train<-train
original.test<-test

M=20


neighbor.train<-NNMatrix(original.train.coord[,c('horizontal','vertical')],M,1,'cb')
train.coord<-original.train.coord[neighbor.train$ord,]
train<-train[order(match(train$zip,train.coord$zip)),]

nzip<-dim(train.coord)[1]
nzip.test<-dim(test.coord)[1]


X <- cbind(1,train$cov1,train$cov2)
y <- train$nclaims

X.test<-cbind(1,test$cov1,test$cov2)
y.test<-test$nclaims

###########SPLIT BLOCKS#####
require(ggplot2)
nb<-10
v_cut <- as.numeric(cut(train.coord$vertical, nb))
h_cut <- as.numeric(cut(train.coord$horizontal, nb))

block.ind<-interaction(v_cut,h_cut)
levels(block.ind)<-1:nb^2
train.coord <- cbind(train.coord,block.ind)
train.coord<-train.coord[order(train.coord$block.ind),]
dist_0<-aggregate(train.coord$horizontal^2+train.coord$vertical^2,list(train.coord$block.ind),mean)
train.coord$block.ind<-as.factor(rep(match(dist_0$x,sort(dist_0$x)),table(train.coord$block.ind)))

block.center<-cbind(1:nb^2,aggregate(train.coord$horizontal,list(train.coord$block.ind),mean)[,2],
                    aggregate(train.coord$vertical,list(train.coord$block.ind),mean)[,2])

M=4 ####Number of nearest blocks

neighbor.block<-NNMatrix(block.center[,2:3],M,1,'cb',ord = order(block.center[,2]+block.center[,3]))

train.coord<-train.coord[order(match(train.coord$block.ind,neighbor.block$ord)),]
train.coord$block.ind<-rep(1:nb^2,table(train.coord$block.ind)[neighbor.block$ord])

train<-train[order(match(train$zip,train.coord$zip)),]

D_bN<-list()
D_NN<-list()
D_b<-list()
for (i in 1:nb^2) {
  if(i==1){
    D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, c("vertical","horizontal")]))

  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[i-1,])
    D_bN[[i]]<-as.matrix(pdist(train.coord[train.coord$block.ind==i, c("vertical","horizontal")],
                               nn.block[, c("vertical","horizontal")]))
    D_NN[[i]]<-as.matrix(dist(nn.block[,c("vertical","horizontal")]))
    D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, c("vertical","horizontal")]))
  }
}




## Calculate hyperparameters
#range <- (8*nu)/2*phis
d=2
prior.range  = c( 0.5*0.15,0.05)
prior.sigma  = c(0.25 ,0.05)
lam11 <- -log(prior.range[2])*prior.range[1]^(d/2)
lam22 <- -log(prior.sigma[2])/prior.sigma[1]

initial.range <- log(prior.range[1]) + 1
initial.sigma <- log(prior.sigma[1]) - 1


block.nngp.model <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                         nb=nb, D_bN=D_bN, D_NN=D_NN, D_b=D_b, lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma, debug = T)

block.nngp.model1 <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                          nb=nb, D_bN=D_bN, D_NN=D_NN, D_b=D_b,lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma,debug = T)

block.nngp.model2 <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                          nb=nb, D_bN=D_bN, D_NN=D_NN, D_b=D_b,lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma,debug = T)

idx<-match(train$zip,train.coord$zip)


svc_components <- ~  -1 + intercept  + cov2 + cov3 +
  spatint(idx, model = block.nngp.model) +
  spatcov1(idx, weights = cov2, model = block.nngp.model1) +
  spatcov2(idx, weights = cov3, model = block.nngp.model2)

# formula, with "." meaning "add all the model components":
svc_formula1 <-  nclaims ~ .

multilevel1 <- bru(svc_components,
                   like(
                     formula = svc_formula1,
                     data =  data.frame(intercept=1,train)
                   ),
                   options = list(
                     control.compute = list(waic = TRUE, cpo = TRUE,config=TRUE),
                     #                     control.inla = list(strategy = "simplified.laplace", int.strategy = "eb"),
                     verbose = FALSE
                   )
)


m.blocknngp<- multilevel1

round(m.blocknngp$summary.fixed[, c(1,4,3, 5)] ,3)

marg.tau<-inla.tmarginal(function(x){1/x},m.blocknngp$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(marg.tau)

marg.sig<-inla.tmarginal(function(x){(exp(x))^2},m.blocknngp$marginals.hyperpar$`Theta1 for spatint`)
inla.zmarginal(marg.sig)
marg.phis<-inla.tmarginal(function(x){2/exp(x)},m.blocknngp$marginals.hyperpar$`Theta2 for spatint`)
inla.zmarginal(marg.phis)

marg.sig2<-inla.tmarginal(function(x){exp(x)^2},m.blocknngp$marginals.hyperpar$`Theta1 for spatcov1`)
inla.zmarginal(marg.sig2)
marg.phis2<-inla.tmarginal(function(x){2/exp(x)},m.blocknngp$marginals.hyperpar$`Theta2 for spatcov1`)
inla.zmarginal(marg.phis2)

marg.sig3<-inla.tmarginal(function(x){exp(x)^2},m.blocknngp$marginals.hyperpar$`Theta1 for spatcov2`)
inla.zmarginal(marg.sig3)
marg.phis3<-inla.tmarginal(function(x){2/exp(x)},m.blocknngp$marginals.hyperpar$`Theta2 for spatcov2`)
inla.zmarginal(marg.phis3)

# posterior means
est.beta <- m.blocknngp$summary.fixed[,1]
estsig1<-inla.emarginal(function(x){(exp(x))^2},m.blocknngp$marginals.hyperpar$`Theta1 for spatint`)
estphi1<-inla.emarginal(function(x){2/exp(x)},m.blocknngp$marginals.hyperpar$`Theta2 for spatint`)
estsig2<-inla.emarginal(function(x){exp(x)^2},m.blocknngp$marginals.hyperpar$`Theta1 for spatcov1`)
estphi2<-inla.emarginal(function(x){2/exp(x)},m.blocknngp$marginals.hyperpar$`Theta2 for spatcov1`)
estsig3<-inla.emarginal(function(x){exp(x)^2},m.blocknngp$marginals.hyperpar$`Theta1 for spatcov2`)
estphi3<-inla.emarginal(function(x){2/exp(x)},m.blocknngp$marginals.hyperpar$`Theta2 for spatcov2`)
esttau<-inla.emarginal(function(x){1/x},m.blocknngp$marginals.hyperpar$`Precision for the Gaussian observations`)


covbeta1=0; if(beta1_true>= inla.hpdmarginal(0.95,m.blocknngp$marginals.fixed$intercept)[1] & beta1_true<= inla.hpdmarginal(0.95,m.blocknngp$marginals.fixed$intercept)[2]) covbeta1 =1
covbeta2=0; if(beta2_true>= inla.hpdmarginal(0.95,m.blocknngp$marginals.fixed$cov2)[1] & beta2_true<= inla.hpdmarginal(0.95,m.blocknngp$marginals.fixed$cov2)[2]) covbeta2 =1
covbeta3=0; if(beta3_true>= inla.hpdmarginal(0.95,m.blocknngp$marginals.fixed$cov3)[1] & beta3_true<= inla.hpdmarginal(0.95,m.blocknngp$marginals.fixed$cov3)[2]) covbeta3 =1
covsig1=0; if(sig2>= inla.hpdmarginal(0.95,marg.sig)[1] & sig2<=inla.hpdmarginal(0.95,marg.sig)[2]) covsig1 =1
covsig2=0; if(sig22>= inla.hpdmarginal(0.95,marg.sig2)[1] & sig22<=inla.hpdmarginal(0.95,marg.sig2)[2]) covsig2 =1
covsig3=0; if(sig23>= inla.hpdmarginal(0.95,marg.sig3)[1] & sig23<=inla.hpdmarginal(0.95,marg.sig3)[2]) covsig3 =1
covphi1=0; if(phis_true>= inla.hpdmarginal(0.95,marg.phis)[1] & phis_true<=inla.hpdmarginal(0.95,marg.phis)[2]) covphi1 =1
covphi2=0; if(phis_true2>= inla.hpdmarginal(0.95,marg.phis2)[1] & phis_true2<=inla.hpdmarginal(0.95,marg.phis2)[2]) covphi2 =1
covphi3=0; if(phis_true3>= inla.hpdmarginal(0.95,marg.phis3)[1] & phis_true3<=inla.hpdmarginal(0.95,marg.phis3)[2]) covphi3 =1
covtau=0; if(tau>= inla.hpdmarginal(0.95,marg.tau)[1] & tau<=inla.hpdmarginal(0.95,marg.tau)[2]) covtau =1

dic <- m.blocknngp$dic$dic
waic <- m.blocknngp$waic$waic
time <- sum(m.blocknngp$bru_timings[,3])
LPML <- sum(log(m.blocknngp$cpo$cpo))

names.mean.estparam <- c("est.beta1", "est.beta2", "est.beta3","estsig1", "estphi1", "estsig2", "estphi2", "estsig3", "estphi3", "esttau")
mean.estparam <- round(c(est.beta, estsig1, estphi1, estsig2, estphi2, estsig3, estphi3, esttau),3)
names.cov.estparam <- c("cov.beta1", "cov.beta2", "cov.beta3","covsig1", "covphi1", "covsig2", "covphi2", "covsig3", "covphi3", "covtau")
cov.param <- c(covbeta1, covbeta2, covbeta3, covsig1, covsig2, covsig3, covphi1, covphi2, covphi3, covtau)

#RMSE
orig.y <- train$nclaims
fit <- m.blocknngp$summary.fitted.values[1:length(orig.y ),1]
RMSEy = sqrt(sum((orig.y- fit)^2)/length(orig.y ))


rz1<-m.blocknngp$summary.random$spatint$mean
rz2<-m.blocknngp$summary.random$spatcov1$mean
rz3<-m.blocknngp$summary.random$spatcov2$mean
unik <- !duplicated(idx)
RMSEz1 = sqrt(sum((rz1-train$z1[unik])^2)/length(rz1))
RMSEz2 = sqrt(sum((rz2-train$z2[unik])^2)/length(rz2))
RMSEz3 = sqrt(sum((rz3-train$z3[unik])^2)/length(rz3))



X <- as.matrix(cbind(1,train[,c( "cov2", "cov3")]))

LP1 <-  (X[,1]*est.beta[1] + rep(rz1,train.coord$nhouse))
LP2 <-  (X[,2]*est.beta[2] + rep(rz2,train.coord$nhouse))
LP3 <-  (X[,3]*est.beta[3] + rep(rz3,train.coord$nhouse))
trueLP1 <- (X[,1]*beta1_true + rep(train$z1[unik],train.coord$nhouse))
trueLP2 <- (X[,2]*beta2_true + rep(train$z2[unik],train.coord$nhouse))
trueLP3 <- (X[,3]*beta3_true + rep(train$z3[unik],train.coord$nhouse))
RMSELP1 = sqrt(sum((trueLP1-LP1)^2)/length(LP1))
RMSELP2 = sqrt(sum((trueLP2-LP2)^2)/length(LP2))
RMSELP3 = sqrt(sum((trueLP3-LP3)^2)/length(LP3))


# save true values
name.true.param <-  c("beta1_true", "beta2_true", "beta3_true" , "sig1", "phis_true1", "sig2", "phis_true2", "sig3", "phis_true3", "tau")
true.param <-  c(beta1_true, beta2_true, beta3_true , sig2, phis_true, sig22, phis_true2, sig23, phis_true3, tau)
final.true <- data.frame(res=name.true.param, val=true.param )
save(res=final.true, file="multilevel_sim1_trueparam.Rdata")

##### end estimation ###

# Prediction assessment
n.sample <- 100 #1e3
X1 <- data.frame(X[,-1])
set.seed(seed)
# pred version 1 uses generate... but it makes biased predictions...
#multilevel1.pos <- generate(m.blocknngp, X1, n.samples = n.sample)
# set number of neighbor blocks for predictions
# it can be Mp=1 or Mp = M, for other Mp you need to implement it...

# for intercept process
predict.spatial1<- function(j="iterations",i="ith location"){
  Mp=M  
  sigma2<- estsig1
  phis<- estphi1
  
  b.i<-train.coord$block.ind[which.min(spDistsN1(as.matrix(train.coord[,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))]
  
  # Mp =1 or Mp=M
  if(Mp==1|b.i==1){
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))
    C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% rz1[which(train.coord$block.ind==b.i)],sqrt(F.i))
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[b.i-1,])
    b.ii<- c(b.i, unique(nn.block$block.ind))
    indblockpred <- NULL
    for(k in 1:length(b.ii)){
      indp <- which(train.coord$block.ind==b.ii[k])
      indblockpred <- c(indblockpred,indp)
    }
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[indblockpred,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))
    D_b1<-as.matrix(dist(train.coord[indblockpred, c("vertical","horizontal")]))
    C_N.inv<-solve(sigma2*exp(-phis*D_b1))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% rz1[indblockpred],sqrt(F.i))
  }
  return(z.i)
}



pre_Z.1<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.1)<-test.coord$location
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res1<-mclapply(1:n.sample,predict.spatial1,i=i,mc.cores = 4)
  pre_Z.1[,i]<-unlist(res1)
}

# for beta1 process
predict.spatial2<- function(j="iterations",i="ith location"){
  Mp=M
  sigma2<- estsig2
  phis<- estphi2
  
  b.i<-train.coord$block.ind[which.min(spDistsN1(as.matrix(train.coord[,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))]
  
  # Mp =1 or Mp=M
  if(Mp==1|b.i==1){
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))
    C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% rz2[which(train.coord$block.ind==b.i)],sqrt(F.i))
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[b.i-1,])
    b.ii<- c(b.i, unique(nn.block$block.ind))
    indblockpred <- NULL
    for(k in 1:length(b.ii)){
      indp <- which(train.coord$block.ind==b.ii[k])
      indblockpred <- c(indblockpred,indp)
    }
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[indblockpred,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))
    D_b1<-as.matrix(dist(train.coord[indblockpred, c("vertical","horizontal")]))
    C_N.inv<-solve(sigma2*exp(-phis*D_b1))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% rz2[indblockpred],sqrt(F.i))
  }
  return(z.i)
}


pre_Z.2<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.2)<-test.coord$location
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res2<-mclapply(1:n.sample,predict.spatial2,i=i,mc.cores = 4)
  pre_Z.2[,i]<-unlist(res2)
}

#for b2 process
predict.spatial3<- function(j="iterations",i="ith location"){
  Mp=M
  sigma2<- estsig3
  phis<- estphi3
  
  b.i<-train.coord$block.ind[which.min(spDistsN1(as.matrix(train.coord[,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))]
  
  # Mp =1 or Mp=M
  if(Mp==1|b.i==1){
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))
    C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% rz3[which(train.coord$block.ind==b.i)],sqrt(F.i))
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[b.i-1,])
    b.ii<- c(b.i, unique(nn.block$block.ind))
    indblockpred <- NULL
    for(k in 1:length(b.ii)){
      indp <- which(train.coord$block.ind==b.ii[k])
      indblockpred <- c(indblockpred,indp)
    }
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[indblockpred,c("horizontal","vertical")]), as.matrix(test.coord[i,c("horizontal","vertical")])))
    D_b1<-as.matrix(dist(train.coord[indblockpred, c("vertical","horizontal")]))
    C_N.inv<-solve(sigma2*exp(-phis*D_b1))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% rz3[indblockpred],sqrt(F.i))
  }
  return(z.i)
}

pre_Z.3<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.3)<-test.coord$location
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res3<-mclapply(1:n.sample,predict.spatial3,i=i,mc.cores = 4)
  pre_Z.3[,i]<-unlist(res3)
}

plot(test$z1, rep(colMeans(pre_Z.1),test.coord$nhouse),xlim=c(-2,2), ylim=c(-2,2))
abline(0,1, col="red")
plot(test$z2, rep(colMeans(pre_Z.2),test.coord$nhouse),xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
abline(0,1, col="red")
plot(test$z3, rep(colMeans(pre_Z.3),test.coord$nhouse),xlim=c(-1,1), ylim=c(-1,1))
abline(0,1, col="red")

sqrt(mean((test$z1- rep(colMeans(pre_Z.1),test.coord$nhouse))^2))
sqrt(mean((test$z2- rep(colMeans(pre_Z.2),test.coord$nhouse))^2))
sqrt(mean((test$z3- rep(colMeans(pre_Z.3),test.coord$nhouse))^2))


X.test <- as.matrix(cbind(1,test[,c("cov2","cov3")]))

get.pre<-function(i){
  precision<- esttau # nugget effect
  coef <- est.beta
  mu <- X.test%*%coef+
              rep(pre_Z.1[i,],test.coord$nhouse)+
              test$cov2*rep(pre_Z.2[i,],test.coord$nhouse)+
              test$cov3*rep(pre_Z.3[i,],test.coord$nhouse)
  res2<-rnorm(n=dim(X.test)[1], mu , sqrt(precision) )
  return(res2)
}


y.pre<-mclapply(1:n.sample, get.pre, mc.cores = 4)
y.pre<-do.call(rbind, y.pre)
pre.block.nngp2<-colMeans(y.pre)
y.test <- test$nclaims

plot(y.test,pre.block.nngp2 )
abline(0,1, col="grey30")

plot0 <- data.frame(x= y.test, y= pre.block.nngp2)

plopt_pred = ggplot(data = plot0,mapping = aes(x=x,y=y))+geom_point(alpha = 1/10)+ 
  geom_abline(colour = "grey50")+ 
  labs(x = "y", y = "y_pred")#+  coord_flip( xlim = c(0, 250),ylim = c(0,250))

jpeg(file="pred100_4normalNEW2.jpeg")
plopt_pred
dev.off()

RMSP<- sqrt(mean((pre.block.nngp2-y.test)^2))


# return final res
name.crit <- c("dic", "waic", "LPML", "RMSEy", "RMSEz1", "RMSEz2", "RMSEz3", "RMSELP1", "RMSELP2", "RMSELP3","time(sec)", "RMSPy")
crit.model<- round(c( dic, waic, LPML, RMSEy, RMSEz1, RMSEz2, RMSEz3, RMSELP1, RMSELP2, RMSELP3, time, RMSP),3)
name.full.results <- c(t(names.mean.estparam), t(names.cov.estparam), t(name.crit))
full.results <- c(t(mean.estparam), t(cov.param), t(crit.model) )

save(res=m.blocknngp, file="blocknngp1.Rdata")

return(full.results)

}

## simulation niter repetitions
names.mean.estparam <- c("est.beta1", "est.beta2", "est.beta3","estsig1", "estphi1", "estsig2", "estphi2", "estsig3", "estphi3", "esttau")
names.cov.estparam <- c("cov.beta1", "cov.beta2", "cov.beta3","covsig1", "covphi1", "covsig2", "covphi2", "covsig3", "covphi3", "covtau")
name.crit <- c("dic", "waic", "LPML", "RMSEy", "RMSEz1", "RMSEz2", "RMSEz3", "RMSELP1", "RMSELP2", "RMSELP3","time(sec)", "RMSPy")
name.full.results <- c(t(names.mean.estparam), t(names.cov.estparam), t(name.crit))

final.res <- data.frame(res=name.full.results)
niter <- 30
for(indseed in (1:niter)){
  full.results <- multilevelsim1(indseed)
  final.res <- cbind(final.res, full.results)
  save(res=final.res, file=paste("multilevel_sim1_",indseed,".Rdata"))
  final.res=NULL
  full.results=NULL
}

## plots for paper


setwd("/set_your_directory/")
for(indseed  in(1:30)){
  load(paste("multilevel_sim1_",indseed,".Rdata"))
  if(indseed==1) {res=final.res}
  if(indseed>1) {
    res = (cbind(res, final.res))
  }
}

save(res=res, file="ALLmultilevel_sim1_100_4_Normal.Rdata")


## plots for paper
library(ggplot2)

setwd("/set_your_directory/")
load(paste("multilevel_sim1_trueparam.Rdata"))

load("ALLmultilevel_sim1_100_4_Normal.Rdata")
niter = dim(res)[2]-1

niter <- 30
final.cov <- data.frame(cov=rowSums(res[11:20,2:31])/niter)
final.cov <- cbind(res[11:20,1], round(final.cov,3))
names(final.cov ) <-  c("param","cov" )
save(res=final.cov, file=paste("coverage.Rdata"))


niter = dim(res)[2]-1

#boxplot beta1, beta2, beta3
df <- data.frame(values=(rbind(t(res[1,2:(niter+1)]), t(res[2,2:(niter+1)]), t(res[3,2:(niter+1)]))),
                 category2=(c(rep("gamma11",niter), rep("gamma12",niter), rep("gamma13",niter))),
                 param=c(rep(final.true[1,2], niter), rep(final.true[2,2],niter), rep(final.true[3,2],niter)))
names(df) <- c("Posterior_means", "cat", "param")


df$cat <- factor(df$cat,
                 labels = c(bquote(gamma[11]),
                            bquote(gamma[12]),
                            bquote(gamma[13])))

p0 <- ggplot(df, aes(x=1, y = Posterior_means)) +
  geom_boxplot() +
  facet_wrap(~cat, nrow=1,
             labeller = label_parsed, scales = "free_y") +
  scale_x_discrete( "") + #,   labels ="blockNNGP"
  geom_hline(data = df,
             aes(yintercept = param), linewidth = 1, color="red",linetype = "dashed")+
  theme_bw(base_size=12)

p0

#boxplot phi1, phi2, phi3
df1 <- data.frame(values=(rbind(t(res[5,2:(niter+1)]), t(res[7,2:(niter+1)]), t(res[9,2:(niter+1)]))),
                  category2=(c(rep(1,niter), rep(2,niter), rep(3,niter))),
                  param=c(rep(final.true[5,2], niter), rep(final.true[7,2],niter), rep(final.true[9,2],niter)))
names(df1) <- c("Posterior_means", "cat", "param")

df1$cat <- factor(df1$cat,
                  labels = c(bquote(phi[1]),
                             bquote(phi[2]),
                             bquote(phi[3])))

p1 <- ggplot(df1, aes(x=1, y = Posterior_means)) +
  geom_boxplot() +
  facet_wrap(~cat, nrow=1,
             labeller = label_parsed) +
  scale_x_discrete( "") + #,   labels ="blockNNGP"
  geom_hline(data = df1,
             aes(yintercept = param), linewidth = 1, color="red",linetype = "dashed")+
  theme_bw(base_size=12)

#boxplot sigma1, sigma2, sigma3

df2 <- data.frame(values=(rbind(t(res[4,2:(niter+1)]), t(res[6,2:(niter+1)]), t(res[8,2:(niter+1)]))),
                  category2=(c(rep(1,niter), rep(2,niter), rep(3,niter))),
                  param=c(rep(final.true[4,2], niter), rep(final.true[6,2],niter), rep(final.true[8,2],niter)))
names(df2) <- c("Posterior_means", "cat", "param")


df2$cat <- factor(df2$cat,
                  labels = c(bquote(sigma[1]^2),
                             bquote(sigma[2]^2),
                             bquote(sigma[3]^2)))

p2 <- ggplot(df2, aes(x=1, y = Posterior_means)) +
  geom_boxplot() +
  facet_wrap(~cat, nrow=1,
             labeller = label_parsed) +
  scale_x_discrete( "") + #,   labels ="blockNNGP"
  geom_hline(data = df2,
             aes(yintercept = param), linewidth = 1, color="red",linetype = "dashed")+
  theme_bw(base_size=12)



# boxplot nugget
df4 <- data.frame(values=t(res[10,2:(niter+1)]), category2 = c(rep(1,niter)))
names(df4) <- c("Posterior_means", "cat")

df4$cat <- factor(df4$cat,
                  labels = c(bquote(tau^2)))

p4 <- ggplot(df4, aes(x=1, y = Posterior_means)) +
  geom_boxplot() +
  facet_wrap(~cat, nrow=1,
             labeller = label_parsed, scales = "free_y") +
  scale_x_discrete( "")   + #,   labels ="blockNNGP"
  geom_hline(data = df4,
             aes(yintercept = final.true[10,2]), size = 1, color="red",linetype = "dashed")+
  theme_bw(base_size=12)

library(ggpubr)

#boxplot of parameter posterior mean estimates
jpeg(file="sim1_plot0.jpeg")
ggarrange(p0, p11, p2,
          ggarrange(p4, ncol = 2),
          nrow = 4)
dev.off()

#boxplot of parameter posterior mean estimates
jpeg(file="sim1_plot1.jpeg")
ggarrange(p0, p1, p2,
          ggarrange(p4, ncol = 2),
          nrow = 4)
dev.off()

# boxplot RMSE spatial random effects
df3 <- data.frame(values=(rbind(t(res[25,2:(niter+1)]), t(res[26,2:(niter+1)]), t(res[27,2:(niter+1)]))),
                  category2=(c(rep(1,niter), rep(2,niter), rep(3,niter))))
names(df3) <- c("RMSE", "cat")

df3$cat <- factor(df3$cat,
                  labels = c(bquote(alpha[1]),
                             bquote(alpha[2]),
                             bquote(alpha[3])))

p3 <- ggplot(df3, aes(x=1, y = RMSE)) +
  geom_boxplot() +
  facet_wrap(~cat, nrow=1,
             labeller = label_parsed) +
  scale_x_discrete( "") + #,   labels ="blockNNGP"
  theme_bw(base_size=20)

jpeg(file="sim1_plot2.jpeg")
p3
dev.off()

# plot WAIC, LPML, time, RMSEy

# criteria assessment
df4 <- data.frame(values=(rbind(t(res[22,2:(niter+1)]),
                                t(res[23,2:(niter+1)]))),
                  category2=(c(rep(1,niter), rep(2,niter))))
names(df4) <- c("values", "cat")

df4$cat <- factor(df4$cat,
                  labels = c("WAIC",
                             "LPML"))

p4 <- ggplot(df4, aes(x=1, y = values)) +
  geom_boxplot() +
  facet_wrap(~cat, nrow=1,
             labeller = label_parsed, scales = "free_y") +
  scale_x_discrete( "") + #,   labels ="blockNNGP"
  theme_bw(base_size=20)


df5 <- data.frame(values=(rbind( t(res[31,2:(niter+1)]),
                                 t(res[24,2:(niter+1)]))),
                  category2=(c(rep(1,niter), rep(2,niter))))
names(df5) <- c("values", "cat")

df5$cat <- factor(df5$cat,
                  labels = c("Time(sec)",
                             "RMSE_y"))

p5 <- ggplot(df5, aes(x=1, y = values)) +
  geom_boxplot() +
  facet_wrap(~cat, nrow=1,
             labeller = label_parsed, scales = "free_y") +
  scale_x_discrete( "") + #,   labels ="blockNNGP"
  theme_bw(base_size=20)

jpeg(file="sim1_plot3.jpeg")
ggarrange(p4, p5, ncol = 1, nrow = 2)
dev.off()


