
gc()
rm(list=ls())


library("INLA")
library(inlabru)
library(sp)
library(ggplot2)
library(parallel)

setwd("/set_your directory/")

load("/set_your directory/NEWbeijing.RData")

dir.save=getwd()

####################
# Fig.9 
####################
# for map of China
library(rgdal)
library(data.table)
library(mapproj)

source("~/mapChina.R")

train2 <- aggregate(train$totalPrice, list(train$location), FUN=length)
train3 <- aggregate(train$totalPrice, list(train$location), FUN=mean) #sum or mean?
train4 <- aggregate(train$Lng, list(train$location), FUN=mean)
train5 <- aggregate(train$Lat, list(train$location), FUN=mean)

finaltrain <- data.frame( totalPrice=train3[,2], 
                         Lng=train4[,2], Lat=train5[,2], count=train2[,2])


jpeg(paste(dir.save,"/Price_building.jpeg",sep=""))
ggplot(data = finaltrain)+
  geom_path(data = dt_Beijing, aes(long, lat, group = group),
            size = 0.3, color = "grey60") +
  geom_point(aes(Lng, Lat, color = log(totalPrice), size=count))+ 
  labs(x = "Longitude", y = "Latitude")+
  coord_fixed(xlim = c(116, 117), 
              ylim = c(39.6, 40.3))
dev.off()

# histogram of total price data
jpeg(paste(dir.save,"/Price_buildingHistogram.jpeg",sep=""))
ggplot(train, aes(x=totalPrice))  + 
  geom_histogram(color="grey", fill="white")
dev.off()


####################
## Run gamma models
####################

rm(list=ls())

load("/set_your directory/NEWbeijing.RData")

dir.save=getwd()

idx<-match(train$location,train.coord$location)
idx1<-match(test$location,test.coord$location)


#############################################
############## Non-spatial  ##############
#############################################

# multilevel gamma model in intercept (non-spatial)

svc_formula0 <- totalPrice~ -1 + intercept +  area  + livingRoom + 
  gamma_int_elev +   gamma_int_subway +    gamma_int_age 


multilevelNS_gamma <- bru( svc_formula0,
                           family = "gamma",
                           data = train,
                           options = list(
                             control.compute = list(waic = TRUE, cpo = TRUE),
                             verbose = FALSE
                           )
)  

## save non-spatial model
save(multilevelNS_gamma, file = "multilevelNS_gamma.RData")

round(multilevelNS_gamma$summary.hyperpar[, c(1, 5,3, 4)],3)

round(multilevelNS_gamma$summary.fixed[, c(1, 5,3, 4)] ,3)

dic= multilevelNS_gamma$dic$dic
waic = multilevelNS_gamma$waic$waic
LPML = sum(log(multilevelNS_gamma$cpo$cpo))
time = sum(multilevelNS_gamma$bru_timings[,3])


est0 <- multilevelNS_gamma$summary.fitted.values

nzip<-dim(train.coord)[1]


y<- (train$totalPrice)
est0 <- multilevelNS_gamma$summary.fitted.values$mean[1:length(y)]
RMSE<- sqrt(mean((est0-y)^2))

plot0 <- data.frame(location=train$location, totalPrice=train$totalPrice, 
                    Lng=train$Lng, Lat=train$Lat, pred=fit0[,1] )
dir.save=getwd()
jpeg(paste(dir.save,"/est_gamma_multilevelNS",".jpeg",sep=""))
ggplot(data = plot0,mapping = aes(x=totalPrice,y=pred))+geom_point()+
  geom_abline(colour = "grey50")+ 
  labs(x = "y", y = "y_est")
dev.off()



X <- as.matrix(cbind(1,train[,c("area","livingRoom","gamma_int_elev",
                                "gamma_int_subway", 
                                "gamma_int_age")]))
X0 <- data.frame(X[,-1])


n.sample <- 100 #1e4
set.seed(5)
multilevel0.pos <- generate(multilevelNS_gamma, X0, n.samples = n.sample)
multilevel0.pos.hyper <- inla.hyperpar.sample(n = n.sample, result=multilevelNS_gamma)

filename2= paste("multilevelNS_gamma_posterior.Rdata",sep="")
save(multilevel0.pos, file=filename2)
filename2= paste("multilevelNS_gamma_posteriorprec.Rdata",sep="")
save(multilevel0.pos.hyper, file=filename2)



get.est0<-function(i){
  precision<-multilevel0.pos.hyper[[i]]
  coef <- matrix(unlist(multilevel0.pos[[i]][2:7]), 6,1)
  mu <- exp(X%*%coef)
  b         <- mu / precision
  res10<-rgamma(n=dim(X)[1], shape = precision ,  scale = b)
  return(res10)
}
y.rep0<-mclapply(1:n.sample,get.est0,mc.cores = 4)
y.rep0<-do.call(cbind,y.rep0)



X.test <- as.matrix(cbind(1,test[,c("area","livingRoom","gamma_int_elev",
                                    "gamma_int_subway", 
                                    "gamma_int_age")]))


get.pre<-function(i){
  precision<-multilevel0.pos.hyper[[i]]#$Precision_parameter_for_the_Gamma_observations
  coef <- matrix(unlist(multilevel0.pos[[i]][2:7]), 6,1)
  mu <- exp(X.test%*%coef)
  b         <- mu / precision
  res2<-rgamma(n=dim(X.test)[1], shape = precision ,  scale = b)
  return(res2)
}


y.pre1<-mclapply(1:n.sample, get.pre, mc.cores = 4)
y.pre<-do.call(cbind, y.pre1)

pre.block.nngp2<-rowMeans(y.pre)

est.block.nngp<-rowMeans(y.rep0)

RMSEy<- sqrt(mean((est.block.nngp-y)^2))

RMSP<- sqrt(mean((pre.block.nngp2-y.test)^2))

# save final res
name.crit <- c("dic", "waic", "LPML", "RMSEy","time(sec)", "RMSPy")
crit.model<- round(c( dic, waic, LPML, RMSE, time, RMSP),3)
res=data.frame(name.crit, val=crit.model)
save(res, file="multilevelNS_gamma_summary.Rdata")



######################################################
##############Block-NNGP MODELS ######################
#####################################################

rm(list=ls())

load("/set_your directory/NEWbeijing.RData")
dir.save=getwd()

# Set the number of blocks L and the number of neighbors M
L=8 #2^L blocks
M=8 #number of neighbors

idx<-match(train$location,train.coord$location)
idx1<-match(test$location,test.coord$location)

##Blocking
#L=8 #2^L blocks

train.coord$block.ind<-1
train.coord$ind<-1
for (i in 1:L) {
  if(i%%2==1){
    d<- 1
  }else{d<- 2}
  K<- 2^(i-1)
  for (k in 1:K) {
    cut.point<-median(train.coord[train.coord$ind==k ,d])
    train.coord[train.coord$ind==k &train.coord[,d]<=cut.point,"block.ind"]<- 2*k-1
    train.coord[train.coord$ind==k &train.coord[,d]>cut.point,"block.ind"]<-2*k
  }
  train.coord$ind<- train.coord$block.ind
}
train.coord$ind<-NULL
ggplot(data = train.coord,mapping = aes(x=Lng,y=Lat,color=as.factor(block.ind)))+geom_point()+
  geom_text(aes(label=as.factor(block.ind)),hjust=0, vjust=0)+ theme(legend.position = "none")


block.center<-cbind(aggregate(train.coord$Lng,list(train.coord$block.ind),mean)[,2],
                    aggregate(train.coord$Lat,list(train.coord$block.ind),mean)[,2],
                    1:2^L)
block.center<-as.data.frame(block.center)

NN_ind<-matrix(rep(0, (dim(block.center)[1]-1)*M),nrow = dim(block.center)[1]-1,ncol = M)
NN_distM<-matrix(rep(0, (dim(block.center)[1]-1)*M*(M-1)/2),nrow = dim(block.center)[1]-1,ncol = M*(M-1)/2)
NN_dist<-matrix(rep(0, (dim(block.center)[1]-1)*M),nrow = dim(block.center)[1]-1,ncol = M)

for (i in 2:dim(block.center)[1]) {
  #compute great circle distance
  distance<-spDistsN1(as.matrix(block.center[1:(i-1),1:2]),as.matrix(block.center[i,1:2]),longlat = T) 
  near.distance<-sort(distance,decreasing = FALSE,index.return=TRUE)
  dist.NN<-near.distance$x[1:min(i-1,M)] #dist to nearest neighbors
  ind.NN<-near.distance$ix[1:min(i-1,M)] #index of nearest neighbors
  distM.NN<-spDists(as.matrix(block.center[ind.NN,1:2]),longlat = T)
  
  NN_ind[i-1,1:min(i-1,M)]<-ind.NN
  if(i>=3){
    NN_distM[i-1,1:(min(i-1,M)*(min(i-1,M)-1)/2)]<-distM.NN[lower.tri(distM.NN,diag = FALSE)]
  }
  NN_dist[i-1,1:min(i-1,M)]<-dist.NN
}

neighbor.block<-list(NN_ind=NN_ind, NN_distM=NN_distM, NN_dist=NN_dist)



train.coord<-train.coord[order(train.coord$block.ind),]
train<-train[order(match(train$location,train.coord$location)),]

D_bN<-list()
D_NN<-list()
D_b<-list()
for (i in 1:2^L) {
  if(i==1){
    D_b[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),longlat = T)
    
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[i-1,])
    D_bN[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),as.matrix(nn.block[, c("Lng","Lat")]),longlat = T)
    D_NN[[i]]<-spDists(as.matrix(nn.block[, c("Lng","Lat")]),longlat = T)
    D_b[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),longlat = T)
  }
}


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
    res.BF<-mclapply(1:2^nb,build.BF,mc.cores = 4)
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
    
    Bs<-mclapply(1:2^nb,construct.Bs,mc.cores = 4)
    Bs<-do.call(cbind,Bs)
    
    F.inv<-bdiag(mclapply(1:2^nb,function(i){ solve(res.BF[[i]]$F)}))
    
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


library("INLA")


## Calculate hyperparameters
#range <- (8*nu)/2*phis
d=2
prior.range  = c(0.1,0.05)
prior.sigma  = c(1 ,0.05)
lam11 <- -log(prior.range[2])*prior.range[1]^(d/2)
lam22 <- -log(prior.sigma[2])/prior.sigma[1]

initial.range <- log(prior.range[1]) + 1
initial.sigma <- log(prior.sigma[1]) - 1


block.nngp.model <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                         nb=L, D_bN=D_bN, D_NN=D_NN, D_b=D_b, lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma, debug = T)

block.nngp.model1 <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                          nb=L, D_bN=D_bN, D_NN=D_NN, D_b=D_b,lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma,debug = T)

block.nngp.model2 <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                          nb=L, D_bN=D_bN, D_NN=D_NN, D_b=D_b,lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma,debug = T)



################################################################
############## Multilevel spatial in intercept ##############
###############################################################

svc_formula1 <- totalPrice ~ -1 + intercept +  area  + livingRoom + 
  gamma_int_elev +   gamma_int_subway + 
  gamma_int_age + 
  spatintercept(idx, model = block.nngp.model)


multilevelS_gamma1 <- bru( svc_formula1,
              family = "gamma",
              data = train,
            options = list(
              control.compute = list(waic = TRUE, cpo = TRUE),
            )
)  


## save spatial model
n.blocks = 2^L
filename=paste('multilevelIntS_gamma',n.blocks,'_',M,".RData",sep="")
save(multilevelS_gamma1, file = filename)

round(multilevelS_gamma1$summary.fixed[, c(1, 4,3, 5)] ,3)

round(multilevelS_gamma1$summary.hyperpar[, c(1, 4,3, 5)],3)

dic = multilevelS_gamma1$dic$dic
waic = multilevelS_gamma1$waic$waic
LPML = sum(log(multilevelS_gamma1$cpo$cpo))
time = sum(multilevelS_gamma1$bru_timings[,3])

y<- (train$totalPrice)
est0 <- multilevelS_gamma1$summary.fitted.values$mean[1:length(y)]
RMSE<- sqrt(mean((est0-y)^2))

marg.sig1<-inla.tmarginal(function(x){(exp(x))^2},multilevelS_gamma1$marginals.hyperpar$`Theta1 for spatintercept`)
inla.zmarginal(marg.sig1)
marg.phis1<-inla.tmarginal(function(x){2/exp(x)},multilevelS_gamma1$marginals.hyperpar$`Theta2 for spatintercept`)
inla.zmarginal(marg.phis1)
marg.ranges1<-inla.tmarginal(function(x){exp(x)},multilevelS_gamma1$marginals.hyperpar$`Theta2 for spatintercept`)
inla.zmarginal(marg.ranges1)

require(fields)
require(akima)
rh<-train.coord$Lng
rv<-train.coord$Lat
rz1<-multilevelS_gamma1$summary.random$spatintercept$mean

library(MBA)
int.elev1 <- mba.surf(cbind(rh,rv,rz1), 200, 200, extend=TRUE)$xyz.est


n.blocks <- 2^L
dir.save = getwd()
jpeg(paste(dir.save,'/Fulmultilevel_gammaInt_L',n.blocks,'_',M,"_etas.jpeg",sep=""))
image.plot(int.elev1 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
dev.off()


X <- as.matrix(cbind(1,train[,c("area","livingRoom","gamma_int_elev",
                                "gamma_int_subway", 
                                "gamma_int_age")]))


X1 <- data.frame(X[,-1])

y<- (train$totalPrice)
nzip<-dim(train.coord)[1]

n.sample <- 100 
set.seed(5)
multilevel0.pos <- generate(multilevelS_gamma1, X1, n.samples = n.sample)
multilevel0.pos.hyper <- inla.hyperpar.sample(n = n.sample, result=multilevelS_gamma1)

filename2= paste("gammaInt_posterior",n.blocks,"_",M,".Rdata",sep="")
save(multilevel0.pos, file=filename2)
filename2= paste("gammaInt_porsteriorprec",n.blocks,"_",M,".Rdata",sep="")
save(multilevel0.pos.hyper, file=filename2)


# for intercept process
predict.spatial1<- function(j="iterations",i="ith location"){
  
  #  Mp=M   
  theta1<-multilevel0.pos[[j]]$Theta1_for_spatintercept
  theta2<-multilevel0.pos[[j]]$Theta2_for_spatintercept
  
  sigma2<-exp(theta1)^2
  phis<-2/exp(theta2)
  
  b.i<-train.coord$block.ind[which.min(spDistsN1(as.matrix(train.coord[,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))]
  
  # Mp =1 or Mp=M
  if(Mp==1|b.i==1){
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% multilevel0.pos[[j]]$spatintercept[which(train.coord$block.ind==b.i)],sqrt(F.i))
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[b.i-1,])
    b.ii<- c(b.i, unique(nn.block$block.ind))
    indblockpred <- NULL
    for(k in 1:length(b.ii)){
      indp <- which(train.coord$block.ind==b.ii[k])
      indblockpred <- c(indblockpred,indp)
    }
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[indblockpred,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    D_b1<-spDists(as.matrix(train.coord[indblockpred, c("Lng","Lat")]),longlat = TRUE)
    C_N.inv<-solve(sigma2*exp(-phis*D_b1))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% multilevel0.pos[[j]]$spatintercept[indblockpred],sqrt(F.i))
  }
  return(z.i)
}


Mp=M

pre_Z.0<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.0)<-test.coord$location
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res1<-mclapply(1:n.sample,predict.spatial1,i=i,mc.cores = 4)
  pre_Z.0[,i]<-unlist(res1)
}



y<- (train$totalPrice)
nzip<-dim(train.coord)[1]

get.est<-function(i){
  precision<-multilevel0.pos.hyper[[i]]#$Precision_parameter_for_the_Gamma_observations
  coef <- matrix(unlist(multilevel0.pos[[i]][2:7]), 6,1)
  mu <- exp(X%*%coef+
              rep(multilevel0.pos[[i]]$spatintercept[1:nzip],train.coord$Freq))
  b         <- mu / precision
  res1<-rgamma(n=dim(X)[1], shape = precision ,  scale = b)
  return(res1)
}
y.rep1<-mclapply(1:n.sample,get.est,mc.cores = 4)
y.rep<-do.call(cbind,y.rep1)



X.test <- as.matrix(cbind(1,test[,c("area","livingRoom","gamma_int_elev",
                                    "gamma_int_subway", 
                                    "gamma_int_age")]))

y.test<-(test$totalPrice)


nzip<-dim(train.coord)[1]



get.pre<-function(i){
  precision<-multilevel0.pos.hyper[[i]]#$Precision_parameter_for_the_Gamma_observations
  coef <- matrix(unlist(multilevel0.pos[[i]][2:7]), 6,1)
  mu <- exp(X.test%*%coef+
              rep(pre_Z.0[i,],test.coord$Freq))
  b         <- mu / precision
  res2<-rgamma(n=dim(X.test)[1], shape = precision ,  scale = b)
  return(res2)
}

y.pre1<-mclapply(1:n.sample, get.pre, mc.cores = 4)
y.pre<-do.call(cbind, y.pre1)

pre.block.nngp0<-rowMeans(y.pre)

est.block.nngp<-rowMeans(y.rep)

RMSEy<- sqrt(mean((est.block.nngp-y)^2))

RMSP<- sqrt(mean((pre.block.nngp0-y.test)^2))

plot3 <- data.frame(location=test$location, totalPrice=test$totalPrice, 
                    Lng=test$Lng, Lat=test$Lat, pred=pre.block.nngp0 )

jpeg(paste(dir.save,"/PREDrandomint_gamma_",n.blocks, "-",M,"gamma.jpeg",sep=""))
ggplot(data = plot3,mapping = aes(x=totalPrice,y=pred))+geom_point()+ geom_abline(colour = "grey50")
dev.off()

# save final res
name.crit <- c("dic", "waic", "LPML", "RMSEy","time(sec)", "RMSPy")
crit.model<- round(c( dic, waic, LPML, RMSE, time, RMSP),3)
res=data.frame(name.crit, val=crit.model)
filename2= paste("randomint_gamma_summary_",n.blocks, "-",M,".Rdata",sep="")
save(res, file=filename2)



# random intercept process beta=w_s*gamma + eta_s
get.est5<-function(i){
  coef <- matrix(unlist(multilevel0.pos[[i]][2:7]), 6,1)
  ind <- c(1,4,5,6)
  mu <-  (X[,ind]%*%coef[ind]+ rep(multilevel0.pos[[i]]$spatintercept[1:nzip],train.coord$Freq))
  return(mu)
}
y.rep5<-mclapply(1:n.sample,get.est5,mc.cores = 4)
y.rep5<-do.call(cbind,y.rep5)
rz5<-rowMeans(y.rep5)

rh5<-train$Lng
rv5<-train$Lat
int.elev5 <- mba.surf(cbind(rh5,rv5,rz5), 200, 200, extend=TRUE)$xyz.est

jpeg(paste(dir.save,"/randomint_gamma_",n.blocks, "-",M,"gamma.jpeg",sep=""))
image.plot(int.elev5 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
dev.off()


################################################################
############## Multilevel spatial in all coeff ##############
###############################################################


rm(list=ls())

load("/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/Fast Spatial Multilevel Model (Spatial Statistics)/FINALapp_pcPrior/NEWbeijing.RData")

# Set the number of blocks L and the number of neighbors M
L=8 #2^L blocks
M=8 #number of neighbors

dir.save=getwd()

idx<-match(train$location,train.coord$location)
idx1<-match(test$location,test.coord$location)

train.coord$block.ind<-1
train.coord$ind<-1
for (i in 1:L) {
  if(i%%2==1){
    d<- 1
  }else{d<- 2}
  K<- 2^(i-1)
  for (k in 1:K) {
    cut.point<-median(train.coord[train.coord$ind==k ,d])
    train.coord[train.coord$ind==k &train.coord[,d]<=cut.point,"block.ind"]<- 2*k-1
    train.coord[train.coord$ind==k &train.coord[,d]>cut.point,"block.ind"]<-2*k
  }
  train.coord$ind<- train.coord$block.ind
}
train.coord$ind<-NULL
ggplot(data = train.coord,mapping = aes(x=Lng,y=Lat,color=as.factor(block.ind)))+geom_point()+
  geom_text(aes(label=as.factor(block.ind)),hjust=0, vjust=0)+ theme(legend.position = "none")


block.center<-cbind(aggregate(train.coord$Lng,list(train.coord$block.ind),mean)[,2],
                    aggregate(train.coord$Lat,list(train.coord$block.ind),mean)[,2],
                    1:2^L)
block.center<-as.data.frame(block.center)

NN_ind<-matrix(rep(0, (dim(block.center)[1]-1)*M),nrow = dim(block.center)[1]-1,ncol = M)
NN_distM<-matrix(rep(0, (dim(block.center)[1]-1)*M*(M-1)/2),nrow = dim(block.center)[1]-1,ncol = M*(M-1)/2)
NN_dist<-matrix(rep(0, (dim(block.center)[1]-1)*M),nrow = dim(block.center)[1]-1,ncol = M)

for (i in 2:dim(block.center)[1]) {
  #compute great circle distance
  distance<-spDistsN1(as.matrix(block.center[1:(i-1),1:2]),as.matrix(block.center[i,1:2]),longlat = T) 
  near.distance<-sort(distance,decreasing = FALSE,index.return=TRUE)
  dist.NN<-near.distance$x[1:min(i-1,M)] #dist to nearest neighbors
  ind.NN<-near.distance$ix[1:min(i-1,M)] #index of nearest neighbors
  distM.NN<-spDists(as.matrix(block.center[ind.NN,1:2]),longlat = T)
  
  NN_ind[i-1,1:min(i-1,M)]<-ind.NN
  if(i>=3){
    NN_distM[i-1,1:(min(i-1,M)*(min(i-1,M)-1)/2)]<-distM.NN[lower.tri(distM.NN,diag = FALSE)]
  }
  NN_dist[i-1,1:min(i-1,M)]<-dist.NN
}

neighbor.block<-list(NN_ind=NN_ind, NN_distM=NN_distM, NN_dist=NN_dist)



train.coord<-train.coord[order(train.coord$block.ind),]
train<-train[order(match(train$location,train.coord$location)),]

D_bN<-list()
D_NN<-list()
D_b<-list()
for (i in 1:2^L) {
  if(i==1){
    D_b[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),longlat = T)
    
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[i-1,])
    D_bN[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),as.matrix(nn.block[, c("Lng","Lat")]),longlat = T)
    D_NN[[i]]<-spDists(as.matrix(nn.block[, c("Lng","Lat")]),longlat = T)
    D_b[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),longlat = T)
  }
}


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
    res.BF<-mclapply(1:2^nb,build.BF,mc.cores = 4)
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
    
    Bs<-mclapply(1:2^nb,construct.Bs,mc.cores = 4)
    Bs<-do.call(cbind,Bs)
    
    F.inv<-bdiag(mclapply(1:2^nb,function(i){ solve(res.BF[[i]]$F)}))
    
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


library("INLA")


## Calculate hyperparameters
d=2
prior.range  = c(0.1,0.05)
prior.sigma  = c(1 ,0.05)
lam11 <- -log(prior.range[2])*prior.range[1]^(d/2)
lam22 <- -log(prior.sigma[2])/prior.sigma[1]

initial.range <- log(prior.range[1]) + 1
initial.sigma <- log(prior.sigma[1]) - 1


block.nngp.model <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                         nb=L, D_bN=D_bN, D_NN=D_NN, D_b=D_b, lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma, debug = T)

block.nngp.model1 <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                          nb=L, D_bN=D_bN, D_NN=D_NN, D_b=D_b,lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma,debug = T)

block.nngp.model2 <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                          nb=L, D_bN=D_bN, D_NN=D_NN, D_b=D_b,lam1=lam11, lam2=lam22,initial.range=initial.range,initial.sigma=initial.sigma,debug = T)



svc_formula1 <- totalPrice ~ -1 + intercept +  area  + livingRoom +
  gamma_int_elev +  gamma_area_elev +  gamma_livingRoom_elev +
  gamma_int_subway + gamma_area_subway   + gamma_livingRoom_subway+
  gamma_int_age + gamma_area_age   + gamma_livingRoom_age +
  spatintercept(idx, model = block.nngp.model) +
  spatareacov(idx, weights = area, model = block.nngp.model1) +
  spatlivingRoomcov(idx, weights = livingRoom, model = block.nngp.model2)


multilevel1 <- bru( svc_formula1,
                     family = "gamma",
                     data = train,
                   options = list(
                     control.compute = list(waic = TRUE, cpo = TRUE,config=TRUE),
                     verbose = FALSE
                   )
)  

## save non-spatial model
n.blocks = 2^L
save(multilevel1, file = paste(dir.save,"/multilevelSpat_gamma_",n.blocks, "-",M,".RData",sep=""))


round(multilevel1$summary.hyperpar[, c(1, 4, 3, 5)] ,3) 

round(multilevel1$summary.fixed[, c(1, 4, 3, 5)] ,3) 

dic = multilevel1$dic$dic
waic = multilevel1$waic$waic
LPML = sum(log(multilevel1$cpo$cpo))
time = sum(multilevel1$bru_timings[,3])


est1 <- multilevel1$summary.fitted.values

y<- (train$totalPrice)
nzip<-dim(train.coord)[1]


y<- (train$totalPrice)
est0 <- multilevel1$summary.fitted.values$mean[1:length(y)]
RMSE<- sqrt(mean((est0-y)^2))


plot2 <- data.frame(location=train$location, totalPrice=train$totalPrice, 
                    Lng=train$Lng, Lat=train$Lat, pred=fit[,1] )
n.blocks <- 2^L
jpeg(paste(dir.save,"/est_gamma_multilevel1", n.blocks, "-",M,".jpeg",sep=""))
ggplot(data = plot2,mapping = aes(x=totalPrice,y=pred))+geom_point()+
  geom_abline(colour = "grey50")+ 
  labs(x = "y", y = "y_est")
dev.off()


marg.sig<-inla.tmarginal(function(x){(exp(x))^2},multilevel1$marginals.hyperpar$`Theta1 for spatintercept`)
inla.zmarginal(marg.sig)
marg.phis<-inla.tmarginal(function(x){2/exp(x)},multilevel1$marginals.hyperpar$`Theta2 for spatintercept`)
inla.zmarginal(marg.phis)
marg.ranges<-inla.tmarginal(function(x){exp(x)},multilevel1$marginals.hyperpar$`Theta2 for spatintercept`)
inla.zmarginal(marg.ranges)

marg.sig2<-inla.tmarginal(function(x){exp(x)^2},multilevel1$marginals.hyperpar$`Theta1 for spatareacov`)
inla.zmarginal(marg.sig2)
marg.phis2<-inla.tmarginal(function(x){2/exp(x)},multilevel1$marginals.hyperpar$`Theta2 for spatareacov`)
inla.zmarginal(marg.phis2)
marg.ranges2<-inla.tmarginal(function(x){exp(x)},multilevel1$marginals.hyperpar$`Theta2 for spatareacov`)
inla.zmarginal(marg.ranges2)

marg.sig3<-inla.tmarginal(function(x){exp(x)^2},multilevel1$marginals.hyperpar$`Theta1 for spatlivingRoomcov`)
inla.zmarginal(marg.sig3)
marg.phis3<-inla.tmarginal(function(x){2/exp(x)},multilevel1$marginals.hyperpar$`Theta2 for spatlivingRoomcov`)
inla.zmarginal(marg.phis3)
marg.ranges3<-inla.tmarginal(function(x){exp(x)},multilevel1$marginals.hyperpar$`Theta2 for spatlivingRoomcov`)
inla.zmarginal(marg.ranges3)


require(fields)
require(akima)
rh<-train.coord$Lng
rv<-train.coord$Lat
rz1<-multilevel1$summary.random$spatintercept$mean

library(MBA)
int.elev1 <- mba.surf(cbind(rh,rv,rz1), 200, 200, extend=TRUE)$xyz.est

rz2<-multilevel1$summary.random$spatareacov$mean
int.elev2 <- mba.surf(cbind(rh,rv,rz2), 200, 200, extend=TRUE)$xyz.est

rz3<-multilevel1$summary.random$spatlivingRoomcov$mean
int.elev3 <- mba.surf(cbind(rh,rv,rz3), 200, 200, extend=TRUE)$xyz.est

n.blocks <- 2^L
dir.save = getwd()
jpeg(paste(dir.save,'/Fulmultilevel_gamma_L',n.blocks,'_',M,"_.jpeg",sep=""))
par(mfrow=c(1,3))
image.plot(int.elev1 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev2 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev3 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
dev.off()


X <- as.matrix(cbind(1,train[,c("area","livingRoom","gamma_int_elev","gamma_area_elev",
                                "gamma_livingRoom_elev","gamma_int_subway","gamma_area_subway",
                                "gamma_livingRoom_subway", 
                                 "gamma_int_age", "gamma_area_age", 
                                "gamma_livingRoom_age")]))


X1 <- data.frame(X[,-1])

y<- (train$totalPrice)
nzip<-dim(train.coord)[1]

n.sample <- 100#1e4
set.seed(5)
multilevel1.pos <- generate(multilevel1, X1, n.samples = n.sample)
multilevel1.pos.hyper <- inla.hyperpar.sample(n = n.sample, result=multilevel1)

filename2= paste("FULLMULTILEVELrandomcoef_gamma_posterior_",n.blocks, "-",M,".Rdata",sep="")
save(multilevel1.pos, file=filename2)
filename2= paste("FULLMULTILEVELrandomcoef_gamma_posteriorprec_",n.blocks, "-",M,".Rdata",sep="")
save(multilevel1.pos.hyper, file=filename2)


Mp = M

# for intercept process
predict.spatial1<- function(j="iterations",i="ith location"){
  
  theta1<-multilevel1.pos[[j]]$Theta1_for_spatintercept
  theta2<-multilevel1.pos[[j]]$Theta2_for_spatintercept
  
  sigma2<-exp(theta1)^2
  phis<-2/exp(theta2)
  
   b.i<-train.coord$block.ind[which.min(spDistsN1(as.matrix(train.coord[,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))]
  
  # Mp =1 or Mp=M
  if(Mp==1|b.i==1){
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% multilevel1.pos[[j]]$spatintercept[which(train.coord$block.ind==b.i)],sqrt(F.i))
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[b.i-1,])
    b.ii<- c(b.i, unique(nn.block$block.ind))
    indblockpred <- NULL
    for(k in 1:length(b.ii)){
      indp <- which(train.coord$block.ind==b.ii[k])
      indblockpred <- c(indblockpred,indp)
    }
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[indblockpred,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    D_b1<-spDists(as.matrix(train.coord[indblockpred, c("Lng","Lat")]),longlat = TRUE)
    C_N.inv<-solve(sigma2*exp(-phis*D_b1))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% multilevel1.pos[[j]]$spatintercept[indblockpred],sqrt(F.i))
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
  
  theta1<-multilevel1.pos[[j]]$Theta1_for_spatareacov
  theta2<-multilevel1.pos[[j]]$Theta2_for_spatareacov
  
  sigma2<-exp(theta1)^2
  phis<-2/exp(theta2)
  
  b.i<-train.coord$block.ind[which.min(spDistsN1(as.matrix(train.coord[,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))]
  
  # Mp =1 or Mp=M
  if(Mp==1|b.i==1){
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% multilevel1.pos[[j]]$spatareacov[which(train.coord$block.ind==b.i)],sqrt(F.i))
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[b.i-1,])
    b.ii<- c(b.i, unique(nn.block$block.ind))
    indblockpred <- NULL
    for(k in 1:length(b.ii)){
      indp <- which(train.coord$block.ind==b.ii[k])
      indblockpred <- c(indblockpred,indp)
    }
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[indblockpred,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    D_b1<-spDists(as.matrix(train.coord[indblockpred, c("Lng","Lat")]),longlat = TRUE)
    C_N.inv<-solve(sigma2*exp(-phis*D_b1))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% multilevel1.pos[[j]]$spatareacov[indblockpred],sqrt(F.i))
  }
  return(z.i)
}



pre_Z.2<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.2)<-test.coord$location
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res1<-mclapply(1:n.sample,predict.spatial2,i=i,mc.cores = 4)
  pre_Z.2[,i]<-unlist(res1)
}

#for b2 process
predict.spatial3<- function(j="iterations",i="ith location"){
  
  theta1<-multilevel1.pos[[j]]$Theta1_for_spatlivingRoomcov
  theta2<-multilevel1.pos[[j]]$Theta2_for_spatlivingRoomcov
  
  sigma2<-exp(theta1)^2
  phis<-2/exp(theta2)
  
  b.i<-train.coord$block.ind[which.min(spDistsN1(as.matrix(train.coord[,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))]
  
  # Mp =1 or Mp=M
  if(Mp==1|b.i==1){
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
     z.i<-rnorm(1, B.i%*% multilevel1.pos[[j]]$spatlivingRoomcov[which(train.coord$block.ind==b.i)],sqrt(F.i))
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[b.i-1,])
    b.ii<- c(b.i, unique(nn.block$block.ind))
    indblockpred <- NULL
    for(k in 1:length(b.ii)){
      indp <- which(train.coord$block.ind==b.ii[k])
      indblockpred <- c(indblockpred,indp)
    }
    C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[indblockpred,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
    D_b1<-spDists(as.matrix(train.coord[indblockpred, c("Lng","Lat")]),longlat = TRUE)
    C_N.inv<-solve(sigma2*exp(-phis*D_b1))
    B.i<-C_iN %*% C_N.inv
    F.i<-sigma2-B.i%*%C_iN
    z.i<-rnorm(1, B.i%*% multilevel1.pos[[j]]$spatlivingRoomcov[indblockpred],sqrt(F.i))
  }
  return(z.i)
}



pre_Z.3<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.3)<-test.coord$location
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res1<-mclapply(1:n.sample,predict.spatial3,i=i,mc.cores = 4)
  pre_Z.3[,i]<-unlist(res1)
}


X <- as.matrix(cbind(1,train[,c("area","livingRoom","gamma_int_elev","gamma_area_elev",
                                "gamma_livingRoom_elev","gamma_int_subway","gamma_area_subway",
                                "gamma_livingRoom_subway", "gamma_int_age", "gamma_area_age", 
                                "gamma_livingRoom_age")]))

y<- (train$totalPrice)
nzip<-dim(train.coord)[1]

get.est<-function(i){
  precision<-multilevel1.pos.hyper[[i]]
  coef <- matrix(unlist(multilevel1.pos[[i]][2:13]), 12,1)
  mu <- exp(X%*%coef+
              rep(multilevel1.pos[[i]]$spatintercept[1:nzip],train.coord$Freq)+
              +train$area*rep(multilevel1.pos[[i]]$spatareacov[1:nzip],train.coord$Freq)+
              +train$livingRoom*rep(multilevel1.pos[[i]]$spatlivingRoomcov[1:nzip],train.coord$Freq))
  b         <- mu / precision
  res1<-rgamma(n=dim(X)[1], shape = precision ,  scale = b)
  return(res1)
}
y.rep<-mclapply(1:n.sample,get.est,mc.cores = 4)
y.rep<-do.call(rbind,y.rep)


X.test <- as.matrix(cbind(1,test[,c("area","livingRoom","gamma_int_elev","gamma_area_elev",
                                    "gamma_livingRoom_elev","gamma_int_subway","gamma_area_subway",
                                    "gamma_livingRoom_subway", "gamma_int_age", "gamma_area_age", 
                                    "gamma_livingRoom_age")]))

y.test<-(test$totalPrice)
nzip<-dim(train.coord)[1]

#test$totalPrice 
rh2<-test.coord$Lng
rv2<-test.coord$Lat
rz11<-colMeans(pre_Z.1)
rz22<-colMeans(pre_Z.2)
rz33<-colMeans(pre_Z.3)


get.pre<-function(i){
  precision<- multilevel1.pos.hyper[[i]]
  coef <- matrix(unlist(multilevel1.pos[[i]][2:13]), 12,1)
  mu <- exp(X.test%*%coef+
              rep(rz11,test.coord$Freq)+
              test$area*rep(rz22,test.coord$Freq)+
              test$livingRoom*rep(rz33,test.coord$Freq)
  )
  b         <- mu / precision
  res2<-rgamma(n=dim(X.test)[1], shape = precision ,  scale = b)
  return(res2)
}

y.pre1<-mclapply(1:n.sample, get.pre, mc.cores = 4)
y.pre<-do.call(rbind, y.pre1)

pre.block.nngp2<-colMeans(y.pre)

est.block.nngp<-colMeans(y.rep)

RMSEy<- sqrt(mean((est.block.nngp-y)^2))

RMSP<- sqrt(mean((pre.block.nngp2-y.test)^2))

plot3 <- data.frame(location=test$location, totalPrice=test$totalPrice, 
                    Lng=test$Lng, Lat=test$Lat, pred=pre.block.nngp2 )


#Fulmultilevel_gamma_L256_4_PRED
jpeg(paste(dir.save,"/Fulmultilevel_gamma_L",n.blocks,"_",M,"_PRED.jpeg",sep=""))
ggplot(data = plot3,mapping = aes(x=totalPrice,y=pred))+geom_point()+ geom_abline(colour = "grey50")
dev.off()

# save final res
name.crit <- c("dic", "waic", "LPML", "RMSEy","time(sec)", "RMSPy")
crit.model<- round(c( dic, waic, LPML, RMSE, time, RMSP),3)
res=data.frame(name.crit, val=crit.model)
filename2= paste("FULLMULTILEVELrandomcoef_gamma_summary_",n.blocks, "-",M,".Rdata",sep="")
save(res, file=filename2)



# random intercept process beta=w_s*gamma + eta_s
get.est5<-function(i){
  coef <- matrix(unlist(multilevel1.pos[[i]][2:13]), 12,1)
  ind <- c(1, 4, 7, 10)
  mu <-  (X[,ind]%*%coef[ind]+ rep(multilevel1.pos[[i]]$spatintercept[1:nzip],train.coord$Freq))
  return(mu)
}
y.rep5<-mclapply(1:n.sample,get.est5,mc.cores = 4)
y.rep5<-do.call(cbind,y.rep5)
rz5<-rowMeans(y.rep5)

rh5<-train$Lng
rv5<-train$Lat
int.elev5 <- mba.surf(cbind(rh5,rv5,rz5), 200, 200, extend=TRUE)$xyz.est

# random beta1 process beta1=w_s*gamma + eta_s
get.est51<-function(i){
  coef <- matrix(unlist(multilevel1.pos[[i]][2:13]), 12,1)
  ind <- c(2, 5, 8, 11)
  mu <- (X[,ind]%*%coef[ind]+ train$area*rep(multilevel1.pos[[i]]$spatareacov[1:nzip],train.coord$Freq))
  return(mu)
}
y.rep51<-mclapply(1:n.sample,get.est51,mc.cores = 4)
y.rep51<-do.call(cbind,y.rep51)
rz51<-rowMeans(y.rep51)

int.elev51 <- mba.surf(cbind(rh5,rv5,rz51), 200, 200, extend=TRUE)$xyz.est

# random beta2 process beta2=w_s*gamma + eta_s
get.est52<-function(i){
  coef <- matrix(unlist(multilevel1.pos[[i]][2:13]), 12,1)
  ind <- c(3, 6, 9, 12)
  mu <- (X[,ind]%*%coef[ind]+ train$livingRoom*rep(multilevel1.pos[[i]]$spatlivingRoomcov[1:nzip],train.coord$Freq))
  return(mu)
}
y.rep52<-mclapply(1:n.sample,get.est52,mc.cores = 4)
y.rep52<-do.call(cbind,y.rep52)
rz52<-rowMeans(y.rep52)


int.elev52 <- mba.surf(cbind(rh5,rv5,rz52), 200, 200, extend=TRUE)$xyz.est


jpeg(paste(dir.save,"/FULLMULTILEVELrandomcoef_gamma_",n.blocks, "-",M,"gamma.jpeg",sep=""))
par(mfrow=c(1,3))
image.plot(int.elev5 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev51 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev52 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
dev.off()



jpeg(paste(dir.save,"/FULLMULTILEVEL_plot10gamma_",n.blocks, "-",M,".jpeg",sep=""))
par(mfrow=c(2,3))
image.plot(int.elev1 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev2 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev3 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev5 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev51 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
image.plot(int.elev52 , #main='spatial random effects',
           xaxs = 'r', yaxs = 'r',
           xlab='Longitude', ylab='Latitude')
dev.off()


