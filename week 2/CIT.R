
## CIT

#
# Script to run the CIT Model on data from the CAS Loss Reserve Database
# To run the LIT model, set prior for rho ~ dunif(-.00001,00001) in JAGS script
# by Glenn Meyers
#
rm(list = ls())      # clear workspace"

#
# user inputs
#
insurer.data="comauto_pos.csv"
#insurer.data="ppauto_pos.csv"
#insurer.data="wkcomp_pos.csv"
#insurer.data="othliab_pos.csv"
#insurer.data="prodliab_pos.csv"
#insurer.data="medmal_pos.csv"
grpcode="353"
outfile="outCIT.csv"
nburn=50000
#
# JAGS script
#
modelString = "model {
mu[1]<-alpha[w[1]]+beta[d[1]]+tau*(w[1]+d[1]-1)   
logz[1]~dnorm(mu[1],1/sig2[1])
z[1]<-exp(logz[1])
loss[1]~dnorm(z[1],deltam2)  
for(i in 2:length(w)){
mu[i]<-alpha[w[i]]+beta[d[i]]+tau*(w[i]+d[i]-1)
logz[i]~dnorm(mu[i],1/sig2[i])
z[i]<-exp(logz[i])+rho*(loss[i-1]-z[i-1])*exp(tau)*wne1[i]
loss[i]~dnorm(z[i],deltam2)
}
#
# set up sig2
#
for (i in 1:length(w)){
sig2[i]<-sigd2[d[i]]
}
sigd2[1]~dunif(.000001,0.5)
for (j in 2:10){
sigd2[j]~dunif(sigd2[j-1],sigd2[j-1]+.1)         # control growth of sigma
}
#
# specify priors
#
for (i in 1:numlev){
alpha[i]~dnorm(log(premium[i])+logelr,0.1)       # std dev of alpha = 1/sqrt(.1) = 3.16
}
for (i in 1:4){
beta[i]~dunif(0,10)
}
for (i in 5:9){
beta[i]~dunif(0,beta[i-1])                       # force beta to decrease for d > 4
}
beta[10]<-0
rho~dunif(-1,1)
# rho~dunif(-.00001,.00001) # Use for LIT model
logelr~dunif(-5,1)
tau~dnorm(0,1000)                                    # std dev of tau = 1/sqrt(1000) = 0.0316
delta~dunif(0,sum(premium)/10)
deltam2<-1/delta^2
}" 
#
# get data
#
a=read.csv(insurer.data)
#
# function to get Schedule P triangle data given ins group and line of business
#
ins.line.data=function(g.code){
  b=subset(a,a$GRCODE==g.code)
  name=b$GRNAME
  grpcode=b$GRCODE
  w=b$AccidentYear
  d=b$DevelopmentLag
  cum_incloss=b[,6]
  cum_pdloss=b[,7]
  bulk_loss=b[,8]
  dir_premium=b[,9]
  ced_premium=b[,10]
  net_premium=b[,11]
  single=b$Single
  posted_reserve97=b[,13]  
  # get incremental paid losses - assume data is sorted by ay and lag
  inc_pdloss=numeric(0)
  for (i in unique(w)){
    s=(w==i)
    pl=c(0,cum_pdloss[s])
    ndev=length(pl)-1  
    il=rep(0,ndev)
    for (j in 1:ndev){            
      il[j]=pl[j+1]-pl[j]
    }
    inc_pdloss=c(inc_pdloss,il)
  }
  data.out=data.frame(grpcode,w,d,net_premium,dir_premium,ced_premium,
                      cum_pdloss,cum_incloss,bulk_loss,inc_pdloss,single,posted_reserve97)
  return(data.out)
}
#
# read and aggregate the insurer data and 
# set up training and test data frames
cdata=ins.line.data(grpcode)
set.seed(12345)
w=cdata$w-1987
d=cdata$d
# 
# sort the data in order of d, then w within d
#
o1=100*d+w
o=order(o1)
w=w[o]
d=d[o]
premium=cdata$net_premium[o]
ipdloss=cdata$inc_pdloss[o]
cpdloss=cdata$cum_pdloss[o]
wne1=ifelse(w==1,0,1)
adata=data.frame(grpcode,w,d,premium,ipdloss,cpdloss,wne1)
rdata=subset(adata,(adata$w+adata$d)<12)
numw=length(unique(rdata$w))
rdata=subset(adata,(adata$w+adata$d)<12)
rloss=rdata$ipdloss
aloss=adata$ipdloss
Premium=rdata$premium[1:10]
#
# Initialize JAGS model
#
inits1=list(.RNG.name= "base::Wichmann-Hill",
            .RNG.seed= 12341)
inits2=list(.RNG.name= "base::Marsaglia-Multicarry",
            .RNG.seed= 12342)
inits3=list(.RNG.name= "base::Super-Duper",
            .RNG.seed= 12343)
inits4=list(.RNG.name= "base::Mersenne-Twister",
            .RNG.seed= 12344)
data.for.jags=list('premium'= Premium,
                   'loss'   = rloss,
                   'numlev' = numw,
                   'w'      = rdata$w,
                   'wne1'   = rdata$wne1,
                   'd'      = rdata$d)
#
# run the model
#
library(runjags)
library(coda)
nthin=2
maxpsrf=2
while (maxpsrf>1.05){
  nthin=nthin*2
  print(paste("nthin =",nthin))
  jagout=run.jags(model=modelString,
                  monitor=c("alpha","beta[1:9]","sigd2","rho","delta","tau","logelr"),
                  data=data.for.jags,n.chains=4,method="parallel",
                  inits=list(inits1,inits2,inits3,inits4),thin=nthin,silent.jags=F,
                  plots=TRUE,burnin=nburn,sample=2500,psrf.target=1.05)
  gelman=gelman.diag(jagout)
  maxpsrf=max(gelman$psrf[,1])
  print(paste("maxpsrf =",maxpsrf))
} 
#
# optional diagnostics
#
library(coda)
#crosscorr.plot(jagout$mcmc)
# traceplot(jagout$mcmc)
# print(gelman.diag(jagout))
# gelman.plot(jagout)
print(summary(jagout)[[1]][,1:2])
#
# extract information from jags output to process in R
#
b=as.matrix(jagout$mcmc)
alpha=b[,1:10]
beta=cbind(b[,11:19],rep(0,dim(b)[1]))
delta=b[,20]
logelr=b[,21]
rho=b[,22]
sigd2=b[,23:32]
tau=b[,33]
#
# simulate loss statistics by accident year reflecting total risk for the JAGS model
#
set.seed(12345)
at.wd10=matrix(0,dim(b)[1],10)
for (w in 1:10){
  latest=sum(subset(rdata$ipdloss,rdata$w==w))
  at.wd10[,w]=rep(latest,dim(b)[1])
}
library(ChainLadder)
itriangle=as.triangle(rdata,origin="w",dev="d",value="ipdloss")
z=matrix(0,10,10)
logz=matrix(0,10,10)
for (i in 1:dim(b)[1]){
  m=matrix(0,10,10)
  for (d in 2:9){
    w=1
    m[w,d]=alpha[i,w]+beta[i,d]+tau[i]*(w+d-1)
    for (w in 2:(11-d)){
      m[w,d]=alpha[i,w]+beta[i,d]+tau[i]*(w+d-1)
      z[w,d]=exp(m[w,d])+rho[i]*(itriangle[w-1,d]-z[w-1,d])*exp(tau[i])
    }
  }
  d=10
  w=1
  m[w,d]=alpha[i,w]+beta[i,d]+tau[i]*(w+d-1)
  #
  for (d in 2:10){
    for (w in (12-d):10){
      m[w,d]=alpha[i,w]+beta[i,d]+tau[i]*(w+d-1)
      logz[w,d]=rnorm(1,m[w,d],sigd2[i,d])
      z[w,d]=exp(m[w,d])+
        rho[i]*(itriangle[w-1,d]-z[w-1,d])*exp(tau[i])
      itriangle[w,d]=rnorm(1,z[w,d],delta[i])
      at.wd10[i,w]=at.wd10[i,w]+itriangle[w,d]
    }
  }
}
ss.wd10=rep(0,10)
ms.wd10=rep(0,10)
ms.wd10[1]=mean(at.wd10[,1])
for (w in 2:10){
  ms.wd10[w]=mean(at.wd10[,w])
  ss.wd10[w]=sd(at.wd10[,w])
}
Pred.CIT=rowSums(at.wd10)
ms.td10=mean(Pred.CIT)
ss.td10=sd(Pred.CIT)
CIT.Estimate=round(ms.wd10)
CIT.SE=round(ss.wd10)
CIT.CV=round(CIT.SE/CIT.Estimate,4)
Outcome=subset(adata$cpdloss,adata$d==10)[1:10]
CIT.Pct=sum(Pred.CIT<=sum(Outcome))/length(Pred.CIT)*100
#
# put accident year statistics into a data frame
#
W=c(1:10,"Total")
CIT.Estimate=c(CIT.Estimate,round(ms.td10))
CIT.SE=c(CIT.SE,round(ss.td10))
CIT.CV=c(CIT.CV,round(ss.td10/ms.td10,4))
Premium=c(Premium,sum(Premium))
Outcome=c(Outcome,sum(Outcome))
Group=rep(grpcode,11)
CIT.Pct=c(rep(NA,10),CIT.Pct)
risk=data.frame(W,Premium,CIT.Estimate,CIT.SE,CIT.CV,Outcome,CIT.Pct)
print(risk)
write.csv(risk,file=outfile,row.names=F)
par(mfrow=c(2,1))
hist(rho,main="",xlab=expression(rho))
hist(tau,main="",xlab=expression(tau))
par(mfrow=c(1,1))
hist(Pred.CIT,main="Predictive Distribution of Outcomes",xlab="")
#
# Plot to show some sample correlations See the "crosscorr.plot output above.
#
par(mfrow=c(2,2))
plot(alpha[,6]+beta[,1],tau,xlab=expression(paste(alpha[6],"+",beta[1])),
     ylab=expression(tau),pch=19)
plot(alpha[,6]+beta[,3],tau,xlab=expression(paste(alpha[6],"+",beta[3])),
     ylab=expression(tau),pch=19)
plot(alpha[,6]+beta[,5],tau,xlab=expression(paste(alpha[6],"+",beta[5])),
     ylab=expression(tau),pch=19)
plot(alpha[,6]+beta[,7],tau,xlab=expression(paste(alpha[6],"+",beta[7])),
     ylab=expression(tau),pch=19)

