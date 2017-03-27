
## Mack

# Script to run the Mack Model on data from the CAS Loss Reserve Database
# by Glenn Meyers
# This script uses the "MackChainLadder" function as it existed on 9/1/2014
#
rm(list = ls())      # clear workspace"
setwd("F:/Github - LabtopLee/Seminar-2017/week 2/triangle data")
#
# user inputs
#
insurer.data = "comauto_pos.csv"
#insurer.data = "ppauto_pos.csv"
#insurer.data = "wkcomp_pos.csv"
#insurer.data = "othliab_pos.csv"
#insurer.data = "prodliab_pos.csv"
#insurer.data = "medmal_pos.csv"
grpcode = "353"
losstype = "incloss"  #"incloss" if incurred loss or "cpdloss" if paid loss
a = read.csv(insurer.data)
library(data.table)
a = data.table(a)
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
  single=b[,12]
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
#
cdata=ins.line.data(grpcode)
w=cdata$w-1987
d=cdata$d
premium=cdata$net_premium
cpdloss=cdata$cum_pdloss
incloss=cdata$cum_incloss-cdata$bulk_loss
cpdloss=pmax(1,cpdloss)
incloss=pmax(1,incloss)
adata=data.frame(grpcode,w,d,premium,cpdloss,incloss)
rdata=subset(adata,(adata$w+adata$d)<12)
maxlevel=log(2*max(rdata$premium))
numw=length(unique(rdata$w))
if(losstype=="incloss") rloss=rdata$incloss else rloss=rdata$cpdloss
if(losstype=="incloss") aloss=adata$incloss else aloss=adata$cpdloss
#
# run Mack model using ChainLadder
#
library(ChainLadder)
# 将数据框形式转化为流量三角形形式
rtriangle = as.triangle(rdata,origin="w",dev="d",value=losstype)     
mcl=MackChainLadder(rtriangle,est.sigma="Mack")
#
# Calculate summary statistic for Mack model
#
Mack.Estimate=round(summary(mcl)$ByOrigin[,3])
Mack.S.E=round(summary(mcl)$ByOrigin[,5])
Mack.CV=round(Mack.S.E/Mack.Estimate,4)
Mack.total.ult=round(sum(summary(mcl)$ByOrigin[1:10,3]))
Mack.total.se=round(summary(mcl)$Totals[5,1])
Mack.total.cv=round(Mack.total.se/Mack.total.ult,4)
Mack.sig=sqrt(log(1+Mack.total.cv^2))
Mack.mu=log(Mack.total.ult)-Mack.sig^2/2
act10d=subset(aloss,adata$d==10)
acttot=sum(act10d)
pct.Mack=round(plnorm(acttot,Mack.mu,Mack.sig)*100,2)
W=c(1:10,"Total")
Mack.Estimate=c(Mack.Estimate,Mack.total.ult)
Mack.S.E=c(Mack.S.E,Mack.total.se)
Mack.CV=c(Mack.CV,Mack.total.cv)
Actual=c(act10d,acttot)
Percentile=c(rep("",10),pct.Mack)
risk=data.frame(W,Mack.Estimate,Mack.S.E,Mack.CV,Actual,Percentile)
print(risk)
#
# write files for the monograph
#
write.csv(risk,file="Mack Ill Incurrred.csv",row.names=F)
atriangle=as.triangle(adata,origin="w",dev="d",value=losstype)
write.csv(atriangle,file="Illustrated Incurred.csv",row.names=F)

