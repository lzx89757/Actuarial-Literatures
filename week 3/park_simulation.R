## GPD拟合尾部样本，尺度和形状参数估计
myparaesti <- function(x,q,scale.ini=0.1,shape.ini=0.01){
  # 输入：x 样本
  #       q 阈值对应分位数
  #       scale.ini 尺度参数初始值
  #       shape.ini 形状参数初始值
  # 输出：尺度和形状参数估计值
  
  n <- length(x)
  u <- quantile(x,q) # 阈值
  x.u <- x[x>u]
  x.u <- sort(x.u,decreasing = TRUE) # 降序排列尾部样本
  
  # GPD分布函数
  mygpdf <- function(x,u,sigma,chsi){
    res <- 1-(1+chsi*(x-u)/sigma)^(-1/chsi)
    return(res)
  }
  
  # 截断经验分布函数值(降序)
  myecdf <- function(n,q){
    res <- (seq(1,(q*n+1)/n,-1/n)-q)/(1-q)
    return(res)
  }
  
  # 权函数值(降序)
  myweight <- function(n,q){
    res <- (n+1)^2*(n+2)/(seq(1, n*(1-q),1)*seq(n, n*q+1,-1))
    return(res)
  }  
  
  # 初始值估计
  rss1 <- function(param,x.u){
    res <- sum((log(1-myecdf(n,q)+1/(10*n))-log(1-mygpdf(x.u,u,param[1],param[2])+1/(10*n)))^2) # 避免出现log(0)
    return(res)
  }
  res.ini <- optim(c(scale.ini,shape.ini),rss1,x.u=x.u)$par
  
  # pot-NLS 估计
  rss2 <- function(param,x.u){
    res <- sum((myecdf(n,q)-mygpdf(x.u,u,param[1],param[2]))^2)
    return(res)
  }
  res.nls <- optim(res.ini,rss2,x.u=x.u)$par
  
  # pot-WNLS 估计
  rss3 <- function(param,x.u){
    res <- sum(myweight(n,q)*(myecdf(n,q)-mygpdf(x.u,u,param[1],param[2]))^2)
    return(res)
  }
  res.wnls <- optim(res.ini,rss3,x.u=x.u)$par 
  
  return(c(res.ini,res.nls,res.wnls))
}

## VaR估计式
myvar <- function(sigma,chsi,u,nu,p){
  res <- u+sigma/chsi*((n/nu*(1-p))^(-chsi)-1)
  return(res)
}

## CTE估计式
mycte <- function(sigma,chsi,u,nu,p){
  res <- u+sigma/chsi*(1/(1-chsi)*(n*(1-p)/nu)^(-chsi)-1)
  return(res)
}

####################################### 模拟设置
nn <- 10^7                            # 总体总量
n <- 10000	                          # 样本量
sim.n <- 200	                        # 重复模拟次数
q <- 0.95                             # 阈值对应分位数 0.95 0.98 0.99
nu <- n*(1-q)	                        # 拟合样本量
p <- c(0.98,0.99,0.999,0.9999)        # 需估计分位数
library(POT)
sigma <- 1                            # GPD尺度参数
chsi <- 0.75                          # GPD形状参数
xx <- rgpd(nn,0,sigma,chsi)           # 产生GPD总体  
VaR.real <- qgpd(p,0,sigma,chsi)      # GPD的真实VaR
VaR.real.mat <-matrix(VaR.real,nrow=sim.n,ncol=length(p),byrow=T)
CTE.real <- (sigma+VaR.real)/(1-chsi) # GPD的真实CTE
CTE.real.mat <-matrix(CTE.real,nrow=sim.n,ncol=length(p),byrow=T)

VaR.NLS <- matrix(0,nrow=sim.n,ncol=length(p))
VaR.WNLS <- matrix(0,nrow=sim.n,ncol=length(p))
VaR.pick <- matrix(0,nrow=sim.n,ncol=length(p))

CTE.NLS <- matrix(0,nrow=sim.n,ncol=length(p))
CTE.WNLS <- matrix(0,nrow=sim.n,ncol=length(p))
CTE.pick <- matrix(0,nrow=sim.n,ncol=length(p))

failure <- matrix(0,nrow=1,ncol=3)

for (j in 1:sim.n){
  # 有放回抽样
  index <- floor(runif(n,1,nn))
  x <- xx[index]
  
  # 估计GPD参数
  u <- quantile(x,q) # 阈值
  temp <- myparaesti(x,q)
  para.NLS <- temp[3:4]
  para.WNLS <- temp[5:6]
  para.pick <- fitgpd(x,u,"pickands")$param
  failure <- failure+c(para.NLS[2]>0.99,para.WNLS[2]>0.99,para.pick[2]>0.99) # 失效次数

  # 估计VaR
  VaR.NLS[j,] <- myvar(para.NLS[1],para.NLS[2],u,nu,p)
  VaR.WNLS[j,] <- myvar(para.WNLS[1],para.WNLS[2],u,nu,p)
  VaR.pick[j,] <- myvar(para.pick[1],para.pick[2],u,nu,p)
  
  # 估计CTE
  CTE.NLS[j,] <- mycte(para.NLS[1],para.NLS[2],u,nu,p)
  CTE.WNLS[j,] <- mycte(para.WNLS[1],para.WNLS[2],u,nu,p)
  CTE.pick[j,] <- mycte(para.pick[1],para.pick[2],u,nu,p)
  
  print(j)
}

# VaR 估计 RMSE
rmse.VaR <- rbind(sqrt(colMeans((VaR.NLS-VaR.real.mat)^2)),sqrt(colMeans((VaR.WNLS-VaR.real.mat)^2)),
              sqrt(colMeans((VaR.pick-VaR.real.mat)^2)))

# VaR 估计 ARB
arb.VaR <- rbind(colMeans(abs(VaR.NLS-VaR.real.mat))/VaR.real,colMeans(abs(VaR.WNLS-VaR.real.mat))/VaR.real,
             colMeans(abs(VaR.pick-VaR.real.mat))/VaR.real)

# CTE 估计 RMSE
rmse.CTE <- rbind(sqrt(colMeans((CTE.NLS-CTE.real.mat)^2)),sqrt(colMeans((CTE.WNLS-CTE.real.mat)^2)),
              sqrt(colMeans((CTE.pick-CTE.real.mat)^2)))

# CTE 估计 ARB
arb.CTE <- rbind(colMeans(abs(CTE.NLS-CTE.real.mat))/CTE.real,colMeans(abs(CTE.NLS-CTE.real.mat))/CTE.real,
              colMeans(abs(CTE.pick-CTE.real.mat))/CTE.real)

print(failure)