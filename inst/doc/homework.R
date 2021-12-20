## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## -----------------------------------------------------------------------------
library(knitr)#加载knitr package

#生成Rayleigh分布随机数
rrayleigh = function (size, sigma) {
  stopifnot(sigma > 0)#参数sigma>0
  u <-  runif(size)#均匀分布随机数
  sqrt(-2*sigma^2*log(u))#通过指数随机数生成Rayleigh随机数，其中-log(u)为指数随机数
}

#Rayleigh分布密度函数
drayleigh = function (x, sigma) {
  stopifnot(sigma > 0)#参数sigma>0
  y <-  x / (sigma^2) * exp(- x^2 / (2 * sigma^2))#Rayleigh分布密度函数
  y[x < 0] <-  0#x<0时密度为0
  y
}

#对不同的sigma=1, 2, 3, 4, 5分别生成Rayleigh随机数
size <- 1500#样本量
sigmas <- c(1, 2, 3, 4, 5)#参数sigma取值
samples <- sapply(sigmas, function(x) rrayleigh(size, x))#利用sapply对不同sigma生成Rayleigh随机数
samples <- data.frame(samples)#转为数据框
dimnames(samples)[[2]] <- paste("sigma =", sigmas)#重命名各列
kable(head(samples), caption = "表1.1：sigma=1, 2, 3, 4, 5时生成的Rayleigh随机数")#绘制数据样本表格

## ----fig.cap="图2.1：p1=0.75时复合分布的样本直方图和密度曲线"-----------------
#生成复合分布随机数
rmix = function(size, p){
  x1 <- rnorm(size, 0, 1)#生成x1
  x2 <- rnorm(size, 3, 1)#生成x2
  u <- runif(size)#生成均匀分布随机数
  t <- as.integer(u<p)#t以概率p取1，以概率1-p取0
  x <- t*x1 + (1-t)*x2#生成复合分布随机数x
}

n <- 1000#样本量1000
p1 <- 0.75#p1取0.75
x <- rmix(n, p1)#生成1000个复合分布随机数

#绘制直方图和密度曲线
hist(x, freq = F, ylim=c(0, 0.3), main="", xlab="", ylab="")
lines(density(x), col="red")

## ----fig.cap="图2.2：不同p1取值下复合分布的样本直方图和密度曲线"--------------
# par(mfrow = c(3,3))#以下生成3*3图形矩阵
# 
# for(p1 in seq(0.1, 0.9, 0.1)){
#   x <- rmix(n, p1)#对不同p1生成复合分布随机数
#   #绘制直方图和密度曲线
#   hist(x, freq = F, xlim=c(-6, 6), ylim=c(0, 0.4), main=paste("p1 = ", p1), xlab="", ylab="")
#   lines(density(x), col="red")
# }

## -----------------------------------------------------------------------------
#泊松过程参数lambda, t；Gamma分布参数shape, scale
#生成复合泊松-伽玛过程的算法
PoiGammaProc = function(lambda, t, shape, scale){
  Tn <- rexp(1000, lambda)#不知道需要多少Tn，因此生成充分多的Tn。Tn为泊松过程的间隔时间服从指数分布
  Sn <- cumsum(Tn)#计算泊松过程到达时间
  n <- min(which(Sn > t)) - 1#[0, t]内的到达次数即为泊松过程随机数
  Y <- rgamma(n, shape = shape, scale = 1/scale)#生成N(t)个iid伽玛随机数Gamma(shape, scale)
  X <- sum(Y)#对N(t)个iid伽玛随机数求和，得到复合泊松-伽玛过程随机数
}

## -----------------------------------------------------------------------------
t <- 10#取t=10
lambda <- c(1:3)#变化泊松过程强度参数lambda=1,2,3
shape <- 4#伽玛分布形状参数取4
scale <- 2#伽玛分布尺度参数取2
size <- 1000#样本量为1000
#对不同参数lambda分别生成复合泊松-伽玛过程
Xt1 <- replicate(size, PoiGammaProc(lambda[1], t, shape, scale))#lambda=1
Xt2 <- replicate(size, PoiGammaProc(lambda[2], t, shape, scale))#lambda=2
Xt3 <- replicate(size, PoiGammaProc(lambda[3], t, shape, scale))#lambda=3
mat <- cbind(Xt1, Xt2, Xt3)#将三组样本合并为样本矩阵

#均值的对比
mean.s <- apply(mat, 2, mean, trim=.5)#利用apply函数计算模拟样本均值
mean.th <- sapply(lambda, function(x){
  x*t*shape/scale
})#计算理论复合泊松-伽玛过程的均值
mean.mat <- rbind(mean.s, mean.th)
kable(mean.mat, caption = "表3.1：复合泊松-伽玛过程模拟样本均值和理论均值对比，变化泊松过程强度参数lambda=1, 2, 3")

#方差的对比
var.s <- (1-1/size)*apply(mat, 2, var)#利用apply函数计算模拟样本方差
var.th <- sapply(lambda, function(x){
  x*t*(shape^2+shape)/scale^2
})#计算理论复合泊松-伽玛过程的方差
var.mat <- rbind(var.s, var.th)
kable(var.mat, caption = "表3.2：复合泊松-伽玛过程模拟样本方差和理论方差对比，变化泊松过程强度参数lambda=1, 2, 3")

## -----------------------------------------------------------------------------
#利用MC方法近似计算Beta分布函数
cdfBeta <- function(x, alpha=3, beta=3, n=10000){ #x为需计算的点，alpha, beta为Beta分布参数，n为MC方法模拟次数
  stopifnot(x>0 && x<1)#排除定义域(0, 1)外的取值
  u <- runif(n, 0, x)#10000个[0,x]上的均匀分布随机数
  g <- (1/beta(alpha, beta)) * u^(alpha-1) * (1-u)^(beta-1)#计算g(x)
  cdf <- (x-0) * mean(g)#求均值得到cdf
  return(min(1, cdf))#返回cdf取值，由于cdf函数值不超过1，所以如果模拟结果大于1则取cdf为1
}

## -----------------------------------------------------------------------------
library(knitr)

set.seed(504)
x <- seq(0.1, 0.9, 0.1)#生成x=0.1, ..., 0.9
cdf_MC <- sapply(x, function(z) cdfBeta(z))#利用sapply对每个x利用MC方法计算对应的cdf
cdf_pbeta <- pbeta(x, 3, 3)#利用R内置函数pbeta计算cdf

result <- cbind(x, cdf_MC, cdf_pbeta)#合并计算结果
dimnames(result)[[2]] <- c("x", "Monte Carlo", "pbeta")#重命名各列
kable(result, caption = "表1.1：Monte Carlo方法与R内置函数pbeta分别计算Beta(3, 3)分布函数的结果对比")#利用kable绘制计算结果表格

## -----------------------------------------------------------------------------
#对偶变量法生成Rayleigh分布随机数
rrayleigh = function(sigma, n=10000, antithetic=TRUE){#sigma为Rayleigh分布参数，n模拟次数，antithetic是否应用对偶变量，默认为是
  u <- runif(n)#生成n/2个均匀随机数
  if(!antithetic) v <- runif(n) else
    v <- 1-u#若应用对偶变量，则令后一半均匀随机数v为1-u，否则v为n个均匀随机数且与u独立
  X1 <- sqrt(-2*sigma^2*log(u))#利用u生成一组Rayleigh分布随机数X1
  X2 <- sqrt(-2*sigma^2*log(v))#利用v生成一组Rayleigh分布随机数X2
  X <- (X1+X2)/2#取两组样本均值得到最终估计
}

## -----------------------------------------------------------------------------
sigmas <- c(1:5)#初始化参数sigma
var_reduction <- matrix(nrow = 1, ncol = length(sigmas))#存储结果矩阵
set.seed(509)
for(i in 1:length(sigmas)){
  simple <- rrayleigh(sigmas[i], antithetic = FALSE)#利用简单方法生成Rayleigh分布随机数，X1与X2独立
  antithetic <- rrayleigh(sigmas[i])#利用对偶变量生成Rayleigh分布随机数，X1与X2负相关
  var_reduction[1, i] <- (var(simple)-var(antithetic))/var(simple)#对偶变量法相比简单方法的方差减少比例
}
dimnames(var_reduction)[[2]] <- paste("sigma =", sigmas)#重命名各列
kable(var_reduction, caption = "表2.1：对偶变量法生成Rayleigh随机数相比简单方法生成Rayleigh随机数方差减少的比例")#绘制结果表格

## ----fig.cap="图3.1：g/f的函数图像", fig.height=7-----------------------------
g = function (x) {
  x ^ 2 / sqrt(2*pi) * exp(-x^2/2) * (x>1)
}

x <- seq(1, 10, 0.1)#取点
y_g <- g(x)#g(x)函数值
y_f1 <- 1 / x^2#f1函数值
y_f2 <- exp(1-x)#f2函数值

#绘制g/f的函数图像
ylim <- max(c(y_g/y_f1, y_g/y_f2))
gs <- c(expression(f[1](x)==1/x^2),
        expression(f[2](x)==e^{1-x}))
w <- 2
plot(x, y_g/y_f1, type = "l", ylab = "",
        ylim = c(0, ylim), lwd = w, lty = 2, col=2)
lines(x, y_g/y_f2, lty = 3, lwd = w,col=3)
legend("topright", legend = gs,
       lty = 2:3, lwd = w, inset = 0.02,col=2:3)

## -----------------------------------------------------------------------------
n <- 10000
theta.hat <- se <- numeric(2)#存储估计值theta.hat和标准差se
set.seed(513)

#importance function取f1，逆变换法
u <- runif(n)
x <- 1 / (1-u)
fg <- g(x) / (1/x^2)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)

#importance function取f2，逆变换法
u <- runif(n)
x <- 1-log(u)
fg <- g(x) / exp(1-x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)

result <- cbind(theta.hat, se)
dimnames(result)[[1]] <- c("f1", "f2")
kable(result, caption = "表3.1：利用f1和f2的importance sampling所生成的估计和标准差对比")#kable绘制结果表格

## -----------------------------------------------------------------------------
n <- 10000
theta.hat <- se <- numeric(1)#存储估计值theta.hat和标准差se
set.seed(514)

#importance function取f，逆变换法
u <- runif(n)
x <- 1 - log(u)
fg <- g(x) / exp(1-x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)

result <- cbind(theta.hat, se)
kable(result, caption = "表4.1：利用importance sampling所生成的MC估计和标准差")#kable绘制结果表格

## -----------------------------------------------------------------------------
library(knitr)
set.seed(605)
calcCI <- function(n, alpha){
  x <- rchisq(n, 2)
  c1 <- mean(x) - sd(x)*qt(1-alpha/2, n-1)/sqrt(n)
  c2 <- mean(x) + sd(x)*qt(1-alpha/2, n-1)/sqrt(n)
  y <- as.integer(2>=c1 && 2<=c2)
  return(y)
}

n <- 20
alpha <- 0.05
m <- 5000
ys <- replicate(m, expr = calcCI(n, alpha))
kable(mean(ys), caption = "表1.1：样本来自chi^2(2)时t置信区间的覆盖概率")

## -----------------------------------------------------------------------------
set.seed(605)
n <- 20
alpha <- 0.05
m <- 5000
ys <- replicate(m, expr = {
  x <- rchisq(n, 2)
  c <- (n-1)*var(x) / qchisq(alpha, df=n-1)
  y <- as.integer(c>4)
  return(y)
})
kable(mean(ys), caption = "表1.2：样本来自chi^2(2)时方差置信区间的覆盖概率")

## -----------------------------------------------------------------------------
set.seed(0)
n <- 20
alpha <- 0.05
mu0 <- 1

m <- 10000
ps <- replicate(m, expr = {
  x <- rchisq(n, mu0)
  ttest <- t.test(x, alternative = "two.sided", mu=mu0)
  p <- ttest$p.value
  return(p)
})

p.hat <- mean(ps < alpha)
kable(p.hat, caption = "表2.1：样本来自chi^2(1)时均值t检验的第I类错误概率")

## -----------------------------------------------------------------------------
set.seed(0)
n <- 20
alpha <- 0.05
mu0 <- 1

m <- 10000
ps <- replicate(m, expr = {
  x <- runif(n, 0, 2)
  ttest <- t.test(x, alternative = "two.sided", mu=mu0)
  p <- ttest$p.value
  return(p)
})

p.hat <- mean(ps < alpha)
kable(p.hat, caption = "表2.2：样本来自U[0, 2]时均值t检验的第I类错误概率")

## -----------------------------------------------------------------------------
set.seed(0)
n <- 20
alpha <- 0.05
mu0 <- 1

m <- 10000
ps <- replicate(m, expr = {
  x <- rexp(n, 1)
  ttest <- t.test(x, alternative = "two.sided", mu=mu0)
  p <- ttest$p.value
  return(p)
})

p.hat <- mean(ps < alpha)
kable(p.hat, caption = "表2.3：样本来自Exp(1)时均值t检验的第I类错误概率")

## -----------------------------------------------------------------------------
mat <- matrix(nrow = 3, ncol = 3)
mat[1,1:3] <- c("a","b","a + b")
mat[2,1:3] <- c("c","d","c + d")
mat[3,1:3] <- c("a + c","b + d","10000")
colnames(mat) <- c("Test2 reject","Test2 accpet","Row total")
rownames(mat) <- c("Test1 reject","Test1 accpet","Column total")
kable(mat)

## -----------------------------------------------------------------------------
library(knitr)
library(MASS)

#计算多元偏度统计量的函数
mvsk <- function(x){
  xx <- scale(x, center=T, scale=F)#将样本矩阵x中心化
  n <- nrow(x)#获取样本量
  sk <- 1/n^2 * sum((xx %*% solve(cov(x)) %*% t(xx))^3)#检验统计量
  return(sk)
}

ns <- c(10, 20, 30, 50, 100, 500)#样本量
d <- 2#维数
alpha <- 0.05#显著性水平
cv <- qchisq(1-alpha, d*(d+1)*(d+2)/6)#检验临界值

p.reject <- numeric(length(ns))#存储第I类错误估计值
m <- 1e4#模拟次数
set.seed(99)
# for(i in 1:length(ns)){
#   mvsktests <- numeric(m)#存储每次试验的检验结果
#   for(j in 1:m){
#     x <- mvrnorm(ns[i], mu = rep(0, d), Sigma = diag(d))#生成n个多元正态分布样本，维数为2
#     mvsktests[j] <- as.integer(ns[i]*mvsk(x)/6 >= cv)#拒绝H0取1，否则取0
#   }
#   p.reject[i] <- mean(mvsktests)#拒绝H0的比例
# }
# 
# #绘制结果表格
# result <- matrix(c(ns, p.reject), nrow = 2, ncol = length(ns), byrow = T)
# dimnames(result)[[1]] <- c("n", "Type I error rate")
# kable(result, caption = "表1：Mardia多元偏度检验第I类错误概率估计")


## ---- fig.cap="图1：在受污染的多元正态总体下，Mardia多元偏度检验势的估计"-----
alpha <- 0.1#显著性水平
n <- 500#样本量
m <- 2500#模拟次数
d <- 2#维数
epsilon <- c(seq(0, 0.15, 0.01), seq(0.15, 1, 0.05))#混合分布比例序列e
N <- length(epsilon)
pwr <- numeric(N)#存储每次模拟检验的势
cv <- qchisq(1-alpha, d*(d+1)*(d+2)/6)#检验临界值

set.seed(100)
# for(i in 1:N){ #对每个epsilon
#   mvsktests <- numeric(m)#存储每次试验结果，拒绝H0取1，否则取0
#   e <- epsilon[i]#混合分布比例e
#   for(j in 1:m){ #对每次试验
#     x1 <- mvrnorm(n, rep(0, d), diag(d))#生成多元正态样本N(0, I)
#     x2 <- mvrnorm(n, rep(0, d), 100*diag(d))#生成多元正态样本N(0, 100*I)
#     emat <- diag(as.integer(runif(n) <= e))#为生成混合分布构造系数矩阵，对角线元素以概率e取1，以概率1-e取0
#     x <- (diag(n)-emat) %*% x1 + emat %*% x2#生成混合分布样本，每行以概率1-e取N(0, I)，以概率e取N(0, 100*I)
#     mvsktests[j] <- as.integer(n*mvsk(x)/6 >= cv)#检验结果拒绝H0取1，否则取0
#   }
#   pwr[i] <- mean(mvsktests)#估计检验的势
# }
# 
# #绘制power vs epsilon曲线
# plot(epsilon, pwr, type = "b",
#      xlab = bquote(epsilon), ylim = c(0,1))
# abline(h = alpha, lty = 3)
# se <- sqrt(pwr * (1-pwr) / m) #添加标准差
# lines(epsilon, pwr+se, lty = 3)
# lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(knitr)
library(bootstrap)
library(boot)

Sigma <- cov(scor)#计算协方差矩阵的MLE
lambda_hat <- eigen(Sigma)$values#特征值
theta_hat <- lambda_hat[1]/sum(lambda_hat)#theta的样本估计

#计算theta bootstrap估计的函数，以传给boot的参数statistic
theta_i <- function(dat, inds){
  Sigma <- cov(dat[inds,])
  lambda <- eigen(Sigma)$values
  return(lambda[1]/sum(lambda))
}

set.seed(77)
obj_boot <- boot(data=scor, statistic=theta_i, R=200)#利用boot进行bootstrap抽样
theta_boot <- as.numeric(obj_boot$t)#bootstrap估计值
bias_boot <- mean(theta_boot) - theta_hat#bootstrap的偏差估计
se_boot <- sd(theta_boot)#bootstrap的标准差估计

#绘制结果表格
result <- matrix(c(theta_hat, bias_boot, se_boot), nrow = 1)
dimnames(result)[[2]] <- c("theta_hat", "bias_boot", "se_boot")
kable(result, caption = "表1.1：theta的样本估计和bootstrap的偏差和标准差估计结果")

## -----------------------------------------------------------------------------
Sigma <- cov(scor)#计算协方差矩阵的MLE
lambda_hat <- eigen(Sigma)$values#特征值
theta_hat <- lambda_hat[1]/sum(lambda_hat)#theta的样本估计

n <- nrow(scor)
theta_jack <- numeric(n)
for(i in 1:n){
  dat <- scor[-i,]
  Sigma <- cov(dat)
  lambda <- eigen(Sigma)$values
  theta_jack[i] <- lambda[1]/sum(lambda)
}
bias_jack <- (n-1)*(mean(theta_jack) - theta_hat)
se_jack <- (n-1)*sqrt(var(theta_jack)/n)

#绘制结果表格
result <- matrix(c(bias_jack, se_jack), nrow = 1)
dimnames(result)[[2]] <- c("bias_jack", "se_jack")
kable(result, caption = "表2.1：jackknife的偏差和标准差估计结果")

## -----------------------------------------------------------------------------
set.seed(79)
obj_boot <- boot(data=scor, statistic=theta_i, R=2000)#利用boot函数进行bootstrap抽样
boot.ci(obj_boot, type = c("perc", "bca"))#利用boot.ci函数计算95% percentile和BCa置信区间

## -----------------------------------------------------------------------------
#计算样本偏度
sk_i <- function(dat, inds){
  x <- dat[inds]
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  sk <- m3 / m2^(3/2)
  return(sk)
}

alpha <- 0.05#显著性水平
n <- 200#样本量
m <- 5e3#MC模拟次数
sk <- 0#正态分布偏度系数
y <- matrix(nrow = m, ncol = 3)#存储CI是否覆盖参数
CI_miss <- matrix(nrow = m, ncol = 6)#存储CI位置，1-3列为参数是否落在3种CI左边，4-6列为参数是否落在3种CI右边

set.seed(0)
# for(i in 1:m){
#   x <- rnorm(n)#生成正态样本
#   obj_boot <- boot(data = x, statistic = sk_i, R=2000)#利用boot函数进行bootstrap抽样
#   obj_ci <- boot.ci(obj_boot, conf = 1-alpha, type = c("norm","basic", "perc"))#利用boot.ci函数计算95%置信区间
#   LCLs <- c(obj_ci[["normal"]][1,2], obj_ci[["basic"]][1,4], obj_ci[["percent"]][1,4])#置信区间下界
#   UCLs <- c(obj_ci[["normal"]][1,3], obj_ci[["basic"]][1,5], obj_ci[["percent"]][1,5])#置信区间上界
#   y[i,] <- as.integer(sk > LCLs & sk < UCLs)#CI是否覆盖参数
#   CI_miss[i,] <- as.integer(c(sk < LCLs, sk > UCLs))#CI位置，1-3列为参数是否落在3种CI左边，4-6列为参数是否落在3种CI右边
# }
# 
# #绘制结果表格
# cp <- apply(y, 2, mean)
# cp <- matrix(cp, nrow = 1)
# dimnames(cp)[[2]] <- c("standard normal bootstrap", "basic bootstrap", "percentile")
# kable(cp, caption = "表4.1：三种bootstrap置信区间覆盖概率的MC估计（正态样本）")
# 
# CI_left <- apply(CI_miss[,1:3], 2, mean)
# CI_right <- apply(CI_miss[,4:6], 2, mean)
# miss_prop <- rbind(CI_left, CI_right)
# dimnames(miss_prop)[[2]] <- c("standard normal bootstrap", "basic bootstrap", "percentile")
# kable(miss_prop, caption = "表4.2：偏度位于三种bootstrap置信区间左边/右边的比例（正态样本）")


## -----------------------------------------------------------------------------
n <- 200#样本量
m <- 5e3#MC模拟次数
sk <- sqrt(8/5)#卡方分布chi^2(5)偏度系数
y <- matrix(nrow = m, ncol = 3)#存储CI是否覆盖参数
CI_miss <- matrix(nrow = m, ncol = 6)#存储CI位置，1-3列为参数是否落在3种CI左边，4-6列为参数是否落在3种CI右边

set.seed(0)
# for(i in 1:m){
#   x <- rchisq(n, 5)#生成卡方样本
#   obj_boot <- boot(data = x, statistic = sk_i, R=2000)#利用boot函数进行bootstrap抽样
#   obj_ci <- boot.ci(obj_boot, conf = 1-alpha, type = c("norm","basic", "perc"))#利用boot.ci函数计算95%置信区间
#   LCLs <- c(obj_ci[["normal"]][1,2], obj_ci[["basic"]][1,4], obj_ci[["percent"]][1,4])#置信区间下界
#   UCLs <- c(obj_ci[["normal"]][1,3], obj_ci[["basic"]][1,5], obj_ci[["percent"]][1,5])#置信区间上界
#   y[i,] <- as.integer(sk > LCLs & sk < UCLs)#CI是否覆盖参数
#   CI_miss[i,] <- as.integer(c(sk < LCLs, sk > UCLs))#CI位置，1-3列为参数是否落在3种CI左边，4-6列为参数是否落在3种CI右边
# }
# 
# #绘制结果表格
# cp <- apply(y, 2, mean)
# cp <- matrix(cp, nrow = 1)
# dimnames(cp)[[2]] <- c("standard normal bootstrap", "basic bootstrap", "percentile")
# kable(cp, caption = "表4.3：三种bootstrap置信区间覆盖概率的MC估计（卡方样本）")
# 
# CI_left <- apply(CI_miss[,1:3], 2, mean)
# CI_right <- apply(CI_miss[,4:6], 2, mean)
# miss_prop <- rbind(CI_left, CI_right)
# dimnames(miss_prop)[[2]] <- c("standard normal bootstrap", "basic bootstrap", "percentile")
# kable(miss_prop, caption = "表4.4：偏度位于三种bootstrap置信区间左边/右边的比例（卡方样本）")


## -----------------------------------------------------------------------------
library(knitr)
library(boot)

data("iris")
setosa <- as.matrix(iris[1:50, -5])
x <- setosa[,1:2]#sepal
y <- setosa[,3:4]#petal

#计算rho permutation估计的函数，以传给boot的参数statistic
rho_i <- function(dat, inds){
  x <- dat[,1:2]
  y <- dat[inds, 3:4]
  rho <- cor.test(x, y, method = "spearman")$estimate
  return(rho)
}

set.seed(802)
obj_boot <- boot(data=setosa, statistic=rho_i, sim="permutation", R=999)#permutation
rho_perm <- c(obj_boot$t, obj_boot$t0)
phat <- mean(rho_perm > obj_boot$t0)#显著性水平的permutation估计
phat


## -----------------------------------------------------------------------------
cor.test(x, y, method = "spearman")$p.value

## -----------------------------------------------------------------------------
# library(RANN)
# library(energy)
# library(Ball)

# #NN统计量
# Tn <- function(dat, inds, sizes, k){
#   n1 <- sizes[1]
#   n2 <- sizes[2]
#   n <- n1 + n2
#   if(is.vector(dat)) 
#     dat <- data.frame(dat)
#   dat <- dat[inds, ]
#   NN <- nn2(data=dat, k=k+1)
#   block1 <- NN$nn.idx[1:n1,-1]
#   block2 <- NN$nn.idx[(n1+1):n,-1]
#   i1 <- sum(block1 <= n1)
#   i2 <- sum(block2 > n1)
#   return((i1 + i2) / (k * n))
# }
# 
# #NN假设检验
# eqdist.nn <- function(dat, sizes, k){
#   obj_boot <- boot(data=dat, statistic=Tn, R=R, sim = "permutation", sizes = sizes, k=k)
#   ts <- c(obj_boot$t0, obj_boot$t)
#   p.value <- mean(ts>=ts[1])
#   return(list(statistic=ts[1], p.value=p.value))
# }
# 
# m <- 1000#模拟次数
# d <- 2#维数
# k <- 3#NN检验参数取3
# n1 <- n2 <- 50#样本量
# R <- 999#permutation个数
# n <- n1 + n2
# N <- c(n1, n2)
# p.values <- matrix(nrow = m, ncol = 3)#存储每次模拟p值
# alpha <- 0.05
# 
# #Unequal variances and equal expectations
# sigma <- 1.5
# set.seed(11041)
# for(i in 1:m){
#   x <- matrix(rnorm(n1*d, 0, 1), ncol = d)
#   y <- matrix(rnorm(n2*d, 0, sigma), ncol = d)
#   dat <- rbind(x, y)
#   p.values[i, 1] <- eqdist.nn(dat, N, k)$p.value#NN
#   p.values[i, 2] <- eqdist.etest(dat, sizes = N, R=R)$p.value#energy
#   p.values[i, 3] <- bd.test(x=x, y=y, num.permutations = R, seed = i*11041)$p.value#ball
# }
# pow <- matrix(colMeans(p.values < alpha), ncol = 3)#模拟得到三种检验的势
# dimnames(pow)[[2]] <- c("NN", "energy", "ball")
# kable(pow, caption = "表2.1：NN, energy和ball三种检验方法的势(Unequal variances and equal expectations)")

## -----------------------------------------------------------------------------
mu <- 0.5
sigma <- 1.5
p.values <- matrix(nrow = m, ncol = 3)#存储每次模拟p值
set.seed(11042)
# for(i in 1:m){
#   x <- matrix(rnorm(n1*d, 0, 1), ncol = d)
#   y <- matrix(rnorm(n2*d, mu, sigma), ncol = d)#Unequal variances and unequal expectations
#   dat <- rbind(x, y)
#   p.values[i, 1] <- eqdist.nn(dat, N, k)$p.value#NN
#   p.values[i, 2] <- eqdist.etest(dat, sizes = N, R=R)$p.value#energy
#   p.values[i, 3] <- bd.test(x=x, y=y, num.permutations = R, seed = i*11041)$p.value#ball
# }
# pow <- matrix(colMeans(p.values < alpha), ncol = 3)#模拟得到三种检验的势
# dimnames(pow)[[2]] <- c("NN", "energy", "ball")
# kable(pow, caption = "表2.2：NN, energy和ball三种检验方法的势(Unequal variances and unequal expectations)")

## -----------------------------------------------------------------------------
#生成bimodel distribution随机数
rbinorm <- function(n, mu=2, sigma=2){
  x <- numeric(n)
  u <- runif(n)
  inds <- which(u < 0.5)
  x[inds] <- rnorm(sum(u < 0.5))
  x[-inds] <- rnorm(n - sum(u < 0.5), mu, sigma)
  return(x)
}

p.values <- matrix(nrow = m, ncol = 3)#存储每次模拟p值
set.seed(11043)
# for(i in 1:m){
#   x <- rt(n1, df=1)#t distribution
#   y <- rbinorm(n2)#bimodel distribution
#   dat <- c(x, y)
#   p.values[i, 1] <- eqdist.nn(dat, N, k)$p.value#NN
#   p.values[i, 2] <- eqdist.etest(dat, sizes = N, R=R)$p.value#energy
#   p.values[i, 3] <- bd.test(x=x, y=y, num.permutations = R, seed = i*11041)$p.value#ball
# }
# pow <- matrix(colMeans(p.values < alpha), ncol = 3)#模拟得到三种检验的势
# dimnames(pow)[[2]] <- c("NN", "energy", "ball")
# kable(pow, caption = "表2.3：NN, energy和ball三种检验方法的势(Non-normal distributions)")
# ```
# 从表2.3可以看到，对于非正态分布情形，energy相比NN和ball检验表现更好。
# 
# (4) Unbalanced samples (say, 1 case versus 10 controls)
# ```{r}
# n1 <- 10
# n2 <- 10 * n1#Unbalanced samples
# n <- n1 + n2
# N <- c(n1, n2)
# mu <- 0.5
# p.values <- matrix(nrow = m, ncol = 3)#存储每次模拟p值
# set.seed(11044)
# for(i in 1:m){
#   x <- rnorm(n1, 0, 1)
#   y <- rnorm(n2, mu, 1)
#   dat <- c(x, y)
#   p.values[i, 1] <- eqdist.nn(dat, N, k)$p.value#NN
#   p.values[i, 2] <- eqdist.etest(dat, sizes = N, R=R)$p.value#energy
#   p.values[i, 3] <- bd.test(x=x, y=y, num.permutations = R, seed = i*11041)$p.value#ball
# }
# pow <- matrix(colMeans(p.values < alpha), ncol = 3)#模拟得到三种检验的势
# dimnames(pow)[[2]] <- c("NN", "energy", "ball")
# kable(pow, caption = "表2.4：NN, energy和ball三种检验方法的势(Unbalanced samples)")

## -----------------------------------------------------------------------------
#柯西分布密度
f <- function(x){
  fx <- 1 / (pi*(1+x^2))
  return(fx)
}

set.seed(1193)
m <- 1e4
x <- numeric(m)
x[1] <- rnorm(1)
k <- 0#记录拒绝次数
u <- runif(m)

#生成马氏链
for(i in 2:m){
  xt <- x[i-1]
  y <- rnorm(1, xt, 1)
  num <- f(y) * dnorm(xt, y, 1)
  den <- f(xt) * dnorm(y, xt, 1)
  r <- num / den
  if(u[i] <= r){
    x[i] <- y
  }else{
    x[i] <- xt
    k <- k+1
  }
}

#比较MH和真实柯西分布的分位数
index <- c(1001:m)
yMH <- x[index]
ps <- seq(0.1, 0.9, 0.1)
result <- rbind(quantile(yMH, ps), qcauchy(ps))#绘制结果表格
dimnames(result)[[1]] <- c("Metropolis-Hastings", "True")
kable(result, caption = "表1.1：分位数对比结果")


## ---- fig.cap="图2.1：Gibbs sampler生成的二元随机数"--------------------------
#Gibbs sampler function, 固定n=10, a=1, b=2
Gibbs <- function(N, x0, y0, n=10, a=1, b=2){
  X <- matrix(nrow = N, ncol = 2)
  
  #初始化
  x <- x0
  y <- y0
  X[1,] <- c(x, y)
  
  for(i in 2:N){
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, y)#x|y为Binomial(n,y)
    x <- X[i, 1]
    X[i, 2] <- rbeta(1, x+a, n-x+b)#y|x为Beta(x+a, n-x+b)
  }
  return(X)
}

set.seed(1198)
n <- 10
a <- 1
b <- 2
N <- 5000
x0 <- 0.5*n
y0 <- (x0+a) / (n+a+b)#初始化
X <- Gibbs(N, x0, y0)
burn <- 1000
XGibbs <- X[(burn+1):N,]
#plot(XGibbs, cex=0.5, xlab="x", ylab="y")#XY散点图

## -----------------------------------------------------------------------------
#生成柯西分布马氏链
cauchy.chain <- function(x0, m){
  x <- numeric(m)
  x[1] <- x0
  k <- 0#记录拒绝次数
  u <- runif(m)
  
  #生成马氏链
  for(i in 2:m){
    xt <- x[i-1]
    y <- rnorm(1, xt, 1)
    num <- dcauchy(y) * dnorm(xt, y, 1)
    den <- dcauchy(xt) * dnorm(y, xt, 1)
    r <- num / den
    if(u[i] <= r){
      x[i] <- y
    }else{
      x[i] <- xt
      k <- k+1
    }
  }
  
  return(x)
}

#Gelman-Rubin方法计算Rhat
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i,1:j]) for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)

  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

set.seed(113)
k <- 4#chain的个数
m <- 15000#chain的长度
burn <- 1000
x0 <- c(-10, -5, 5, 10)#初始值

# X <- matrix(nrow = k, ncol = m)
# for(i in 1:k)
#   X[i,] <- cauchy.chain(x0[i], m)
# 
# #计算检验统计量psi
# psi <- t(apply(X, 1, cumsum))
# for(i in 1:nrow(psi))
#   psi[i,] <- psi[i,] / (1:ncol(psi))
# 
# #绘制psi
# par(mfrow=c(2, 2))
# for(i in 1:k)
#   plot(psi[i, (burn+1):m], type="l", xlab=i, ylab=bquote(psi))
# 
# #绘制Rhat
# par(mfrow=c(1, 1))
# rhat <- numeric(m)
# for(j in (burn+1):m)
#   rhat[j] <- Gelman.Rubin(psi[,1:j])
# plot(rhat[(burn+1):m], type="l", xlab="", ylab="R")
# abline(h=1.2, lty=2)


## -----------------------------------------------------------------------------
set.seed(98)
k <- 3#chain的个数
m <- 15000#chain的长度
burn <- 1000
n <- 10
a <- 1
b <- 2
x0 <- c(0, 5, 10)#初始值
y0 <- (x0+a) / (n+a+b)#初始化

X <- matrix(nrow = k, ncol = m)
Y <- matrix(nrow = k, ncol = m)
# for(i in 1:k){
#   XY <- t(Gibbs(m, x0[i], y0[i]))
#   X[i,] <- XY[1,]
#   Y[i,] <- XY[2,]
# }
  
#对X的边缘分布利用Gelman-Rubin方法观测收敛性
#计算检验统计量psi
psi <- t(apply(X, 1, cumsum))
for(i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

#绘制psi
# par(mfrow=c(2, 2))
# for(i in 1:k)
#   plot(psi[i, (burn+1):m], type="l", xlab=i, ylab=bquote(psi))

## -----------------------------------------------------------------------------
# #绘制Rhat
# par(mfrow=c(1, 1))
# rhat <- numeric(m)
# for(j in (burn+1):m)
#   rhat[j] <- Gelman.Rubin(psi[,1:j])
# plot(rhat[(burn+1):m], type="l", xlab="", ylab="R")
# abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
# #对Y的边缘分布利用Gelman-Rubin方法观测收敛性
# #计算检验统计量psi
# psi <- t(apply(Y, 1, cumsum))
# for(i in 1:nrow(psi))
#   psi[i,] <- psi[i,] / (1:ncol(psi))
# 
# #绘制psi
# par(mfrow=c(2, 2))
# for(i in 1:k)
#   plot(psi[i, (burn+1):m], type="l", xlab=i, ylab=bquote(psi))

## -----------------------------------------------------------------------------
# #绘制Rhat
# par(mfrow=c(1, 1))
# rhat <- numeric(m)
# for(j in (burn+1):m)
#   rhat[j] <- Gelman.Rubin(psi[,1:j])
# plot(rhat[(burn+1):m], type="l", xlab="", ylab="R")
# abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
Termk <- function(a, k){
  d <- length(a)
  return(
    (-1)^k * gamma((d+1)/2) / ((2*k+1)*(2*k+2)) * exp(-lgamma(k+1) - k*log(2) + (2*k+2)*log(norm(a, type = "2")) + lgamma(k+3/2) - lgamma(k+d/2+1))
  )
}

## -----------------------------------------------------------------------------
SumTerm <- function(a){
  K <- c(0:200)#求和至k=200
  return(sum(Termk(a, K)))
}

## -----------------------------------------------------------------------------
a <- c(1, 2)
SumTerm(a)

## -----------------------------------------------------------------------------
ck <- function(a, k){
  sqrt((a^2*k) / (k+1-a^2))
}#积分限

g <- function(u, k){
  (1 + u^2/k)^(-(k+1)/2)
}#被积函数

intg <- function(ck, k){
  integrate(g, lower = 0, upper = ck, 
            k = k)
}#积分项

bk <- function(k){
  2/sqrt(pi*k)*exp(lgamma((k+1)/2) - lgamma(k/2))
}#常数项

f <- function(a, k){
  bk(k)*intg(ck(a, k), k)$value - bk(k-1)*intg(ck(a, k-1), k-1)$value
}#转化为求f(a ,k)=0

#求根
eps <- .Machine$double.eps^0.25
root_115 <- function(k){
  uniroot(f, interval = c(eps, 2), k=k)$root
}


## -----------------------------------------------------------------------------
S <- function(a,k){
  pt(q=sqrt(a^2*k / (k+1-a^2)), df=k, lower.tail = F)
}

obj <- function(a, k){
  S(a, k) - S(a, k-1)
}

root_114 <- function(k){
  uniroot(obj, interval = c(eps, 2), k=k)$root
}

## -----------------------------------------------------------------------------
ks <- c(4:25, 100, 500, 1000)
result <- matrix(nrow = length(ks), ncol = 3)
result[,1] <- ks
for(i in 1:length(ks)){
  result[i, 2:3] <- c(root_115(ks[i]), root_114(ks[i]))
}
dimnames(result)[[2]] <- c("k", "ex11.5", "ex11.4")
kable(result, caption = "表2.1：ex11.5与ex11.4求根结果对比")

## -----------------------------------------------------------------------------
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
lambda0 <- 0.5

EMlambda <- function(y, lambda){
  N <- 10000
  tol <- .Machine$double.eps^0.5
  n <- length(y)
  nt <- sum(y==1)
  
  for(i in 1:N){
    total <- sum(y)
    lambda_old <- lambda
    lambda <- (total + nt*lambda_old) / n
    
    if(abs(lambda - lambda_old)/lambda_old <tol) break
  }
  
  return(lambda)
}

EMlambda(y, lambda0)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

fit3 <- lapply(formulas, lm, data = mtcars)
sapply(fit3, rsq)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

fit4 <- lapply(bootstraps, lm, formula = mpg ~ disp)
sapply(fit4, rsq)

## -----------------------------------------------------------------------------
dt <- as.data.frame(matrix(1:9, nrow = 3))
vapply(dt, sd, numeric(1))

## -----------------------------------------------------------------------------
#以iris数据集为例
vapply(iris[vapply(iris, is.numeric, logical(1))], sd, numeric(1))

## -----------------------------------------------------------------------------
# library(parallel)
# cores <-  detectCores()
mcsapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE){
  FUN <- match.fun(FUN)
  answer <- mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
      names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer))
      simplify2array(answer, higher = (simplify == "array"))
  else answer
}

## -----------------------------------------------------------------------------
#Gibbs sampler function, 固定n=10, a=1, b=2
#R function for ex9.8
gibbsR <- function(N, x0, y0, n=10, a=1, b=2){
  X <- matrix(nrow = N, ncol = 2)
  
  #初始化
  x <- x0
  y <- y0
  X[1,] <- c(x, y)
  
  for(i in 2:N){
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, y)#x|y为Binomial(n,y)
    x <- X[i, 1]
    X[i, 2] <- rbeta(1, x+a, n-x+b)#y|x为Beta(x+a, n-x+b)
  }
  return(X)
}

## -----------------------------------------------------------------------------
#Cpp function for ex9.8
library(Rcpp)
cppFunction('NumericMatrix gibbsC(int N, double x0, double y0, int n=10, int a=1, int b=2){
  NumericMatrix mat(N, 2);
  double x = x0, y = y0;
  mat(0, 0) = x;
  mat(0, 1) = y;

  for(int i = 1; i < N; i++){
    y = mat(i-1, 1);
    mat(i, 0) = rbinom(1, n, y)[0];
    x = mat(i, 0);
    mat(i, 1) = rbeta(1, x+a, n-x+b)[0];
  }
  return mat;
}')


## -----------------------------------------------------------------------------
set.seed(1198)
#初始化
n <- 10
a <- 1
b <- 2
N <- 5000
x0 <- 0.5*n
y0 <- (x0+a) / (n+a+b)

Chain_R <- gibbsR(N, x0, y0)
Chain_Cpp <- gibbsC(N, x0, y0)

##qqplot compare for x
par(mfrow=c(1,2))
burn <- 1000
x_R <- as.numeric(Chain_R[(burn+1):N,1])
x_Cpp <- as.numeric(Chain_Cpp[(burn+1):N,1])
p <- ppoints(100)
qqplot(quantile(x_R, p), quantile(x_Cpp, p), xlab = "R", ylab = "Rcpp")
##qqplot compare for y
y_R <- as.numeric(Chain_R[(burn+1):N,2])
y_Cpp <- as.numeric(Chain_Cpp[(burn+1):N,2])
qqplot(quantile(y_R, p), quantile(y_Cpp, p), xlab = "R", ylab = "Rcpp")

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(gR=gibbsR(N, x0, y0), gC=gibbsC(N, x0, y0))
summary(ts)[,c(1,3,5,6)]

