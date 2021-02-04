# Housingkeeping
options(repr.plot.width=14,repr.plot.antialias='subpixel',repr.plot.res=218)
sessionInfo()
update.packages()

# Browser spending analysis
# 1. Read data 
setwd("/Users/nicoleyin88/Documents/1. Big Data/3. Presentation/Sample Code/examples/")
browser <- read.csv("web-browsers.csv")

dim(browser)
head(browser)

# 2. histogram 
par(mai=c(.8,.8,.1,.1))
hist(log(browser$spend), freq=FALSE,         
     xaxt="n", main="", xlab="total online spend", col=8, border="grey90")
lgrid = c(1,10,100,1000,10000,100000)
axis(1, at=log(lgrid), labels=sprintf("%.0e",lgrid))

# 3. Basic stat 
nrow(browser)
mean(browser$spend)
var(browser$spend)/nrow(browser)

# 4. density graph 
xbar <- mean(browser$spend)
xbse <-  sd(browser$spend)/sqrt(nrow(browser))
xx <- seq(1650,2250,length=1000)
par(mai=c(.9,.8,.2,.2))
plot(xx, dnorm(xx, xbar, xbse), type="l", col="royalblue", lwd=1.5,
     xlab="average total online spend", ylab="density")

# 5. Nonparametric bootstrap for x "spend" (use sub-sample) 
set.seed(1234)
B <- 10000
mub <- c()
for (b in 1:B){
  samp_b = sample.int(nrow(browser), replace=TRUE)
  mub <- c(mub, mean(browser$spend[samp_b]))
}

sd(mub)

# 6. Bootstrap sampling distribution 
par(mai=c(.8,.8,.2,.2))
hist(mub, main="", xlab="average total online spend",
     col=8, border="grey90", freq=FALSE)
lines(xx, dnorm(xx, xbar, xbse), col="royalblue", lwd=1.5)

# 7. Parametric Bootstrap for x "spend" (use full sample)
xbar <- mean(browser$spend)
sig2 <-  var(browser$spend)
B <- 10000
mus <- c()
for(b in 1:B){
  xsamp <- rnorm(1e4, xbar, sqrt(sig2))
  mus <- c(mus, mean(xsamp))
}

sd(mus)      # sample 
sqrt(sig2/1e4)    # population 

# 8. Nonparametric bootstrap for CI (use sub-sample)
smallsamp <- browser$spend[sample.int(nrow(browser),100)] 
s <- sd(smallsamp) # sample variance
s
sd(browser$spend)   # full sample variance 
s/sd(browser$spend)  # sample variance /full sample variance ratio 

# 9. Bootstrap CI can also correct for bias
eb <- c()
for (b in 1:B){
  sb <- sd(smallsamp[sample.int(100, replace=TRUE)])
  eb <- c(eb, sb-s)
}

mean(eb)    # mean of error is not close to zero at all
mean(s-eb)    # but this is close to full sample sd 
sd(browswer$spend)   # full sample sd 

tvals <- quantile(eb, c(0.05, 0.95))   
tvals

s-tvals[2:1]     # 90% CI is centered around true sigma

sd(mub)*100

# 10. Regression 1 
summary( glm( log(spend) ~ broadband + anychildren, data=browser) )

# 11. Nonparametric bootstrap for coefficients betas
B <- 1000
betas <- c()
for (b in 1:B){
  samp_b = sample.int(nrow(browser), replace=TRUE)
  reg_b <- glm(log(spend) ~ broadband + anychildren, data=browser[samp_b,])
  betas <- rbind(betas, coef(reg_b))
}

head(betas)

cor(betas[,"broadband"], betas[,"anychildren"])   # little correlation 

# 12. Bootstrap sampling distribution 
xx <- seq(min(betas[,2]),max(betas[,2]),length=100)
par(mai=c(.8,.8,.2,.2))
hist(betas[,2], main="", xlab="broadband coefficient",
     col=8, border="grey90", freq=FALSE)
lines(xx, dnorm(xx, 0.55285, 0.04357), col="royalblue", lwd=1.5)

# 13. Bootstrap sampling distribution for multiplicative effect
par(mai=c(.8,.8,.2,.2))
hist(exp(betas[,2]), main="", xlab="broadband multiplier",
     col=8, border="grey90", freq=FALSE)

# 13. Regression: all variables except "id"
spendy <- glm( log(spend) ~ .-id, data=browser)
round(summary(spendy)$coef,2)

pval <- summary(spendy)$coef[-1,"Pr(>|t|)"]    # extract p-values 

# 14. BH algorithm for FDR
par(mai=c(.8,.8,.2,.2))
plot(sort(pval), bty="n", xlab="rank", ylab="p-values")
abline(a=0, b=.1/9)
points(sort(pval)[1:5], col=2, pch=20)

# 15. Graph for uniform distribution
par(mai=c(.8,.8,.2,.2))
plot(c(-1,0,0,1,1,2), c(0,0,1,1,0,0), ylim=c(0,1.5), xlim=c(-0.1,1.1),
     type="l", bty="n", xlab="U", ylab="probability density", main = "uniform pdf")

# 16. Expected uniform order statistics 
par(mai=c(.8,.8,.2,.2))
plot(1:9, (1:9)/10, ylim=c(0,1),
     pch=16, col="black", bty="n", ylab="p-value",
     xlab="order", main = "p-value order statistics")
points(1:9, sort(pval), pch=17, col=rgb(0,0,1,.5))
legend("topleft", bty="n",
       legend=c("expectation under null","observed"), pch=c(16,17),
       col=c("black",rgb(0,0,1,.5)))

# 17. 
library(data.table)
system.time(lipids <-  fread("jointGwasMc_LDL.txt"))
lipids <- as.data.frame(lipids)

pvals <- as.numeric(lipids[,'P-value'])
names(pvals) <- lipids[,'rsid']

# Draw histogram
hist(pvals, main='', xlab='p-values',freq=FALSE)

names(pvals)[order(pvals)[1:10]]

# FDR cutoff 
fdr_cut <- function(pvals, q){
  pvals <- pvals[!is.na(pvals)]
  N <- length(pvals)
  k <- rank(pvals, ties.method="min")
  alpha <- max(pvals[ pvals<= (q*k/N) ])
  return(alpha)
}

cutoff10 <- fdr_cut(pvals,q=.1)   # q=0.1
print(cutoff10)
print(sum(pvals<=cutoff10))

cutoff1 <- fdr_cut(pvals,q=.01)   # q=0.01
print(cutoff1)
print(sum(pvals<=cutoff1))

cutoff01 <- fdr_cut(pvals,q=.001)    # q=0.001
print(cutoff01)
print(sum(pvals<=cutoff01))

sig <- factor(pvals<=cutoff01)
o <- order(pvals)
N <- length(pvals)
plot(pvals[o], log="xy", col=c("grey60","red")[sig[o]], pch=20,
     ylab="p-values", xlab="tests ordered by p-value", main = 'FDR = 0.1%')
lines(1:N, 0.01*(1:N)/N)

# Bayesian example 
a <- 1
b <- 1
qgrid <- seq(0,1,length=500)
plot(qgrid, dbeta(qgrid, a, b), ylim=c(0,10), col=8,
     type="l", bty="n", lwd=1, xlab="q", ylab="posterior density")
K <- 20
lcols <- terrain.colors(K*1.1)[K:1]
lg <- c(1,5,10,15,20)
for(k in 1:K){
  n <- 5
  x <- rbinom(n, p=1/3, size=1)
  a <- a + sum(x)
  b <- b + n - sum(x)
  print(c(a,b))
  lines(qgrid, dbeta(qgrid,a,b), col=lcols[k], lwd=1)
}
legend("topright", bty="n", title="sample size",
       legend=c(0,n*lg), col=c(8,lcols[lg]), lwd=1)
