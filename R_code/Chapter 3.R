# -------------------------- Intro -------------------------- #
######## 1. Load Semiconductors data 
setwd("/Users/nicoleyin88/Documents/1. Big Data/3. Presentation/2. Data/")

# :SC: 
SC <- read.csv("semiconductor.csv")

# There're n= 1477 rows (observations), 200 signals (variables), 1 dependent variable FAIL
dim(SC)                 

######## 2. Full In-Sample Analysis
# (1) Logistic regression for all in-sample 
# :full: 
full <- glm(FAIL ~ ., data=SC, family=binomial)     

# (2) find in sample R2
R2 <- print(1 - full$deviance/full$null.deviance)     

# (3) Result 
# About 56% variation in failure versus success is explained by the 200 diagnostic signals. 

# (4) Histogram for p-values 
# -1 to drop the intercept, 4 to extract 4th column of p-values 
pvals <- summary(full)$coef[-1,4]         
hist(pvals, xlab="p-value", main="", col="lightblue")

# (5) Result 
# Distribution of p-value is not uniform.

######## 3. Cut In-Sample Analysis (25 Significant signals)
# (1) Plot FDR for the semiconductor p-values
# Create fdr_cut function, set FDR cutoff at 10% (q=0.1)
fdr_cut <- function(pvals, q=0.1){        
  pvals <- sort(pvals[!is.na(pvals)])
  N <- length(pvals)
  k <- rank(pvals, ties.method="min")
  alpha <- max(pvals[ pvals<= (q*k/(N+1)) ])
  plot(pvals, log="xy", xlab="order", main=sprintf("FDR of %g",q),
       ylab="p-value", bty="n", col=c(8,2)[(pvals<=alpha) + 1], pch=20)
  lines(1:N, q*(1:N)/(N+1))
  return(alpha)
}

fdr_cut(pvals)              

# (2) Result 
# cutoff p-value is 0.012, this gives 25 significant signals 

# (3) Logistic Regression for 25 significant signals 
# i) Extract significant p-values by cutoff 
# note 1: if put the line in (), it prints results
# note 2: which() function returns an array
( signif <- which(pvals <= 0.0121704339598325) )          
                                                          
# ii) Create X that contains only 25 significant signals :cutvar:  
cutvar <- c("FAIL", names(signif))
head(cutvar)

# iii) Run logistic regression for 25 significant signal 
# :cut:  
cut <- glm(FAIL ~ ., data = SC[,cutvar], family="binomial") 

# (4) New in-sample R2
(R2_25signals <- 1 - cut$deviance/cut$null.deviance)   

# (5) Result 
# The R2 from 25 signals is much smaller than the R2 from full in sample. 

# -------------------------- K-Fold Out-of-Sample (OSS) Validation -------------------------- #
######## 1. Create Deviance function 
# :deviance: 
deviance <- function(y, pred, family=c("gaussian","binomial")){
  family <- match.arg(family)
  if(family=="gaussian"){
    return( sum( (y-pred)^2 ) )
  }else{
    if(is.factor(y)) y <- as.numeric(y)>1
    return( -2*sum( y*log(pred) + (1-y)*log(1-pred) ) )
  } 
}

######## 2. Create R2 function 
# :R2:  
R2 <- function(y, pred, family=c("gaussian","binomial")){
  fam <- match.arg(family)
  if(fam=="binomial"){
    if(is.factor(y)){ y <- as.numeric(y)>1 }
  }
  dev <- deviance(y, pred, family=fam)
  dev0 <- deviance(y, mean(y), family=fam)
  return(1-dev/dev0)
} 

######## 3. 10-Fold OSS Validation
# (1) Set up folds K=10 
n <- nrow(SC)
K <- 10 

# (2) Create a vector of fold memberships (random order)
# :foldid: 
foldid <- rep(1:K, each=ceiling(n/K))[sample(1:n)]

# assign each observation a fold ID (range 1 ~ 10)
head(foldid)                 

# (2) Create an empty dataframe of results
# :OSS: 
OOS <- data.frame(full=rep(NA,K), cut=rep(NA,K))

# the dataframe has 2 columns: full and cut 
head(OOS)                    

# (3) Run OSS experiment 
# :train: 
# :rfull: 
# :rcut: 
for(k in 1:K){
  
  # train all but fold k
  train <- which(foldid!=k)            
  
  # full  regression
  rfull <- glm(FAIL~., data=SC, subset=train, family=binomial)   
  
  # cut regression
  rcut <- glm(FAIL~., data=SC[,cutvar], subset=train, family=binomial) 
  
  # type = "response" gives prediction in probability
  predfull <- predict(rfull, newdata=SC[-train,], type="response")       
  predcut <- predict(rcut, newdata=SC[-train,], type="response")
  
  # Get R2 for test data
  OOS$full[k] <- R2(y=SC$FAIL[-train], pred=predfull, family="binomial")
  OOS$cut[k] <- R2(y=SC$FAIL[-train], pred=predcut, family="binomial")
  
  # print progress for each fold 
  cat(k, " ")          
}

# (4) OSS gives full and cut model R2 for each K
OOS

# (5) Result
# There're negative R2

# (6) Box Plot for R2 from full and cut model 
par(mai=c(.9,.9,.1,.1))
boxplot(OOS, col="plum", ylab="R2", xlab="model", bty="n")

# (7) Check average R2 from full and cut model 
colMeans(OOS)

# (8) Result 
# It confirms that full model has negative average R2, while cut model has positive one. 
# Negative R2 happens if the fitted model performs worse than null model.
# In this case, cut model is better than full model. 
# Negative R2 happens a lot in OOS.
# OOS R2 matters,not IS R2. We should find other methods. 

# -------------------------- Regularization: Forward Stepwise -------------------------- #
######## 1. Forward Stepwise Regression
# (1) Create a null model
# :null:
# Null model is a starting point: y~1 means with only intercept
null <- glm(FAIL~1, data=SC)
summary(null)

# (2) Forward stepwise regression
# :fwd: 
# system.time() shows running time
# step() is the stepwise function
system.time(fwd <- step(null, scope=formula(full), dir="forward")) 

# (3) Find number of coefficients for fwd regression
# This takes a few minutes
length(coef(fwd))
coef(fwd)

# (4) Result 
# fwd gives 69 coefficients 

# -------------------------- Regularization: Lasso -------------------------- #
library(gamlr)

#### 1. Lasso 
# (1) Create X matrix (Dimension is 1477 x 200)
# [,-1] excludes the intercept 
scx <- sparse.model.matrix(FAIL ~ ., data=SC)[,-1]

# (2) Create y 
scy <- SC$FAIL

# (3) run Lasso
# :sclasso:
sclasso <- gamlr(scx, scy, family="binomial")
summary(sclasso)

# (4) Result 
# from summary, we can see seg24 has lowest AIC
# with 31 df (number of coefficients or parameters)
# with lambda = 0.0110416566, or log(lambda) = -4.50608 

# (4) Plot Lasso Path 
plot(sclasso)
  
#### 2. AICc
# (1) Extract Betas from sclasso
scbeta <- print(coef(sclasso))

# (2) Extract log(lambda)
log(sclasso$lambda[which.min(AICc(sclasso))])

# (3) Find total number of coefficients
sum(scbeta!=0) 

# (4) Result 
# AICc gives log(lambda) = -4.50608, total 31 coefficients 

######## 3. BIC 
# (1) BIC segment = 1
BICseg <- print(which.min(BIC(sclasso)))

# (2) BIC coefficients only has intercept
scb.bic <- print(coef(sclasso, s=BICseg))

# (3) BIC uses only 1 intercept
sum(scb.bic!=0)

# (4) CV-1se 
sccvl <- cv.gamlr(scx, scy, family="binomial", verb=TRUE)
plot(sccvl, bty="n")

# ----- p.261



#### re-run all code with glmnet
library(glmnet)   # similar to glm package? 
library(lars)     # similar to gamlr package? 


