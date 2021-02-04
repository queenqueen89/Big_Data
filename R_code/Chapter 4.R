#### 2. Nearest Neighbors for Forensic Glass
# load the data 
library(MASS)
data(fgl)
head(fgl)   # 9 types of elements (RI,Na, etc) and glass "type"
dim(fgl)

# boxplot 
options(repr.plot.width=14,repr.plot.antialias='subpixel',repr.plot.res=218)

par(mfrow=c(2,3))
plot(RI ~ type, data=fgl, col=c(grey(.2),2:6))
plot(Al ~ type, data=fgl, col=c(grey(.2),2:6))
plot(Na ~ type, data=fgl, col=c(grey(.2),2:6))
plot(Mg ~ type, data=fgl, col=c(grey(.2),2:6))
plot(Ba ~ type, data=fgl, col=c(grey(.2),2:6))
plot(Si ~ type, data=fgl, col=c(grey(.2),2:6))

# Scale elements (exclude 10th column for glass "type")
head(fgl[,1:9])

x <- scale(fgl[,1:9])
head(x)

apply(x,2,sd)

# 1-NN, 5-NN
library(class)

# Create testing set (with 14 samples)
set.seed(1234)
test <- sample(1:214,14)
test

# 1-NN
nearest1 <- knn(train=x[-test,], test=x[test,], cl=fgl$type[-test], k=1)

# 5-NN
nearest5 <- knn(train=x[-test,], test=x[test,], cl=fgl$type[-test], k=5)

# Compare results 
data.frame(fgl$type[test],nearest1,nearest5) 

##### 3. Credit 
# load the data 
setwd("/Users/nicoleyin88/Documents/1. Big Data/3. Presentation/2. Data/")
credit <- read.csv("credit.csv")

dim(credit)
head(credit)

# 2. factor and label the levels of credit history
credit$history = factor(credit$history, 
                        levels=c("A30","A31","A32","A33","A34"))
levels(credit$history) = c("good","good","poor","poor","terrible")

credit$foreign <- factor(credit$foreign, levels=c("A201","A202"),
                         labels=c("foreign","german"))
credit$rent <- factor(credit$housing=="A151")
credit$purpose <- factor(credit$purpose, 
                         levels=c("A40","A41","A42","A43","A44",
                                  "A45","A46","A47","A48","A49","A410")) 
levels(credit$purpose) <- c("newcar","usedcar",
                            rep("goods/repair",4),"edu",NA,"edu","biz","biz")

credit <- credit[,c("Default", "duration", "amount",
                    "installment", "age", "history",
                    "purpose", "foreign", "rent")]

par(mai=c(.8,.8,.1,.1))
plot(factor(Default) ~ history, 
     data=credit, col=c(8,2), ylab="Default")

# Create naref function
xnaref <- function(x){
  if(is.factor(x))
    if(!is.na(levels(x)[1]))
      x <- factor(x,levels=c(NA,levels(x)),exclude=NULL)
    return(x) }

naref <- function(DF){
  if(is.null(dim(DF))) return(xnaref(DF))
  if(!is.data.frame(DF))
    stop("You need to give me a data.frame or a factor")
  DF <- lapply(DF, xnaref)
  return(as.data.frame(DF))
}

# Create sparse matrix 
library(gamlr)
set.seed(1234)
credx <- sparse.model.matrix(Default ~ .^2, data=naref(credit))[,-1]
default <- credit$Default

# Logistic lasso regression 
credscore <- cv.gamlr(credx, default, family="binomial", verb=TRUE)

# Plot  
par(mfrow=c(1,2))
plot(credscore$gamlr)
plot(credscore)

# Compare AIC, AICc, BIC
# AICc
sum(coef(credscore, s="min")!=0) # min sum(coef(credscore$gamlr)!=0) 

# cv-1se 
sum(coef(credscore)!=0)  

# AIC 
sum(coef(credscore$gamlr, s=which.min(AIC(credscore$gamlr)))!=0) 

# BIC
sum(coef(credscore$gamlr, s=which.min(BIC(credscore$gamlr)))!=0) 

# OOS R2
1 - credscore$cvm[credscore$seg.min]/credscore$cvm[1]

# get IS probabilities 
pred <- predict(credscore$gamlr, credx, type="response")

# remove the sparse Matrix formatting
pred <- drop(pred) 
head(pred)

# Boxplot
boxplot(pred ~ default, xlab="default", 
        ylab="prob of default", col=c("pink","dodgerblue"))

# set cutoff (rule) p=1/5
rule <- 1/5
sum((pred>rule)[default==0] )/sum(pred>rule)
sum((pred<rule)[default==1] )/sum(pred<rule)

# set cutoff (rule) p=1/2
rule2 <- 1/2
sum((pred>rule2)[default==0] )/sum(pred>rule2)
sum((pred<rule2)[default==1] )/sum(pred<rule2)

## sensitivity 
sum( (pred>rule)[default==1] )/sum(default==1) 

## specificity
sum( (pred<rule)[default==0] )/sum(default==0) 

##### 4. Multinomial Logistic Regression with Forensic Glass
library(glmnet)

# create sparse matrix 
xfgl <- sparse.model.matrix(type~.*RI, data=fgl)[,-1]
gtype <- fgl$type

# lasso regression 
glassfit <- cv.glmnet(xfgl, gtype, family="multinomial")

plot(glassfit)