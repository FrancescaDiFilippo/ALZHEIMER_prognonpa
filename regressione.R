setwd("C:/Users/E5440/Desktop/esami/Nonparametric statistics/Progetto/Alzheimer")
dataset<- read.csv("oasis_longitudinal.csv")
dataset <- dataset[-which(is.na(dataset$MMSE)),]
library(dplyr)
#EDUC: years of education 
#SES: socio economic status
#MMSE: mini mental state examination (0:30 30=stato mentale perfetto)
#CDR: clinical dementia rating (0:3 3=perdita di memoria grave)
#eTIV: extimated total intracranial volume 
#nWBV: Normalize Whole Brain Volume
#ASF: Atlas Scaling Factor (the volume-scaling factor required to match each individual to the atlas target)



#multidim:  MMSE vs EDUC + nBWV
library(mgcv)
library(rgl)
library(splines)
library(pbapply)


converted <- which(dataset$Group=='Converted')
conver.data <- dataset[converted,]
train <- dataset[-converted,-c(10,15)]
demented <- rep(0, dim(train)[1])
demented[which(train$Group=='Demented')] <- 1
demented <- factor(demented)


model_gam=gam(demented ~ s(EDUC,bs='cr') + s(nWBV,bs='cr') + s(Age, bs='cr')  + s(MMSE, bs='cr') + s(eTIV) , data = train, select = TRUE, family = binomial)
summary(model_gam)

model_gam2 <- gam(demented ~ s(EDUC,bs='cr') + s(nWBV,bs='cr') + s(Age, bs='cr')  + s(MMSE, bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam2)

anova(model_gam2, model_gam)

pred <- predict(model_gam, newdata = conver.data, type = 'response')
pred <- as.data.frame(pred)
pred$type <- rep('Nondemented', dim(pred)[1])
colnames(pred) <- c('prob','type')
pred[which(pred$prob>=0.5),2] <- 'Demented'

true.lab <- rep('Nondemented', dim(pred)[1])
true.lab[which(conver.data$CDR>0)] <- 'Demented'
i.equal <- which(true.lab==pred$type)
n.equal <- length(i.equal)
conver.equal <- conver.data[i.equal,]
n00 <- length(which(conver.equal$CDR==0))
n11 <- n.equal-n00

n01 <- length(which(pred$type=='Nondemented'&true.lab=='Demented')) #classified,true lab
n10 <- length(which(pred$type=='Demented'&true.lab=='Nondemented')) #classified,true lab

sensitivity <- n11/(n01+n11)
specificity <- n00/(n00+n10)

#ROC curve
p0 <- seq(0,1,by=0.001)
spec <- NULL
sens <- NULL
for (i in 1:length(p0)) {
  pred <- predict(mod, newdata = conver.data, type = 'response')
  pred <- as.data.frame(pred)
  pred$type <- rep('Nondemented', dim(pred)[1])
  colnames(pred) <- c('prob','type')
  pred[which(pred$prob>=p0[i]),2] <- 'Demented'
  true.lab <- rep('Nondemented', dim(pred)[1])
  true.lab[which(conver.data$CDR>0)] <- 'Demented'
  i.equal <- which(true.lab==pred$type)
  n.equal <- length(i.equal)
  conver.equal <- conver.data[i.equal,]
  n00 <- length(which(conver.equal$CDR==0))
  n11 <- n.equal-n00
  
  n01 <- length(which(pred$type=='Nondemented'&true.lab=='Demented')) #classified,true lab
  n10 <- length(which(pred$type=='Demented'&true.lab=='Nondemented')) #classified,true lab
  
  sensitivity <- n11/(n01+n11)
  specificity <- n00/(n00+n10)
  spec <- c(spec, specificity)
  sens <- c(sens, sensitivity)
}
x11()
plot(rep(1,length(spec))-spec,sens)
lines(rep(1,length(spec))-spec,sens)
points(seq(0,1,length.out =length(spec)), seq(0,1,length.out =length(spec)), col='red')

i.bar <- which(sens>=0.6 & spec>=0.5)
i.bar <- 122
p <- p0[i.bar]
pred[which(pred$prob>=p),2] <- 'Demented'

#non funziona
mod.red <- stepwise(model_gam, trace = 1, criterion='BIC') #criterion=BIC. quello che è considerato come Aic in realtà è Bic, Aik=205.02



#GLM MODEL
library(tidyverse)
library(caret)
library(leaps)
library(MASS)
library(RcmdrMisc)

model_glm=glm(demented ~ M.F + EDUC + EDUC:M.F + nWBV + nWBV:M.F + Age + Age:M.F  + MMSE + MMSE:M.F +eTIV + eTIV:M.F , data = train, family = binomial)
summary(model_glm)

1-model_glm$deviance/model_glm$null.deviance


#MODEL SELECTION 
mod <- stepwise(model_glm, trace = 1,criterion='BIC' ) #criterion=BIC. quello che è considerato come Aic in realtà è Bic, Aik=205.02
summary(mod)
1-mod$deviance/mod$null.deviance
mod$deviance-model_glm$deviance
mod$deviance/model_glm$deviance 



pred <- predict(mod, newdata = conver.data, type = 'response')
pred <- as.data.frame(pred)
pred$type <- rep('Nondemented', dim(pred)[1])
colnames(pred) <- c('prob','type')
true.lab <- rep('Nondemented', dim(conver.data)[1])
true.lab[which(conver.data$CDR>0)] <- 'Demented'
p0 <- seq(0,1,by=0.001)
spec <- NULL
sens <- NULL
for (i in 1:length(p0)) {
  pred <- predict(mod, newdata = conver.data, type = 'response')
  pred <- as.data.frame(pred)
  pred$type <- rep('Nondemented', dim(pred)[1])
  colnames(pred) <- c('prob','type')
  pred[which(pred$prob>=p0[i]),2] <- 'Demented'
  true.lab <- rep('Nondemented', dim(pred)[1])
  true.lab[which(conver.data$CDR>0)] <- 'Demented'
  i.equal <- which(true.lab==pred$type)
  n.equal <- length(i.equal)
  conver.equal <- conver.data[i.equal,]
  n00 <- length(which(conver.equal$CDR==0))
  n11 <- n.equal-n00
  
  n01 <- length(which(pred$type=='Nondemented'&true.lab=='Demented')) #classified,true lab
  n10 <- length(which(pred$type=='Demented'&true.lab=='Nondemented')) #classified,true lab
  
  sensitivity <- n11/(n01+n11)
  specificity <- n00/(n00+n10)
  spec <- c(spec, specificity)
  sens <- c(sens, sensitivity)
}
x11()
plot(rep(1,length(spec))-spec,sens)
lines(rep(1,length(spec))-spec,sens)
points(seq(0,1,length.out =length(spec)), seq(0,1,length.out =length(spec)), col='red')



######################################################################################################################
######################################################################################################################

