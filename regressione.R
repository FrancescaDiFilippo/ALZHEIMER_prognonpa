setwd("C:/Users/E5440/Desktop/esami/Nonparametric statistics/Progetto/Alzheimer")
dataset<- read.csv("oasis_longitudinal.csv")
library(dplyr)
#EDUC: years of education 
#SES: socio economic status
#MMSE: mini mental state examination (0:30 30=stato mentale perfetto)
#CDR: clinical dementia rating (0:3 3=perdita di memoria grave)
#eTIV: extimated total intracranial volume 
#nWBV: Normalize Whole Brain Volume
#ASF: Atlas Scaling Factor (the volume-scaling factor required to match each individual to the atlas target)

Nondemented <- filter(dataset, Group=="Nondemented")
Demented <- filter(dataset, Group=="Demented")
Converted <- filter(dataset, Group=="Converted")
x11()
par(mfrow = c(1,2))
hist(Nondemented$MMSE)
hist(Demented$MMSE)

FEMMINE <- filter(dataset, M.F=="F")
MASCHI <- filter(dataset, M.F=="M")

x11()
pairs(dataset[,c(5,8:15)])
graphics.off()
x11()
plot(dataset$eTIV, dataset$ASF)
graphics.off()
x11()
plot(dataset$ASF,dataset$MMSE)

#MMSE vs ASF
library(ISLR2)
library(car)
Nondemented <- as.data.frame(scale(Nondemented[,c(11,15)]))
Demented <- as.data.frame(scale(Demented[,c(11,15)]))
x11()
par(mfrow = c(1,2))
plot(Nondemented$ASF,Nondemented$MMSE)
plot(Demented$ASF,Demented$MMSE)
MMSE.n <- Nondemented$MMSE
ASF.n <- Nondemented$ASF
m_list.n=lapply(1:5,function(degree){lm(MMSE.n ~ poly(ASF.n,degree=degree))})
do.call(anova,m_list.n)

asf.grid.n=seq(range(Nondemented$ASF)[1],range(Nondemented$ASF)[2],by=0.005)
preds.n=predict(m_list.n[[2]],list(ASF.n=asf.grid.n),se=T)
se.bands.n=cbind(preds.n$fit +2* preds.n$se.fit ,preds.n$fit -2* preds.n$se.fit)
x11()
plot(ASF.n ,MMSE.n, xlim=range(asf.grid.n) ,cex =.5, col =" darkgrey ",main='Degree 2 Poly - Fit')
lines(asf.grid.n,preds.n$fit ,lwd =2, col =" blue")
matlines (asf.grid.n ,se.bands.n ,lwd =1, col =" blue",lty =3)

MMSE.d <- Demented$MMSE
ASF.d <- Demented$ASF
m_list.d=lapply(1:5,function(degree){lm(MMSE.d ~ poly(ASF.d,degree=degree))})
do.call(anova,m_list.d)

asf.grid.d=seq(range(Demented$ASF)[1],range(Demented$ASF)[2],by=0.005)
preds.d=predict(m_list.d[[3]],list(ASF.d=asf.grid.d),se=T)
se.bands.d=cbind(preds.d$fit +2* preds.d$se.fit ,preds.d$fit -2* preds.d$se.fit)
x11()
plot(ASF.d ,MMSE.d ,xlim=range(asf.grid.d) ,cex =.5, col =" darkgrey ",main='Degree 4 Poly - Fit')
lines(asf.grid.d,preds.d$fit ,lwd =2, col =" blue")
matlines (asf.grid.d ,se.bands.d ,lwd =1, col =" blue",lty =3)

#multidim:  MMSE vs EDUC + nBWV
library(mgcv)
library(rgl)
library(splines)
library(pbapply)

x11()

converted <- which(dataset$Group=='Converted')
conver.data <- dataset[converted,]
train <- dataset[-converted,]
demented <- rep(0:dim(train)[1])
demented[which(train$Group=='Demented')] <- 1
demented <- factor(demented)


model_gam=gam(demented ~ s(EDUC,bs='cr') + s(nWBV,bs='cr') + s(Age, bs='cr')  + s(MMSE, bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam)

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
  pred <- predict(model_gam, newdata = conver.data, type = 'response')
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






x11()
par(mfrow = c(1,2))
plot(model_gam)
graphics.off()

educ.grid=seq(range(dataset$EDUC)[1],range(dataset$EDUC)[2],length.out = 100)
nwbv.grid=seq(range(dataset$nWBV)[1],range(dataset$nWBV)[2],length.out = 100)
grid=expand.grid(educ.grid,nwbv.grid)
colnames(grid) <- c('EDUC', 'nWBV')
pred_gam=predict(model_gam,newdata=grid)
persp3d(educ.grid,nwbv.grid,pred_gam,col='grey30')
with(dataset,points3d(EDUC,nWBV,MMSE,col='black',size=5))
