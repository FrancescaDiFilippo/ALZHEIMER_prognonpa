remove(list=ls())
setwd("C:/Users/E5440/Desktop/esami/Nonparametric statistics/Progetto/Alzheimer")
library(dplyr)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)
library(ISLR2)
library(car)
library(ggplot2)
library(yardstick)

dataset<- read.csv("oasis_longitudinal.csv")
dataset <- filter(dataset, Visit==1)

#EDUC: years of education 
#SES: socio economic status
#MMSE: mini mental state examination (0:30 30=stato mentale perfetto. Anche persone con un iniziale deterioramento cognitivo, ma con un'alta scolarizzazione possono ottenere un punteggio pari a 29 e 30)
#CDR: clinical dementia rating (0:3 3=perdita di memoria grave)
#eTIV: extimated total intracranial volume 
#nWBV: Normalize Whole Brain Volume
#ASF: Atlas Scaling Factor (the volume-scaling factor required to match each individual to the atlas target)

data <- dataset[,c(6,8,9,11,13,14)]
i.converted <- which(dataset$Group=='Converted')
data <- data[-i.converted,]
data.aux <- dataset[-i.converted,]
data$label <- rep('Nondem', dim(data)[1])
i.dem <- which(data.aux$Group=='Demented')
data$label[i.dem] <- 'Dem'
data$label <- factor(data$label, labels = c('Nondem', 'Dem'))
data$M <- ifelse(data$M.F=='M',1,0)

dataset2 <- read.table('C:/Users/E5440/Downloads/test set.csv', header = TRUE, sep=';')
dataset2$label <- ifelse(dataset2$CDR > 0, 'Dem', 'Nondem')
dataset2$label <- factor(dataset2$label, labels = c('Nondem', 'Dem'))
dataset2  <- dataset2[,c(2,4,5,7,9,10,12)]
dataset2$M <- ifelse(dataset2$M.F=='M',1,0)
colnames(dataset2) <- colnames(data)
data <- rbind(data, dataset2) 

train <- data[1:271,]
test <- data[272:dim(data)[1],]

dem <- rep(0, dim(data)[1])
dem[which(data$label=='Dem')] <- 1

x11()
with(data = data, scatterplotMatrix(data.frame(Age, EDUC, MMSE, eTIV, nWBV), groups = label))

x11()
with(data = data, scatterplotMatrix(data.frame(Age, EDUC, MMSE, eTIV, nWBV), groups = M.F))

#se le distribuzioni sono simili nelle diverse label, la dummy non è significativa, se sono shiftate ma la forma è uguale la dummy influeisce solo nell'intercetta

model_gam=gam(label ~ M + s(EDUC,bs='cr') + s(nWBV,bs='cr') + s(Age, bs='cr')  + s(MMSE, bs='cr') + s(eTIV, bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam)

#model_gam2 <- gam(label ~ M + s(EDUC,bs='cr') + s(I(EDUC*M),bs='cr') + s(nWBV,bs='cr') + s(I(nWBV*M),bs='cr') + s(Age, bs='cr') + s(I(Age*M), bs='cr') + s(MMSE, bs='cr') + s(I(MMSE*M), bs='cr') + s(eTIV, bs='cr') + s(I(M*eTIV), bs='cr'), data = train, select = TRUE, family = binomial)
#summary(model_gam2) #non va bene

model_gam2=gam(label ~ M + s(Age, bs='cr') + s(nWBV,bs='cr') + s(MMSE, bs='cr') + s(eTIV, bs='cr'), data = train, select = TRUE, family = binomial)
summary(model_gam2)#senza educ

#sembra che etiv possa essere modellizzata con la dummy M.F
model_gam3=gam(label ~ M + s(Age, bs='cr') + s(nWBV,bs='cr') + s(MMSE, bs='cr') + s(eTIV, bs='cr') + s(I(eTIV*M), bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam3)#senza educ, con dummy per etiv
anova.gam(model_gam2, model_gam3, test = 'Chisq') #perchè il  pv non è uguale al pv del summary?

model_gam4=gam(label ~ M + s(Age, bs='cr') + s(nWBV,bs='cr') + s(MMSE, bs='cr'), data = train, select = TRUE, family = binomial)
summary(model_gam4)#senza educ e etiv
anova.gam(model_gam2, model_gam4, test = 'Chisq')
model_gam2$aic
model_gam4$aic
#il migliore sembra model_gam2
roc.curve <- function(predicted, test.set){
  p0 <- seq(0,1,by=0.001)
  spec <- NULL
  sens <- NULL
  for (i in 1:length(p0)) {
    predicted$class.assigned <- rep('Nondem', dim(predicted)[1])
    colnames(predicted) <- c('prob','class.assigned')
    predicted[which(predicted$prob>=p0[i]),2] <- 'Dem'
    true.lab <- test.set$label
    i.equal <- which(true.lab==predicted$class.assigned)
    n.equal <- length(i.equal)
    test.equal <- test.set[i.equal,]
    n00 <- length(which(test.equal$label=='Nondem'))
    n11 <- n.equal - n00
    
    n01 <- length(which(predicted$class.assigned=='Nondem'&true.lab=='Dem')) #classified,true lab
    n10 <- length(which(predicted$class.assigned=='Dem'&true.lab=='Nondem')) #classified,true lab
    
    sensitivity <- n11/(n01+n11)
    specificity <- n00/(n00+n10)
    spec <- c(spec, specificity)
    sens <- c(sens, sensitivity)
  }
  x.roc <- rep(1,length(spec))-spec
  y.roc <- sens
  data.frame(cbind(x.roc,y.roc))
}

pred <- predict(model_gam2, newdata = test, type = 'response')
pred <- as.data.frame(pred)

#ROC curve
#NB è PIù IMPORTANTE MASSIMIZZARE SENSITIVITà
#model_gam2
data.roc <- roc.curve(pred, test)
x11()
plot.roc <- ggplot(data.roc,aes(x.roc,y.roc))+
  theme(panel.background = element_rect(fill = 'gray90', colour = 'white'))+
  xlim(c(0,1)) + ylim(c(0,1))+
  geom_path(col='dodgerblue4', lwd=1.5)+ #geom_path
  geom_ribbon(aes(ymin = 0, ymax=y.roc),fill = 'seagreen2')+
  geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
  xlab('1 - specificity') + ylab('sensitivity')

roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second') #da aggiungere nel plot

#model_gam4 (quello con solo age, MMSE, nWBV)
pred <- predict(model_gam4, newdata = test, type = 'response')
pred <- as.data.frame(pred)
data.roc <- roc.curve(pred,test)
x11()
ggplot(data.roc,aes(x.roc,y.roc))+
  theme(panel.background = element_rect(fill = 'gray90', colour = 'white'))+
  xlim(c(0,1)) + ylim(c(0,1))+
  geom_line(col='dodgerblue4', lwd=1.9)+
  geom_ribbon(aes(ymin=0, ymax=y.roc), fill = 'seagreen2')+
  geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
  xlab('1 - specificity')+ylab('sensitivity')
roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second')#da aggiungere nel plot

#provo solo con MMSE
model_gam_MMSE <- gam(label ~ M + s(MMSE, bs='cr') + s(I(MMSE*M), bs='cr'), data = train, select = TRUE, family = binomial)
summary(model_gam_MMSE)#senza educ e etiv
pred <- predict(model_gam_MMSE, newdata = test, type = 'response')
pred <- as.data.frame(pred)
data.roc <- roc.curve(pred,test)
x11()
ggplot(data.roc,aes(x.roc,y.roc))+
  theme(panel.background = element_rect(fill = 'gray90', colour = 'white'))+
  xlim(c(0,1)) + ylim(c(0,1))+
  geom_line(col='dodgerblue4', lwd=1.9)+
  geom_ribbon(aes(ymin=0, ymax=y.roc), fill = 'seagreen2')+
  geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
  xlab('1 - specificity')+ylab('sensitivity')
roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second')#da aggiungere nel plot

##############################################################################################################################################################
##############################################################################################################################################################

#GLM MODEL
library(tidyverse)
library(caret)
library(leaps)
library(MASS)
library(RcmdrMisc)

model_glm=glm(label ~ M + Age + nWBV + MMSE + eTIV, family = 'binomial' , data = train)
summary(model_glm)

1-model_glm$deviance/model_glm$null.deviance


#MODEL SELECTION 
# mod <- stepwise(model_glm, trace = 1,criterion='BIC' ) #criterion=BIC. quello che è considerato come Aic in realtà è Bic, Aik=205.02
# summary(mod)
# 1-mod$deviance/mod$null.deviance
# mod$deviance-model_glm$deviance
# mod$deviance/model_glm$deviance 

pred <- predict(model_glm, newdata = test, type = 'response')
pred <- as.data.frame(pred)
data.roc <- roc.curve(pred,test)
x11()
ggplot(data.roc,aes(x.roc,y.roc))+
  theme(panel.background = element_rect(fill = 'gray90', colour = 'white'))+
  xlim(c(0,1)) + ylim(c(0,1))+
  geom_line(col='dodgerblue4', lwd=1.9)+
  geom_ribbon(aes(ymin=0, ymax=y.roc), fill = 'seagreen2')+
  geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
  xlab('1 - specificity')+ylab('sensitivity')
roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second')#da aggiungere nel plot


model_glm2=glm(label ~ Age + nWBV + MMSE, family = 'binomial' , data = train)
summary(model_glm2)

anova(model_glm, model_glm2, test = 'Chisq')

pred <- predict(model_glm2, newdata = test, type = 'response')
pred <- as.data.frame(pred)
data.roc <- roc.curve(pred,test)
x11()
ggplot(data.roc,aes(x.roc,y.roc))+
  theme(panel.background = element_rect(fill = 'gray90', colour = 'white'))+
  xlim(c(0,1)) + ylim(c(0,1))+
  geom_line(col='dodgerblue4', lwd=1.9)+
  geom_ribbon(aes(ymin=0, ymax=y.roc), fill = 'seagreen2')+
  geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
  xlab('1 - specificity')+ylab('sensitivity')
roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second')#da aggiungere nel plot


######################################################################################################################
######################################################################################################################

