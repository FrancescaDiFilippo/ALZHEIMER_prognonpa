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
data$label <- rep('Dem', dim(data)[1])
i.nondem <- which(data.aux$Group=='Nondemented')
data$label[i.nondem] <- 'Nondem'
data$label <- factor(data$label)
data$M <- ifelse(data$M.F=='M',1,0)

dataset2 <- read.table('C:/Users/E5440/Downloads/test set.csv', header = TRUE, sep=';')
dataset2$label <- ifelse(dataset2$CDR > 0, 'Dem', 'Nondem')
dataset2  <- dataset2[,c(2,4,5,7,9,10,12)]
dataset2$M <- ifelse(dataset2$M.F=='M',1,0)
colnames(dataset2) <- colnames(data)
data <- rbind(data, dataset2)
data <- data[,c(1,4,7,8)]

train <- data[1:250,]
test <- data[251:dim(data)[1],]

dem <- rep(0, dim(data)[1])
dem[which(data$label=='Dem')] <- 1

#GLM WITH 1 COVARIATE
mdl <- glm(label ~ MMSE, data = train, family = binomial)
summary(mdl)

x <- data.frame(MMSE = seq(range(train$MMSE)[1], range(train$MMSE)[2], by=0.25))
p.hat <- predict(mdl, newdata = x, type = 'response') #1/(1+exp(-mdl$coefficients[1]-mdl$coefficients[2]*x))
x11()
ggplot()+
  geom_line(aes(x$MMSE,p.hat), col='red')+
  geom_point(aes(train$MMSE, as.numeric(train$label=='Nondem')))


#polynomial regression
mdl_poly=lapply(1:8,function(degree){glm(label ~ poly(MMSE,degree=degree), data = train, family = binomial)})
do.call(what = anova, args = mdl_poly)
mdl_poly[[1]]$aic
mdl_poly[[2]]$aic
mdl_poly[[3]]$aic
mdl_poly[[4]]$aic
mdl_poly[[5]]$aic
mdl_poly[[6]]$aic
mdl_poly[[7]]$aic
mdl_poly[[8]]$aic
summary(mdl_poly[[1]])
summary(mdl_poly[[5]])
auc_poly <- NULL
for (i in 1:8) {
  p.hat <- predict(mdl_poly[[i]], newdata = test, type = 'response')
  auc_poly <- c(auc_poly,roc_auc_vec(truth = test$label, estimate = rep(1,dim(test)[1])-p.hat))
}












