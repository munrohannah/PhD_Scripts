setwd("D:/Rfiles/AI") #how to set my working directory
library(lattice)#To calculate correlation matrix
library(fmsb)#to calculate Negelkerke's R^2
library(bbmle)
library(AICcmodavg)

data<-read.csv("RawGull09.11.csv",header=TRUE) #how to import a file, csv file, and has a header
attach(data)

head(data) # gives the first few lines of data
str(data) #check the data type and varibles
data$Year<-as.factor(data$Year)
data$Pos.Neg<-as.factor(data$Pos.Neg)

CORR <- NULL
CORR <- rbind(CORR,data.frame(Sample.Type = data$Sample.Type))
CORR$Species <- data$Species
CORR$Location <- data$Location
CORR$Year <- data$Year
CORR$Season <- data$Season

head(CORR)
splom(CORR)

Cand.models <- list( )
Cand.models[[1]] <- glm(Pos.Neg ~ Season + Sample.Type, family = binomial, data=data)
Cand.models[[2]] <- glm(Pos.Neg ~ Age ,family = binomial, data=data)
Cand.models[[3]] <- glm(Pos.Neg ~ Year + Sample.Type ,family = binomial, data=data)
Cand.models[[4]] <- glm(Pos.Neg ~ Species, family = binomial, data=data)
Cand.models[[5]] <- glm(Pos.Neg ~ Age+Species, family = binomial, data=data)
Cand.models[[6]] <- glm(Pos.Neg ~ Sample.Type ,family = binomial, data=data)
Cand.models[[7]] <- glm(Pos.Neg ~ Year,family = binomial, data=data)
Cand.models[[8]] <- glm(Pos.Neg ~ Season ,family = binomial, data=data)
Cand.models[[9]] <- glm(Pos.Neg ~ Age + Season ,family = binomial, data=data)



##create a vector of names to trace back models in set
Modnames <- c("GLM_1", "GLM_2", "GLM_3", "GLM_4","GLM_5", "GLM_6" , "GLM_7" , "GLM_8", "GLM_9")

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)


##find parameters for model 2
glm2<- glm(Pos.Neg ~ Season + Sample.Type,family = binomial, data=data)
glm2
anova(glm2)


citation() 
R.Version()