install.packages("MASS")
install.packages("epiDisplay")

library(MASS)
library(epiDisplay)
#----------------------------------------#
########## Question 1 ##########
#----------------------------------------#
#Code needed for Question 1f
questionOneTable <- as.table(matrix(c(80,261,771,1044), ncol=2, byrow=TRUE))
Q1chisq <- chisq.test(questionOneTable)
Q1chisq
#Code needed for Question 1g
pchisq(256.026, df=1,lower.tail=FALSE)
#----------------------------------------#
########## Question 2 ##########
#read grouped and ungrouped datasets
Q2url1 <- url("https://ajmolstad.github.io/docs/Q2Grouped.RDS")
Q2Grouped <- readRDS(Q2url1)
Q2url2 <- url("https://ajmolstad.github.io/docs/Q2Ungrouped.RDS")
Q2Ungrouped <- readRDS(Q2url2)
#Code needed for Question 2 (a)-(d)
Q2UngroupedModel <- glm(y~ x1+x2+(x1:x2)+x3, family=binomial(link = "logit"), data = Q2Ungrouped)
summary(Q2UngroupedModel)
#Created linear model with grouped data
Q2GroupedModel <- glm(cbind(success,failure) ~ x1+x2+x1:x2+x3, family = binomial, data = Q2Grouped)
summary(Q2GroupedModel)
1-pchisq(Q2GroupedModel$deviance,Q2GroupedModel$df.residual)
#Created Models using only x1 and x2 to compare deviances
ungroupedLRTModel <- glm(y~ x1+x2, family = binomial, data = Q2Ungrouped)
groupedLRTModel <- glm(cbind(success,failure)~ x1+x2, family = binomial, data=Q2Grouped)
#conducted likelihood ratio test using the models
ungroupedTestStatistic <- ungroupedLRTModel$deviance - Q2UngroupedModel$deviance
1-pchisq(ungroupedTestStatistic, 2)
groupedTestStatistic <- groupedLRTModel$deviance - Q2GroupedModel$deviance
1-pchisq(groupedTestStatistic, 2)
#----------------------------------------#
########## Question 3 ##########
#----------------------------------------#
inputURL <- url("https://ajmolstad.github.io/docs/Cancer.RDS")
Cancer <- readRDS(inputURL)
#code used for question 3a
SNPmargModel <- glm(y~as.factor(SNPcat), family=binomial, data = Cancer)
summary(SNPmargModel)
q3aTestStat <- SNPmargModel$null.deviance-SNPmargModel$deviance
1-pchisq(q3aTestStat, 2)
#model needed for question 3b-g
threeBModel <- glm(y~as.factor(SNPcat)+age+smoking, family=binomial, data=Cancer)
summary(threeBModel)
#Interval for 3c
q3Int=predict(threeBModel,newdata=data.frame(SNPcat=0, age=80, smoking=0),type="link",se.fit=TRUE)
q3Int.ci <- q3Int$fit+c(-1,1)*qnorm(0.975)*q3Int$se.fit
q3Int.ci
#calculation for change in odds ratio NOTE: did not do change from 50 to 60 by hand
#because it doesn't change between x and x+10 (noted in lectures and done in homework once)
exp(10*(-0.001812013))
#probit regression model
threeFModel <- glm(y~as.factor(SNPcat)+age+smoking, family=binomial(link=probit), data=Cancer)
summary(threeFModel)
q3Fifty <- predict(threeFModel, newdata=data.frame(SNPcat=2, age=50, smoking=10), type="link")
q3Sixty <-predict(threeFModel, newdata=data.frame(SNPcat=2, age=60, smoking=10), type="link")
p1 <- pnorm(q3Sixty)
p2 <- pnorm(q3Fifty)
ProbitOddsRatio <- (p1*(1-p2))/(p2*(1-p1))
ProbitOddsRatio
#logit model using SNPnum
threeGModel <- glm(y~SNPnum+age+smoking, family=binomial, data=Cancer)
summary(threeGModel)
#making models for (SNPcat, age), (SNPnum), and (SNPnum, age)
#since others have previously been made
catAgeModel <- glm(y~as.factor(SNPcat)+age, family=binomial, data=Cancer)
numModel <- glm(y~SNPnum, family=binomial, data=Cancer)
numAgeModel <- glm(y~SNPnum+age, family = binomial, data=Cancer)
AIC(SNPmargModel)
BIC(SNPmargModel)
AIC(catAgeModel)
BIC(catAgeModel)
AIC(threeBModel)
BIC(threeBModel)
AIC(numModel)
BIC(numModel)
AIC(numAgeModel)
BIC(numAgeModel)
AIC(threeGModel)
BIC(threeGModel)
#ROC Curves
#for some reason I keep getting "invalid graphics state errors, I use dev.off()
#to solve it, but I don't know if it's just a problem with my computer or not
ROC1 = lroc(threeGModel,grid=FALSE,graph=TRUE,line.col=2)
title("ROC curve")
ROC1$auc
#Code to conduct CCR
predicted.probs <- predict(threeGModel, type = "response")
pi_0 <- 0.5
hatY <- 1*(predicted.probs > pi_0)
cbind(hatY, predicted.probs)[1:5,]
ccr0.5 <- (sum(hatY == 1 & Cancer$y == 1) + sum(hatY == 0 & Cancer$y == 0))/length(Cancer$y)
ccr0.5
