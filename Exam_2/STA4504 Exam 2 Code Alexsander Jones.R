#Libraries
library(lme4)
library(gee)
library(VGAM)

#Q1data
pathQ1 <- url("https://ajmolstad.github.io/docs/Arthritis.RDS")
RAdata <- readRDS(pathQ1)
#Code for Baseline-Category Logit Model, note that reference group is "Worsened"
Q1blogit <- vglm(Symptoms ~ Events, family=multinomial, data=RAdata, hde.NA= FALSE)
summary(Q1blogit)
#Likelihood Ratio Test
Q1trashBlogit <- vglm(Symptoms ~ 1, family=multinomial, data=RAdata, hde.NA= FALSE)
lrtest(Q1blogit, Q1trashBlogit)
#Confidence Interval
vcov(Q1blogit)
Q1CstandardError <- sqrt(0.1904763+0.2725196-2*0.1904763)
-0.5819+1.96*Q1CstandardError
-0.5819-1.96*Q1CstandardError
#Cumulative Logit Model
Q1Clogit <- vglm(Symptoms ~ Events, family = cumulative(parallel = TRUE),
                 data = RAdata)
summary(Q1Clogit)
Q1trashClogit <- vglm(Symptoms ~ 1, family = cumulative(parallel = TRUE),
                      data = RAdata)
lrtest(Q1Clogit, Q1trashClogit)
#Confidence Interval for clogit model
vcov(Q1Clogit)
Q1HstandardError <- sqrt(0.04546813+0.07375263-2*0.05400321)
0.1829+1.96*Q1HstandardError
0.1829-1.96*Q1HstandardError
#Second clogit model, this one using EventsNum as ordinal predictor
Q1evntClogit <- vglm(Symptoms ~ EventsNum, family = cumulative(parallel = TRUE), data= RAdata)
summary(Q1evntClogit)
lrtest(Q1evntClogit, Q1trashClogit)
#AIC values for each model
AIC(Q1blogit)
AIC(Q1Clogit)
AIC(Q1evntClogit)

#Question 2
pathQ2 <- url("https://ajmolstad.github.io/docs/Arthritis_SideEffects_Apr23.RDS")
SEdata <- readRDS(pathQ2)
#Independent GEE model
indGEE <- gee(side.effects ~ severity, id=patient, family=binomial, corstr = "independence", data=SEdata)
summary(indGEE)
#Because we looking at relationships in 
#observations over time, we will use auto-regressive
atrgGEE <- gee(side.effects ~ severity, id=patient, family=binomial, corstr = "AR-M", data=SEdata)
summary(atrgGEE)
#GLMM model
Q2glmm <- glmer(side.effects ~ severity + (1|patient), family=binomial, data=SEdata)
#Fitted values from GLMM
observed <- tapply(SEdata$side.effects, SEdata$patient, mean)
predVal <- fitted(Q2glmm) 
predSingle <- tapply(predVal, SEdata$patient, function(x){x[1]})
#Plotting
plot(observed, predSingle, pch=20, xlim=c(0,1), ylim=c(0,1))
#LRT
Q2trashGLMM <- glmer(side.effects ~ (1|patient), family=binomial, data=SEdata)
anova(Q2trashGLMM, Q2glmm, test="Chisq")

#Question 3
pathQ3 <- url("https://ajmolstad.github.io/docs/Q3logLinear.RDS")
ACDdata <- readRDS(pathQ3)
Q3logLinModel <- glm(Freq ~ A + C + D + A:C + A:D + C:D, family=poisson(link=log), data=ACDdata)
summary(Q3logLinModel)
#Other models needed for significance testing
Q3indModel <- glm(Freq ~ A + C + D, family=poisson(link=log), data=ACDdata)
Q3satModel <- glm(Freq ~ A + C + D + A:C + A:D + C:D + A:C:D, family=poisson(link=log), data=ACDdata)
#Significance tests
anova(Q3logLinModel, Q3satModel, test="Chisq")
anova(Q3indModel, Q3satModel, test="Chisq")
#Constructing standard error for CI
vcov(Q3logLinModel)
