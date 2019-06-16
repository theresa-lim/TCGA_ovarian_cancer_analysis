require(ggplot2)

#PAK4 PLOTS
qplot(PAK4.z.score, clinical_stage, data = patients_all) 
qplot(PAK4.expression_log2, clinical_stage, data = patients_all)

qplot(PAK4.z.score, days_lived, data = patients_alive) 
qplot(PAK4.expression_log2, days_lived, data = patients_alive)

qplot(days_lived, data = patients_alive) 
qplot(days_to_death, data = patients_dead) 

qplot(age_at_initial_pathologic_diagnosis, PAK4.z.score, data = patients_numeric) 
qplot(PAK4.expression_log2, days_to_death, data = patients_dead)

qplot(PAK4.z.score, vital_status, data = patients_all) 
qplot(PAK4.expression_log2, vital_status, data = patients_all)

#YEAR vs SURVIVAL
qplot(year_of_initial_pathologic_diagnosis, vital_status, data = patients_all)
qplot(age_at_initial_pathologic_diagnosis, days_to_death, data = patients_dead)
qplot(age_at_initial_pathologic_diagnosis, days_lived, data = patients_alive)

#ALIVE
patients_alive <- patients_all[which(patients_all$vital_status == "alive"), ]
ggplot(patients_alive, aes(ROR1.expression_log2)) +
  geom_histogram(binwidth = 0.1)

#DEAD
patients_dead <- patients_all[which(patients_all$vital_status == "dead"), ]
ggplot(patients_dead, aes(ROR1.expression_log2)) +
  geom_histogram(binwidth = 0.1)

qplot(clinical_stage, data = patients_numeric) 
qplot(age_at_initial_pathologic_diagnosis, days_to_death, data= patients_numeric)
qplot(clinical_stage, days_lived, data= patients_alive)

#####HISTOGRAMS####
qplot(BRCA1.expression_log2, data = patients_all)
qplot(BRCA1.z.score, data = patients_all)

qplot(BRCA1.expression_log2, data = patients_dead)
qplot(BRCA1.z.score, data = patients_dead)

####Exploratory regressions, etc####

model <- lm(days_to_death ~ clinical_stage + BRCA1.expression_log2 + BRCA1.z.score + BRCA2.expression_log2 + BRCA2.z.score + PTEN.expression_log2 + PTEN.z.score + CLASP1.expression_log2 + CLASP1.z.score + PAK1.expression_log2 + PAK1.z.score +PAK2.expression_log2 + PAK2.z.score + PAK4.expression_log2 + PAK4.z.score + ROR1.expression_log2 + ROR1.z.score + TP53.expression_log2 + TP53.z.score, data = patients_dead)
summary(model)
vif(model)
library(car)

anova(model2)
par(mfrow=c(2,2))
plot(model2)

library (MASS)

subdat <- patients_numeric[1:301,c("age_at_initial_pathologic_diagnosis", "days_to_death", "clinical_stage", "BRCA1.expression_log2", "BRCA1.z.score", "BRCA2.expression_log2", "BRCA2.z.score", "PTEN.expression_log2", "PTEN.z.score", "CLASP1.expression_log2", "CLASP1.z.score", "PAK4.expression_log2", "PAK4.z.score", "ROR1.expression_log2", "ROR1.z.score", "TP53.expression_log2", "TP53.z.score")]
model3 <- stepAIC(model, direction = "backward")
summary(model3)

require(corrplot)
mcor <- cor(subdat, use = "complete.obs")
mcor
corrplot(mcor, method="shade", shade.col=NA, tl.col="black",tl.cex=0.5)

##### OLS Regression #####

OLS.z.score <- lm(days_to_death ~ PAK4.z.score + clinical_stage + age_at_initial_pathologic_diagnosis, data = patients_numeric)
summary(model5)

OLS.expression_log <- lm(days_to_death ~ PAK4.expression_log2 + clinical_stage + age_at_initial_pathologic_diagnosis, data = patients_numeric)
summary(OLS.expression_log)

#####survival analysis: Kaplan Meier + Cox#####
install.packages("gof")
require(gof)
require(survival)
require(survminer)

surv_object <- Surv(time = patients_numeric$days_to_death, event = patients_numeric$vital_status)
surv_model_ALL <- lm(surv_object ~ clinical_stage + BRCA1.expression_log2 + BRCA1.z.score + BRCA2.expression_log2 + BRCA2.z.score + PTEN.expression_log2 + PTEN.z.score + CLASP1.expression_log2 + CLASP1.z.score + PAK1.expression_log2 + PAK1.z.score +PAK2.expression_log2 + PAK2.z.score + PAK4.expression_log2 + PAK4.z.score + ROR1.expression_log2 + ROR1.z.score + TP53.expression_log2 + TP53.z.score, data = patients_numeric)
surv_model2_3VAR <- lm(surv_object ~ clinical_stage + PAK4.expression_log2 + age_at_initial_pathologic_diagnosis, data = patients_numeric)


cox_model <- coxph(surv_object ~ PAK4.expression_log2+clinical_stage+age_at_initial_pathologic_diagnosis, data = patients_numeric)
ggforest(cox_model, data = patients_numeric)
summary(cox_model)
AIC(cox_model)
visreg(cox_model, "PAK4.expression_log2")
cumres(cox_model)

#####Survival analysis: Weibull#####
install.packages("SurvRegCensCov")
install.packages("visreg")
require(SurvRegCensCov)
require(MASS)
require(visreg)

weibull <- survreg(surv_object ~ PAK4.expression_log2 + clinical_stage + age_at_initial_pathologic_diagnosis, data = patients_numeric, dist = "weibull")
summary(weibull)
AIC(weibull, cox_model, OLS.expression_log)
visreg(weibull, "age_at_initial_pathologic_diagnosis")

WeibullReg(surv_object ~ PAK4.expression_log2 + clinical_stage + age_at_initial_pathologic_diagnosis, data = patients_numeric)
summary(weibull2)

WeibullDiag(surv_object ~ PAK4dummy_var, data = patients_numeric)

#####Survival analysis: logistic#####
qplot(days_to_death, data = patients_numeric, binwidth = 100)

logistic <- glm(lived_longer_than_3000 ~ PAK4.expression_log2 + clinical_stage + age_at_initial_pathologic_diagnosis, data = patients_numeric, family = binomial)
summary(logistic)
glm.probs <- predict(logistic,type = "response")

qplot(age_at_initial_pathologic_diagnosis, data = patients_numeric,  binwidth = 1)

#####Testing logistic assumptions#####
install.packages("tidyverse")
install.packages("broom")
install.packages("magrittr")
require(tidyverse)
require(broom)

probabilities <- predict(glm.fit, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
head(predicted.classes)

mydata <- patients_numeric[c("age_at_initial_pathologic_diagnosis", "clinical_stage", "PAK4.expression_log2")]
predictors <- colnames(mydata)

mydata <- mydata %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() +
  facet_wrap(~predictors, scales = "free_y")

#####determining where to break PAK4 dummy variable #####
qplot(PAK4.expression_log2, data = patients_numeric)
median(patients_numeric$PAK4.expression_log2)

qplot(days_to_death, data = patients_numeric)
