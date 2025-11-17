library(survival)
library(survminer)
library(dplyr) 
library(ggplot2)
library(knitr)
library(broom)
library(tidyr)
library(readxl)
library(fastDummies)

data = read_excel("attrition_dataset.xlsx")
data = subset(data, select = -c(EmployeeCount, EmployeeNumber, Over18, StandardHours, DailyRate,
                                HourlyRate, MonthlyRate))
data = data[-which(data$YearsAtCompany == 0),]
data = data.frame(data)

data$Attrition = ifelse(data$Attrition == "Yes", 1, 0)
data$Gender = ifelse(data$Gender == "Male", 1, 0)
data$OverTime = ifelse(data$OverTime == "Yes", 1, 0)

temp = data[,1]
data[,1] = data[,2]
data[,2] = temp
rm(temp)
colnames(data)[1] <- "Attrition"
colnames(data)[2] <- "Age"

colnames(data)[colnames(data) == "Department_Research & Development"] <- "Department_Research_Development"
colnames(data)[colnames(data) == "JobRole_Laboratory Technician"] <- "JobRole_Laboratory_Technician"
colnames(data)[colnames(data) == "JobRole_Sales Representative"] <- "JobRole_Sales_Representative"

################################################

data$Attrition <- factor(data$Attrition, labels = (c('Censor', 'Event')))
id <- factor(seq(1, dim(data)[1], 1))
data <- data.frame(data, id)

ggplot(data = data[1:15,], aes(x = id,y = YearsAtCompany)) + geom_bar(stat = 'identity', width = 0.2) + 
  geom_point(aes(color = Attrition, shape = Attrition), size = 6) + coord_flip()

fit <- survfit(Surv(YearsAtCompany, Attrition =='Event') ~ 1, data = data)
summary(fit)

library(knitr)
library(broom)
kable(head(tidy(fit),20))

surv_median(fit)

plot(fit, conf.int = T, xlab='Time [years]', ylab = 'Survival Probability', col='red',
     main="Kaplan-Meier Curve for Attrition")

ggsurvplot(fit,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           break.time.by=10,
           title = "Kaplan-Meier Curve for Attrition",
           data = data)

H <- fit$cumhaz

ggsurvplot(fit,
           risk.table = TRUE,
           ggtheme = theme_bw(),
           fun='cumhaz',
           break.time.by=10,
           title="Cumulative Hazard Curve for Attrition",
           data = data)

################################################

fit.gender <- survfit(Surv(YearsAtCompany, Attrition=='Event') ~ Gender, data=data)

ggsurvplot(fit.gender, conf.int = F,
           risk.table = TRUE,
           risk.table.col = "strata", 
           ggtheme = theme_bw(),
           break.time.by=10,
           legend.labs=c("Female","Male"), legend.title="Gender class",  
           palette=c("violet","blue"), 
           title="Kaplan-Meier Curves by gender class for Attrition",
           data = data)

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ Gender, data = data)
log_rank_test
# p = 0.3 the difference is not significant in survival

hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[2]/log_rank_test$exp[2])
hazard_ratio
# 0.87 < 1
# the risk of attrition in female is 0.87 times the risk in male: being female is a protective factor.

################################################

data$Department <- factor(data$Department)

fit.department <- survfit(Surv(YearsAtCompany, Attrition == 'Event') ~ Department, data = data)

ggsurvplot(fit.department, data = data,
           risk.table = TRUE,
           risk.table.col = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_bw(),
           break.time.by=10,
           legend.labs=c("HR","RD","Sales"), legend.title="Department class",  
           palette=c("darkblue","red","forestgreen"), 
           title = "Kaplan-Meier Curves by department class for Attrition")

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ Department, data = data)
log_rank_test
# p = 0.04 so the difference is significant

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ Department,
                          data = data[c(which(data$Department == "Human Resources"),which(data$Department == "Sales")),])
log_rank_test
# p = 1, Sales and HR are the same

data_depart = data
data_depart$Department = ifelse(data$Department == "Sales", "Sales or HR", ifelse(data$Department == "Human Resources", "Sales or HR", "RD"))

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ Department, data = data_depart)
log_rank_test
# p = 0.01, the two groups are different

hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[2]/log_rank_test$exp[2])
hazard_ratio
# 0.7 < 1
# the risk of attrition in Department RD is 0.7 times the risk in the others departments:
# being part of RD is a protective factor

################################################

data$JobLevel <- factor(data$JobLevel)

fit.Joblevel <- survfit(Surv(YearsAtCompany, Attrition == 'Event') ~ JobLevel, data = data)

ggsurvplot(fit.Joblevel, data = data,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           break.time.by=10,
           legend.labs=c("1","2","3","4","5"), legend.title="JobLevel class",  
           palette=c("darkblue","cyan3","green","purple","blue"), 
           title="Kaplan-Meier Curves by JobLevel class for Attrition")

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ JobLevel, data=data)
log_rank_test
# p < 2e-16 so the difference is very significant

# 2 vs 3
log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ JobLevel, 
                          data = data[c(which(data$JobLevel == 3),which(data$JobLevel == 2)),])
log_rank_test
# p = 1, no difference

# 4 vs 5
log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ JobLevel, 
                          data = data[c(which(data$JobLevel == 4),which(data$JobLevel == 5)),])
log_rank_test
# p = 0.9, no substancial difference

data_joblevel = data
data_joblevel$JobLevel = ifelse(data$JobLevel == 1, "Group 1", ifelse(data$JobLevel == 2, "Group 2/3", ifelse(data$JobLevel == 3, "Group 2/3", "Group 4/5"))) 

# 1 vs 2/3 vs 4/5
log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ JobLevel, data = data_joblevel)
log_rank_test
# p < 2e-16, three macro-groups

# 1 vs 2/3
hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[2]/log_rank_test$exp[2])
hazard_ratio
# 3.73 >> 1
# The risk of JL1 is 3.73 times the risk of JL2/3

# 1 vs 4/5
hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[3]/log_rank_test$exp[3])
hazard_ratio
# 11 >> 1
# The risk of JL1 is 11 times the risk of JL4/5

# 2/3 vs 4/5
hazard_ratio <- (log_rank_test$obs[2]/log_rank_test$exp[2])/(log_rank_test$obs[3]/log_rank_test$exp[3])
hazard_ratio
# 2.94 > 1
# The risk of JL2/3 is 2.94 times the risk of JL4/5

fit.groups.Joblevel <- survfit(Surv(YearsAtCompany, Attrition == 'Event') ~ JobLevel, data = data_joblevel)

ggsurvplot(fit.groups.Joblevel, data = data_joblevel,
           risk.table = TRUE, 
           risk.table.col = "strata",
           ggtheme = theme_bw(),
           break.time.by=10,
           legend.labs=c("1","2/3","4/5"), legend.title="JobLevel class",  
           palette=c("darkblue","cyan3","green"), 
           title="Kaplan-Meier Curves by JobLevel class for Attrition")


################################################

data$BusinessTravel <- factor(data$BusinessTravel)
data$Education <- factor(data$Education)
data$EducationField <- factor(data$EducationField)
data$Gender <- factor(data$Gender)
data$EnvironmentSatisfaction <- factor(data$EnvironmentSatisfaction)
data$JobInvolvement <- factor(data$JobInvolvement)
data$JobRole <- factor(data$JobRole)
data$JobSatisfaction <- factor(data$JobSatisfaction)
data$MaritalStatus <- factor(data$MaritalStatus)
data$OverTime <- factor(data$OverTime)
data$PerformanceRating <- factor(data$PerformanceRating)
data$RelationshipSatisfaction <- factor(data$RelationshipSatisfaction)
data$StockOptionLevel <- factor(data$StockOptionLevel)
data$WorkLifeBalance <- factor(data$WorkLifeBalance)

glimpse(data)

# Assumption : covariates time independent, so we removed YearsInCurrentRole, YearsSinceLastPromotion, 
#              YearsSinceLastPromotion, YearsWithCurrManager, TotalWorkingYears 

mod.cox <- coxph(Surv(YearsAtCompany, Attrition == 'Event') ~ 
                   Age + BusinessTravel + DistanceFromHome +
                   EnvironmentSatisfaction + Gender + JobInvolvement + I(JobLevel == 1) +
                   I(JobLevel == 2)  + I(JobLevel == 4)  +
                   I(JobRole == "Research Scientist") + I(JobRole == "Sales Executive") + 
                   I(JobRole == "Sales Representative") + I(JobSatisfaction == 4) + MonthlyIncome + 
                   NumCompaniesWorked + OverTime + RelationshipSatisfaction + 
                   I(StockOptionLevel == 0) + TrainingTimesLastYear + WorkLifeBalance, data = data)
summary(mod.cox)

ggcoxdiagnostics(mod.cox, type = "martingale")

ggcoxdiagnostics(mod.cox, type = "deviance")

ggcoxdiagnostics(mod.cox, type = "schoenfeld")

ggcoxdiagnostics(mod.cox, type = "scaledsch")

test.ph <- cox.zph(mod.cox)
test.ph
