library(fastDummies)
library(MASS)
library(rgl)
library(DepthProc)
library(ggplot2)
library(readxl)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(corrplot)
library(pbapply)
library(ISLR2)

library(mgcv)
library(locfit)
library(kernlab)

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

data = dummy_cols(data, select_columns = c("BusinessTravel", "Department", "EducationField", "JobRole",
                                           "MaritalStatus"), 
                  remove_selected_columns = TRUE,
                  remove_first_dummy = TRUE)

colnames(data)[colnames(data) == "JobRole_Laboratory Technician"] <- "JobRole_Laboratory_Technician"
colnames(data)[colnames(data) == "Department_Research & Development"] <- "Department_Research_Development"
colnames(data)[colnames(data) == "JobRole_Sales Representative"] <- "JobRole_Sales_Representative"

sigmoid = function(x){
  return(1/(1+exp(-x)))
}

################################################
##  "Attrition"                        "Age"                             
##  "DistanceFromHome"                 "EnvironmentSatisfaction"         
##  "JobInvolvement"                   "JobSatisfaction"                 
##  "NumCompaniesWorked"               "OverTime"                        
##  "RelationshipSatisfaction"         "TotalWorkingYears"               
##  "TrainingTimesLastYear"            "WorkLifeBalance"                 
##  "YearsAtCompany"                   "YearsInCurrentRole"              
##  "YearsSinceLastPromotion"          "YearsWithCurrManager"            
##  "BusinessTravel_Travel_Frequently" "BusinessTravel_Travel_Rarely"    
##  "Department_Sales"                 "JobRole_Laboratory_Technician"   
##  "JobRole_Sales_Representative"     "MaritalStatus_Single"

hist(data$Age)
betweenAge = aggregate(Attrition ~ Age, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenAge$Attrition, names.arg = betweenAge$Age)

betweenDFH = aggregate(Attrition ~ DistanceFromHome, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenDFH$Attrition, names.arg = betweenDFH$DistanceFromHome)

betweenES = aggregate(Attrition ~ EnvironmentSatisfaction, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenES$Attrition, names.arg = betweenES$EnvironmentSatisfaction)

betweenJI = aggregate(Attrition ~ JobInvolvement, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenJI$Attrition, names.arg = betweenJI$JobInvolvement)

betweenJS = aggregate(Attrition ~ JobSatisfaction, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenJS$Attrition, names.arg = betweenJS$JobSatisfaction)

betweenNCW = aggregate(Attrition ~ NumCompaniesWorked, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenNCW$Attrition, names.arg = betweenNCW$NumCompaniesWorked)

betweenOT = aggregate(Attrition ~ OverTime, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenOT$Attrition, names.arg = betweenOT$OverTime)

betweenRS = aggregate(Attrition ~ RelationshipSatisfaction, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenRS$Attrition, names.arg = betweenRS$RelationshipSatisfaction)

betweenTWY = aggregate(Attrition ~ TotalWorkingYears, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenTWY$Attrition, names.arg = betweenTWY$TotalWorkingYears)

betweenTTLY = aggregate(Attrition ~ TrainingTimesLastYear, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenTTLY$Attrition, names.arg = betweenTTLY$TrainingTimesLastYear)

betweenWLB = aggregate(Attrition ~ WorkLifeBalance, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenWLB$Attrition, names.arg = betweenWLB$WorkLifeBalance)

betweenYAC = aggregate(Attrition ~ YearsAtCompany, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenYAC$Attrition, names.arg = betweenYAC$YearsAtCompany)

betweenYCR = aggregate(Attrition ~ YearsInCurrentRole, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenYCR$Attrition, names.arg = betweenYCR$YearsInCurrentRole)

betweenYLP = aggregate(Attrition ~ YearsSinceLastPromotion, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenYLP$Attrition, names.arg = betweenYLP$YearsSinceLastPromotion)

betweenYWM = aggregate(Attrition ~ YearsWithCurrManager, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenYWM$Attrition, names.arg = betweenYWM$YearsWithCurrManager)

betweenBTF = aggregate(Attrition ~ BusinessTravel_Travel_Frequently, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenBTF$Attrition, names.arg = betweenBTF$BusinessTravel_Travel_Frequently)

betweenBTR = aggregate(Attrition ~ BusinessTravel_Travel_Rarely, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenBTR$Attrition, names.arg = betweenBTR$BusinessTravel_Travel_Rarely)

betweenDRD = aggregate(Attrition ~ Department_Research_Development, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenDRD$Attrition, names.arg = betweenDRD$Department_Research_Development)

betweenDS = aggregate(Attrition ~ Department_Sales, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenDS$Attrition, names.arg = betweenDS$Department_Sales)

betweenJRT = aggregate(Attrition ~ JobRole_Laboratory_Technician, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenJRT$Attrition, names.arg = betweenJRT$JobRole_Laboratory_Technician)

betweenJRS = aggregate(Attrition ~ JobRole_Sales_Representative, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenJRS$Attrition, names.arg = betweenJRS$JobRole_Sales_Representative)

betweenMSS = aggregate(Attrition ~ MaritalStatus_Single, data = data, FUN = function(x) mean(x == 1) * 100)
barplot(betweenMSS$Attrition, names.arg = betweenMSS$MaritalStatus_Single)

gen_add_model = gam(Attrition ~ s(Age, bs = 'cr') + s(DistanceFromHome, bs = 'cr') +
                      I(EnvironmentSatisfaction == 1) + JobInvolvement + JobSatisfaction + 
                      s(NumCompaniesWorked, bs = 'cr') + OverTime + I(RelationshipSatisfaction == 1) + 
                      s(TotalWorkingYears, bs = 'cr') + I(WorkLifeBalance == 1) + 
                      s(YearsAtCompany, bs = 'cr') + s(YearsSinceLastPromotion, bs = 'cr') + 
                      s(YearsWithCurrManager, bs = 'cr') + BusinessTravel_Travel_Frequently +
                      BusinessTravel_Travel_Rarely + Department_Research_Development + 
                      JobRole_Laboratory_Technician + MaritalStatus_Single , 
                      family = binomial(link=logit), data = data)

summary(gen_add_model)
#plot(gen_add_model)
gam.check(gen_add_model) # for checking the number of basis for splines

output = numeric(dim(data)[1])
sequences = cut(1:dim(data)[1], breaks = 10, labels = FALSE, include.lowest = TRUE)

for(i in 1:length(unique(sequences))){
  indexes = which(sequences == i)
  gen_add = gam(Attrition ~ s(Age, bs = 'cr') + s(DistanceFromHome, bs = 'cr') +
                  I(EnvironmentSatisfaction == 1) + JobInvolvement + JobSatisfaction + 
                  s(NumCompaniesWorked, bs = 'cr') + OverTime + I(RelationshipSatisfaction == 1) + 
                  s(TotalWorkingYears, bs = 'cr') + I(WorkLifeBalance == 1) + 
                  s(YearsAtCompany, bs = 'cr') + s(YearsSinceLastPromotion, bs = 'cr') + 
                  s(YearsWithCurrManager, bs = 'cr') + BusinessTravel_Travel_Frequently +
                  BusinessTravel_Travel_Rarely + Department_Research_Development + 
                  JobRole_Laboratory_Technician + MaritalStatus_Single ,
                  family = binomial(link=logit), data = data[-indexes,])
  pred = predict(gen_add, data[indexes,-1])
  output[indexes] = ifelse(sigmoid(pred) > 0.5, 1, 0)
}

conf_matrix = table(Actual = data[,1], Predicted = output)
conf_matrix
conf_matrix[1,2] + conf_matrix[2,1]
err = (conf_matrix[1,2] + conf_matrix[2,1])/sum(conf_matrix)
err

prob_pred = sigmoid(predict(gen_add_model, data))

hist(prob_pred, breaks = seq(0, 1, by = 0.025), col = "green")
hist(prob_pred[which(data$Attrition == 1)], breaks = seq(0, 1, by = 0.025), col="red", add = TRUE)

