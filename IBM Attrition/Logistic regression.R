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

fit = glm(Attrition ~ . , family = binomial(link=logit), data = data)
summary(fit)

pearson_res = residuals(fit, type = 'pearson')
deviance_res = residuals(fit, type = 'deviance')
working_res = residuals(fit, type = 'working') # == fit$residuals
response_res = residuals(fit, type = 'response')
pearson_std_res = rstandard(fit, type = 'pearson')
deviance_std_res = rstandard(fit, type = 'deviance')

hist(pearson_res)
hist(deviance_res)
hist(working_res)
hist(response_res)
hist(pearson_std_res)
hist(deviance_std_res)

shapiro.test(scale(pearson_res))
shapiro.test(scale(deviance_res))
shapiro.test(scale(working_res))
shapiro.test(scale(response_res))
shapiro.test(scale(pearson_std_res))
shapiro.test(scale(deviance_std_res))


output = numeric(dim(data)[1])
sequences = cut(1:dim(data)[1], breaks = 10, labels = FALSE, include.lowest = TRUE)

for(i in 1:length(unique(sequences))){
  indexes = which(sequences == i)
  fit = glm(Attrition ~ ., family = binomial(link=logit), data = data[-indexes,])
  pred = predict(fit, data[indexes,-1])
  output[indexes] = ifelse(sigmoid(pred) > 0.5, 1, 0)
}

conf_matrix = table(Actual = data[,1], Predicted = output)
conf_matrix
err = max(conf_matrix[1,2], conf_matrix[2,1])
err

prob_pred = sigmoid(predict(fit, data))

hist(prob_pred, breaks = seq(0, 1, by = 0.025), col = "green")
hist(prob_pred[which(data_hat$Attrition == 1)], breaks = seq(0, 1, by = 0.025), col="red", add = TRUE)
