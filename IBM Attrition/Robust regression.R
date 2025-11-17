library(robustbase)
library(np)
library(readxl)
library(fastDummies)
library(MASS)
library(robust)

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

data = dummy_cols(data, select_columns = c("BusinessTravel", "Department", "EducationField", "JobRole",
                                           "MaritalStatus"), 
                  remove_selected_columns = TRUE,
                  remove_first_dummy = TRUE)

sigmoid = function(x){
  return(1/(1+exp(-x)))
}

################################################

data = data[,c(1,2,3,5,7,9,11,12,15,17,18,19,20,21,22,23,24,25,27,34,40,42)]
# data has the selected covariates now

robfit_hat = glmrob(Attrition ~ . , family = binomial(link=logit), data = data,
                method = "Mqle", weights.on.x = "hat")
summary(robfit_hat)

data = data[,-13]
robfit_hat = glmrob(Attrition ~ . , family = binomial(link=logit), data = data,
                    method = "Mqle", weights.on.x = "hat")
summary(robfit_hat)

data = data[,-15]
robfit_hat = glmrob(Attrition ~ . , family = binomial(link=logit), data = data,
                    method = "Mqle", weights.on.x = "hat")
summary(robfit_hat)

output = numeric(dim(data)[1])
sequences = cut(1:dim(data)[1], breaks = 10, labels = FALSE, include.lowest = TRUE)

h_const = seq(0.2, 1.4, by = 0.2)
err = numeric(length(h_const))

for(j in 1:length(h_const)){
  
  for(i in 1:length(unique(sequences))){
    indexes = which(sequences == i)
    fit = glmrob(Attrition ~ . , family = binomial(link=logit), method = "Mqle",
                 weights.on.x = "hat", data = data[-indexes,],
                 control = glmrobMqle.control(tcc = h_const[j], acc = 1e-4))
    pred = predict(fit, data[indexes,-1])
    output[indexes] = ifelse(sigmoid(pred) > 0.5, 1, 0)
  }
  
  conf_matrix = table(Actual = data[,1], Predicted = output)
  err[j] = (conf_matrix[1,2] + conf_matrix[2,1])/sum(conf_matrix)
  
}

h_opt = h_const[which.min(err)]
err_opt = min(err)

opt_fit = glmrob(Attrition ~ . , family = binomial(link=logit), method = "Mqle",
                 weights.on.x = "hat", data = data, 
                 control = glmrobMqle.control(tcc = h_opt, acc = 1e-16))
summary(opt_fit)

opt_cov = as.matrix(opt_fit$cov)[-1,-1]
heatmap(opt_cov)
weights = opt_fit$w.r
hist(weights, breaks = seq(0, 1, by = 0.020))

prob_pred = sigmoid(predict(opt_fit, data))

hist(prob_pred, breaks = seq(0, 1, by = 0.025), col = "green")
hist(prob_pred[which(data$Attrition == 1)], breaks = seq(0, 1, by = 0.025), col="red", add = TRUE)

risk_index = numeric(length(prob_pred))
risk_index[which(prob_pred >= 0.2)] = 1 # at risk
risk_index[which(prob_pred >= 0.6)] = 2 # warning !!! (most likely already resigned)

table = table(data$Attrition, risk_index)
barplot(table, col = c("green","red"), beside = T, border = T, main = "Predicted Risk Level vs Real Attrition",
        legend.text = c("Attrition No", "Attrition Yes"),
        names.arg = c("Seems ok..(<0.2)", "At Risk (0.2,0.6)",
                      "Warning !! (>0.6)"))
