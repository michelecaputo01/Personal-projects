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
library(glmperm)
library(lmtest)

library(pbapply)
library(parallel)

set.seed(2023)
B = 1000

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

data = data[uncontaminated_data,]

colnames(data)[colnames(data) == "Department_Research & Development"] <- "Department_Research_Development"
colnames(data)[colnames(data) == "JobRole_Laboratory Technician"] <- "JobRole_Laboratory_Technician"
colnames(data)[colnames(data) == "JobRole_Sales Representative"] <- "JobRole_Sales_Representative"

data_dumm = dummy_cols(data, select_columns = c("BusinessTravel", "Department", "EducationField", "JobRole",
                                           "MaritalStatus"), 
                  remove_selected_columns = TRUE,
                  remove_first_dummy = TRUE)

################################################

fit = glm(Attrition ~ ., family = binomial(link=logit), data = data)
summary(fit)

sigmoid = function(x){
  return(1/(1+exp(-x)))
}

variable_selection_pval_perm = function(d, index, B, print = TRUE){
  n = dim(d)[1]
  
  fit = glm(Attrition ~ ., family = binomial(link=logit), data = d)
  T0 = abs(summary(fit)$coefficients[index,3])
  
  redfit = glm(Attrition ~ . , family = binomial(link=logit), data = d[,-index])
  res = redfit$residuals
  
  stat = numeric(B)
  newdata = d
  for(i in 1:B){
    perm = sample(n)
    linear_perm = redfit$linear.predictors + res[perm]
    newdata$Attrition = ifelse(sigmoid(linear_perm) >= 0.5, 1, 0)
    stat[i] = abs(summary(glm(Attrition ~ . , family = binomial(link=logit), data = newdata))$coefficients[index,3])
  }
  
  rm(newdata)
  
  if(print == TRUE){
    hist(stat, xlim = c(min(c(stat,T0)), max(stat,T0)))
    abline(v = T0, col = 'green')
  }
  
  result = sum(stat >= T0)/B
  print(result)
  return(result)
}

variable_selection_pval_perm(data, index = 2, B = 1000) 

variable_selection_pval_perm(data, index = 15, B = 1000)

variable_selection_pval_perm(data, index = 17, B = 1000)

system.time({variable_selection_pval_perm(data, index = 15, B = 1000)})

################################################

variable_selection_pval_perm_parallelized = function(d, index, B, print = TRUE) {
  n = dim(d)[1]
  
  fit = glm(Attrition ~ ., family = binomial(link = logit), data = d)
  T0 = abs(summary(fit)$coefficients[index, 3])
  
  redfit = glm(Attrition ~ . , family = binomial(link = logit), data = d[, -index])
  res = redfit$residuals
  
  stat = numeric(B)
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterExport(cl = cl, varlist = c('sigmoid'), envir = .GlobalEnv)
  clusterExport(cl = cl, varlist = c('n', 'redfit', 'res', 'd', 'index'), envir = environment())
  
  stat = unlist(parSapply(cl, 1:B, function(i) {
    perm = sample(n)
    linear_perm = redfit$linear.predictors + res[perm]
    newdata = d
    newdata$Attrition = as.integer(sigmoid(linear_perm) >= 0.5)
    abs(summary(glm(Attrition ~ . , family = binomial(link = logit), data = newdata))$coefficients[index, 3])
  }))
  
  stopCluster(cl)
  
  if(print){
    hist(stat, xlim = c(min(c(stat, T0)), max(stat, T0)))
    abline(v = T0, col = 'green', lwd = 3)
  }
  
  result = sum(stat >= T0)/B
  print(result)
  return(result)
  
}

system.time({variable_selection_pval_perm_parallelized(data, index = 16, B = 1000)}) # about 18 seconds after parallelizing

###############################################################

prr_test_parallelized = function(d, index, B = 1000, f = 1, print = FALSE) {
  
  n = dim(d)[1]
  data1 = d[,-c(1,index)]
  data2 = d[,-index]
  
  model1 = lm(d[,index] ~ ., data = data1)
  r = model1$residuals
  
  model2 = glm(Attrition ~ r + ., family = binomial(link=logit), data = data2)
  p0 = summary(model2)$coefficients[2,4]
  
  stat = numeric(B)
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterExport(cl = cl, varlist = c('n', 'r', 'data2'), envir = environment())
  
  stat = unlist(parSapply(cl, 1:B, function(i) {
    perm = sample(n)
    rstar = r[perm]
    summary(glm(Attrition ~ rstar + ., family = binomial(link=logit), data = data2))$coefficients[2,4]
  }))
  
  stopCluster(cl)
  
  if(print){
    title = paste("P_value of ",colnames(d)[index],",", index, "of", dim(d)[2], "covariates.")
    hist(stat, xlim = c(0,1), main = title)
    abline(v = p0, col = 'green', lwd = 3)
  }
  
  result = sum(stat <= f * p0)/B
  
  se = sqrt( f*p0*(1-f*p0) / B)
  
  if(print){ 
    print(result)
    }
  
  return(result)
  
}

system.time({prr_test_parallelized(data, index = 2, B = 1000, print = TRUE)}) # 20 seconds

#prr_test_parallelized si riferisce al paper, gli altri funzionano ma senza valenza teorica

perm_var_selection = function(d, alpha = 0.05, nsim = 300) {
  
  pvalues = rep(1,dim(d)[2])
  pvalues[1] = 0
  remake = 1
  
  while(remake == 1){
    
    for(i in 2:length(pvalues)){
      pvalues[i] = prr_test_parallelized(d, index = i, B = nsim, print = TRUE)
    }
    
    if(sum(pvalues > alpha) == 0){
      remake = 0
    }
    
    if(remake == 1){
      
      if(sum(pvalues > 0.5) >= 1){
        d = d[,-which(pvalues > 0.5)]
        pvalues = pvalues[-which(pvalues > 0.5)]
      }
      
      else{
        if(sum(pvalues > alpha) >= 1){
          d = d[,-which.max(pvalues)]
          pvalues = pvalues[-which.max(pvalues)]
        }
      }
    }
  }
  
  print(colnames(d))
  return(d)
}

selected_data = perm_var_selection(data, alpha = 0.10, nsim = 300)

selected_data = perm_var_selection(selected_data, alpha = 0.05, nsim = 1000)

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

selected_data = data[,c(1,2,3,5,7,9,11,12,15,17,18,19,20,21,22,23,24,25,27,34,40,42)]

glm_red = glm(Attrition ~ ., family = binomial(link=logit), data = selected_data)
summary(glm_red)

output = numeric(dim(selected_data)[1])
sequences = cut(1:dim(selected_data)[1], breaks = 10, labels = FALSE, include.lowest = TRUE)

for(i in 1:length(unique(sequences))){
  indexes = which(sequences == i)
  fit = glm(Attrition ~ ., family = binomial(link=logit), data = selected_data[-indexes,])
  pred = predict(fit, selected_data[indexes,-1])
  output[indexes] = ifelse(sigmoid(pred) > 0.5, 1, 0)
}

conf_matrix = table(Actual = selected_data[,1], Predicted = output)
conf_matrix
conf_matrix[1,2] + conf_matrix[2,1]
(conf_matrix[1,2] + conf_matrix[2,1])/sum(conf_matrix)