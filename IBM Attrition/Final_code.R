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
library(np)
library(robustbase)
library(glmperm)
library(lmtest)
library(pbapply)
library(parallel)
library(ISLR2)
library(survival)
library(survminer)

sigmoid = function(x){
  return(1/(1+exp(-x)))
}

####  Prepare the dataset ------------------------------------------------------
ds = read_excel("attrition_dataset.xlsx")

# Remove variables that are meaningless or useless for our analysis
ds = subset(ds, select = -c(EmployeeCount, EmployeeNumber, Over18, StandardHours, 
                                DailyRate, HourlyRate, MonthlyRate))
# Remove data with years at company equal to 0 since they have no sens for our analysis
ds = ds[-which(ds$YearsAtCompany == 0),]
ds = data.frame(ds)

ds$Attrition = ifelse(ds$Attrition == "Yes", 1, 0)
ds$Gender = ifelse(ds$Gender == "Male", 1, 0)
ds$OverTime = ifelse(ds$OverTime == "Yes", 1, 0)

# Attrition at the firs position
temp = ds[,1]
ds[,1] = ds[,2]
ds[,2] = temp
rm(temp)
colnames(ds)[1] <- "Attrition"
colnames(ds)[2] <- "Age"

# Save only numerical variable
numericvar = c()
for(i in 1:dim(ds)[2])
  if(is.numeric(ds[1,i])) numericvar = c(numericvar,i)
datanum = ds[,numericvar]

# Transform categorical variables in dummies
data = dummy_cols(ds, select_columns = c("BusinessTravel", "Department", 
                                           "EducationField", "JobRole", "MaritalStatus"), 
                  remove_selected_columns = TRUE,
                  remove_first_dummy = TRUE)

# Changes on the names of the columns in order to remove spaces
colnames(data)[colnames(data) == "Department_Research & Development"] <- "Department_Research_Development"
colnames(data)[colnames(data) == "JobRole_Laboratory Technician"] <- "JobRole_Laboratory_Technician"
colnames(data)[colnames(data) == "JobRole_Sales Representative"] <- "JobRole_Sales_Representative"

# data has dummies, ds is the pure dataset, datanum has only the numeric variables

####  Explorative analysis -----------------------------------------------------

attach(data)

barplot(table(Attrition), col=c(4,3), border=T, ylim=c(0,1233), cex.axis=0.8, cex.names=0.8, 
        main="Employee Attrition", names.arg=c("No", "Yes"))
#unbalanced
yes = which(data$Attrition == 1)   #221
no = which(data$Attrition == 0)   #1205

boxplot(data$Age ~ data$Attrition, main="Boxplot",xlab="Age",ylab = "Attrition",horizontal=TRUE,col=c("pink","lightblue"))
boxplot(data$DistanceFromHome~data$Attrition,main="Boxplot",xlab="DistanceFromHome",ylab = "Attrition",horizontal=TRUE,col=c("pink","lightblue"))
boxplot(data$MonthlyIncome~data$Attrition,main="Boxplot",xlab="MonthlyIncome",ylab = "Attrition(Yes/No)",horizontal=TRUE,col=c("pink","lightblue"))


# Corralation among variables
corrmatrix = cor(datanum)
heatmap(corrmatrix, Rowv = NA, Colv = NA)
corrplot(corrmatrix, type = "upper")

which(corrmatrix >= 0.75)

# We can see there are some variables highly correlated:
corrmatrix["JobLevel", "MonthlyIncome"]               # 0.950043
corrmatrix["JobLevel", "TotalWorkingYears"]           # 0.780533
corrmatrix["MonthlyIncome", "TotalWorkingYears"]      # 0.769697
corrmatrix["PercentSalaryHike", "PerformanceRating"]  # 0.774875
corrmatrix["YearsAtCompany", "YearsWithCurrManager"]  # 0.759394
corrmatrix["YearsAtCompany", "YearsInCurrentRole"]    # 0.748387

scatterplotMatrix(data.frame(MonthlyIncome, JobLevel))
boxplot(MonthlyIncome ~ JobLevel)
anova1 = aov(MonthlyIncome ~ JobLevel)
summary(anova1)
# As expected monthly income is increasing with an increasing Jjob level

# We look also at the relationship between percent salary hike and performance rating
pairs(datanum[,c(13,14)])
model2 = lm(PerformanceRating ~ cut(PercentSalaryHike, 
                                         breaks = c(-Inf,19.5,max(PercentSalaryHike))))
summary(model2)
plot(PercentSalaryHike, PerformanceRating)
points(PercentSalaryHike, predict(model2),lwd =2, col =" blue")
# As expected the percent salary hike increase with the rate of the performance


# Gender analysis
barplot(table(data$Gender), col=c(10,4), border=T, ylim=c(0,1233), cex.axis=0.8, cex.names=0.8, 
        main="Employee Gender", names.arg=c("Female", "Male"))

boxplot(data$Age ~ data$Gender, border=T, col=c(10,4),cex.axis=0.8, cex.names=0.8, 
        names=c("Female", "Male"), xlab="Gender", ylab= "Age")

boxplot(data$JobSatisfaction ~ data$Gender + data$Attrition, col=c(10,4,3,5), border=T, cex.axis=0.8, cex.names=0.8, 
        names=c("Female, No", "Male, No", "Female, Yes","Male, Yes"),xlab="Gender", ylab= "Job satisfaction")

boxplot(data$MonthlyIncome ~ data$Gender, border=T, col = c(10,4),cex.axis=0.8, cex.names=0.8, 
        names=c("Female", "Male"), xlab="Gender", ylab= "Age")

length(which(data$Attrition== 1 & data$Gender==0))  #83 women
length(which(data$Attrition== 1 & data$Gender==1))  #138 man

# costruisco la tabella
table1<-table(Gender, ds$Department)

# costruisco il grafico della tabella
barplot(table1, col=c(10,4),beside=T, border=T, main="Gender division in departments",
        legend.text=c("Female","Male"),
        names.arg=c("HR", "R&D","Sales"))

# costruisco la tabella
table1<-table(Attrition, ds$Department)

# costruisco il grafico della tabella
barplot(table1, col=c(10,4),beside=T, border=T, main="Attrition division in departments",
        legend.text=c("No","Yes"),
        names.arg=c("HR", "R&D","Sales"))


detach(data)

####  Outlier detection  -----------------------------------------------------

data_for_outliers = ds[,-c(1,3,4,7,9,12,14,17,19)]
data_for_outliers = as.data.frame(lapply(data_for_outliers, scale))

set.seed(2023)

fit_MCD = covMcd(x = data_for_outliers, alpha = .75, nsamp = 500, cor = TRUE)
fit_MCD

#plot(fit_MCD)

robust_cov = fit_MCD$cov
robust_cor = fit_MCD$cor
heatmap(robust_cov, Rowv = NA, Colv = NA)
heatmap(robust_cor, Rowv = NA, Colv = NA)

uncontaminated_data = which(mahalanobis(x = data_for_outliers, center = fit_MCD$raw.center,
                                        cov = fit_MCD$raw.cov) <= qchisq(p = .999, df = dim(data_for_outliers)[2]))
rm(data_for_outliers)

data_clean = data[uncontaminated_data,]

####  Permutation tests --------------------------------------------------------

# Permutational univarite tests using trimmed mean - robust method

perm_trim_test=function(x,y,iter=10^5){
  
  T0=abs(mean(x,trim=0.1)-mean(y,trim=0.1))  # define the test statistic
  T_stat=numeric(iter) 
  x_pooled=c(x,y) 
  n=length(x_pooled)
  n1=length(x)
  
  set.seed(2024)
  for(perm in 1:iter){ # loop for conditional MC
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    
    # Test statistic:
    T_stat[perm] <- abs(mean(x1_perm, trim=0.1) - mean(x2_perm, trim=0.1))
  }
  
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  return(p_val)
}


# 1) Permutational test: X = attrition yes vs Y = attrition no
# Do these populations have the same distribution with respect to some variables?

yes = which(data$Attrition == 1)   
no = which(data$Attrition == 0)   

# Wrt age
age_yes <- data[yes,"Age"]
age_no <- data[no,"Age"]

p.value <- perm_trim_test(age_yes, age_no)
p.value
# 0

# Wrt distance from home
dist_yes <- data[yes,"DistanceFromHome"]
dist_no <- data[no,"DistanceFromHome"]

p.value <- perm_trim_test(dist_yes,dist_no)
p.value
# 0.00018

# Wrt monthly income
monthly_yes <- data[yes,"MonthlyIncome"]
monthly_no <- data[no,"MonthlyIncome"]

p.value <- perm_trim_test(monthly_yes,monthly_no)
p.value
# 0

# Wrt total working years
total_yes <- data[yes,"TotalWorkingYears"]
total_no <- data[no,"TotalWorkingYears"]

p.value <- perm_trim_test(total_yes,total_no)
p.value
# 0

# Wrt percent salary hike
perc_yes <- data[yes,"PercentSalaryHike"]
perc_no <- data[no,"PercentSalaryHike"]

p.value <- perm_trim_test(perc_yes,perc_no)
p.value
# 0.5429


# Plot the distributions of these populations, comparing them
par(mfrow=c(2,5))

hist(age_yes, main='Attrition Yes - Age',xlab='')
hist(total_yes, main='Attrition Yes - Total working years', xlab='')
hist(dist_yes, main='Attrition Yes - Distance from home', xlab='', breaks=6)
hist(monthly_yes, main='Attrition Yes - Monthly income', xlab='')
hist(perc_yes, main='Attrition Yes - Percent salary hike', xlab='', xlim=c(10,25), breaks=8)

hist(age_no, main='Attrition No - Age', xlab='')
hist(total_no, main='Attrition No - Total working years', xlab='')
hist(dist_no, main='Attrition No - Distance from home', xlab='', breaks=6)
hist(monthly_no, main='Attrition No - Monthly income', xlab='')
hist(perc_no ,main='Attrition No - Percent salary hike', xlab='', xlim=c(10,25), breaks=8)
graphics.off()


## 1) Permutational test: X = male vs Y = female
# Do these populations have the same distribution with respect to some varibales?

m <- which(data$Gender==1)
f <- which(data$Gender==0)

# Wrt age
age_male <- data[m,"Age"]
age_female <- data[f,"Age"]

p.value <- perm_trim_test(age_male, age_female)
p.value
# 0.1722

# Wrt distance from home
dist_male <- data[m,"DistanceFromHome"]
dist_female <- data[f,"DistanceFromHome"]

p.value <- perm_trim_test(dist_male, dist_female)
p.value
# 0.84179

# Wrt monthly income
monthly_male <- data[m,"MonthlyIncome"]
monthly_female <- data[f,"MonthlyIncome"]

p.value <- perm_trim_test(monthly_male, monthly_female)
p.value
# 0.16967

# Wrt percent salary hike
perc_male <- data[m,"PercentSalaryHike"]
perc_female <- data[f,"PercentSalaryHike"]

p.value <- perm_trim_test(perc_male, perc_female)
p.value
# 0.59514

# Plot the distributions of these populations, comparing them
par(mfrow=c(2,4))

hist(age_male, main='Male - Age',xlab='')
hist(dist_male, main='Male - Distance from home', xlab='')
hist(monthly_male, main='Male - Monthly income', xlab='')
hist(perc_male, main='Male - Percent salary hike', xlab='')

hist(age_female, main='Female - Age', xlab='')
hist(dist_female, main='Female - Distance from home', xlab='')
hist(monthly_female, main='Female - Monthly income', xlab='')
hist(perc_female ,main='Female - Percent salary hike', xlab='')
graphics.off()



####  Variable selection -------------------------------------------------------

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

#prr_test_parallelized si riferisce al paper

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

####  Logistic regression ------------------------------------------------------

fit_glm = glm(Attrition ~ ., family = binomial(link=logit), data = selected_data)
summary(fit_glm)

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

####  GAMs ---------------------------------------------------------------------

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
prob_pred = numeric(dim(data)[1])
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
  for(j in indexes){
    prob_pred[j] = sigmoid(predict(gen_add, data[j,-1]))
  }
}

for(i in 1:dim(data)[1]){
  output[i] = ifelse(prob_pred[i] > 0.5, 1, 0)
}

conf_matrix = table(Actual = data[,1], Predicted = output)
conf_matrix
conf_matrix[1,2] + conf_matrix[2,1]
err = (conf_matrix[1,2] + conf_matrix[2,1])/sum(conf_matrix)
err

quantile(prob_pred, probs = 0.75)
quantile(prob_pred, probs = 0.95)

hist(prob_pred, breaks = seq(0, 1, by = 0.025), col = "green")
hist(prob_pred[which(data$Attrition == 1)], breaks = seq(0, 1, by = 0.025), col="red", add = TRUE)

risk_index = numeric(length(prob_pred))
risk_index[which(prob_pred >= 0.2)] = 1 # at risk
risk_index[which(prob_pred >= 0.6)] = 2 # warning !!! (most likely already resigned)

table = table(data$Attrition, risk_index)
barplot(table, col = c("green","red"), beside = T, border = T, main = "Predicted Risk Level vs Real Attrition",
        legend.text = c("Attrition No", "Attrition Yes"),
        names.arg = c("Low", "Medium", "High"))

####  Robust logistic regression -----------------------------------------------

robfit_hat = glmrob(Attrition ~ . , family = binomial(link=logit), data = selected_data,
                    method = "Mqle", weights.on.x = "hat")
summary(robfit_hat)

selected_data = selected_data[,-13]
robfit_hat = glmrob(Attrition ~ . , family = binomial(link=logit), data = selected_data,
                    method = "Mqle", weights.on.x = "hat")
summary(robfit_hat)

selected_data = selected_data[,-15]
robfit_hat = glmrob(Attrition ~ . , family = binomial(link=logit), data = selected_data,
                    method = "Mqle", weights.on.x = "hat")
summary(robfit_hat)

output = numeric(dim(selected_data)[1])
sequences = cut(1:dim(selected_data)[1], breaks = 10, labels = FALSE, include.lowest = TRUE)

h_const = seq(0.2, 1.4, by = 0.2)
err = numeric(length(h_const))

for(j in 1:length(h_const)){
  
  for(i in 1:length(unique(sequences))){
    indexes = which(sequences == i)
    fit = glmrob(Attrition ~ . , family = binomial(link=logit), method = "Mqle",
                 weights.on.x = "hat", data = selected_data[-indexes,],
                 control = glmrobMqle.control(tcc = h_const[j], acc = 1e-8))
    pred = predict(fit, selected_data[indexes,-1])
    output[indexes] = ifelse(sigmoid(pred) > 0.5, 1, 0)
  }
  
  conf_matrix = table(Actual = selected_data[,1], Predicted = output)
  err[j] = (conf_matrix[1,2] + conf_matrix[2,1])/sum(conf_matrix)
  
}

h_opt = h_const[which.min(err)]
err_opt = min(err)

opt_fit = glmrob(Attrition ~ . , family = binomial(link=logit), method = "Mqle",
                 weights.on.x = "hat", data = selected_data, 
                 control = glmrobMqle.control(tcc = h_opt, acc = 1e-16))
summary(opt_fit)

opt_cov = as.matrix(opt_fit$cov)[-1,-1]
heatmap(opt_cov)
weights = opt_fit$w.r
hist(weights, breaks = seq(0, 1, by = 0.020))

####  Survival analysis --------------------------------------------------------

ds$Attrition <- factor(ds$Attrition, labels = (c('Censor', 'Event')))
id <- factor(seq(1, dim(ds)[1], 1))
ds <- data.frame(ds, id)

ggplot(data = ds[1:15,], aes(x = id,y = YearsAtCompany)) + geom_bar(stat = 'identity', width = 0.2) + 
  geom_point(aes(color = Attrition, shape = Attrition), size = 6) + coord_flip()

fit <- survfit(Surv(YearsAtCompany, Attrition =='Event') ~ 1, data = ds)
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
           data = ds)

H <- fit$cumhaz

ggsurvplot(fit,
           risk.table = TRUE,
           ggtheme = theme_bw(),
           fun='cumhaz',
           break.time.by=10,
           title="Cumulative Hazard Curve for Attrition",
           data = ds)

fit.gender <- survfit(Surv(YearsAtCompany, Attrition =='Event') ~ Gender, data=ds)

ggsurvplot(fit.gender, conf.int = F,
           risk.table = TRUE,
           risk.table.col = "strata", 
           ggtheme = theme_bw(),
           break.time.by=10,
           legend.labs=c("Female","Male"), legend.title="Gender class",  
           palette=c("violet","blue"), 
           title="Kaplan-Meier Curves by gender class for Attrition",
           data = ds)

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition =='Event') ~ Gender, data = ds)
log_rank_test
# p = 0.3 the difference is not significant in survival

hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[2]/log_rank_test$exp[2])
hazard_ratio
# 0.87 < 1
# the risk of attrition in female is 0.87 times the risk in male: being female is a protective factor.



ds$Department <- factor(ds$Department)

fit.department <- survfit(Surv(YearsAtCompany, Attrition == 'Event') ~ Department, data = ds)

ggsurvplot(fit.department, data = ds,
           risk.table = TRUE,
           risk.table.col = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_bw(),
           break.time.by=10,
           legend.labs=c("HR","RD","Sales"), legend.title="Department class",  
           palette=c("darkblue","red","forestgreen"), 
           title = "Kaplan-Meier Curves by department class for Attrition")

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ Department, data = ds)
log_rank_test
# p = 0.04 so the difference is significant

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ Department,
                          data = ds[c(which(ds$Department == "Human Resources"),which(ds$Department == "Sales")),])
log_rank_test
# p = 1, Sales and HR are the same

data_depart = ds
data_depart$Department = ifelse(ds$Department == "Sales", "Sales or HR", ifelse(ds$Department == "Human Resources", "Sales or HR", "RD"))

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ Department, data = data_depart)
log_rank_test
# p = 0.01, the two groups are different

hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[2]/log_rank_test$exp[2])
hazard_ratio
# 0.7 < 1
# the risk of attrition in Department RD is 0.7 times the risk in the others departments:
# being part of RD is a protective factor



ds$JobLevel <- factor(ds$JobLevel)

fit.Joblevel <- survfit(Surv(YearsAtCompany, Attrition == 'Event') ~ JobLevel, data = ds)

ggsurvplot(fit.Joblevel, data = ds,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           break.time.by=10,
           legend.labs=c("1","2","3","4","5"), legend.title="JobLevel class",  
           palette=c("darkblue","cyan3","green","purple","blue"), 
           title="Kaplan-Meier Curves by JobLevel class for Attrition")

log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ JobLevel, data=ds)
log_rank_test
# p < 2e-16 so the difference is very significant

# 2 vs 3
log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ JobLevel, 
                          data = ds[c(which(ds$JobLevel == 3),which(ds$JobLevel == 2)),])
log_rank_test
# p = 1, no difference

# 4 vs 5
log_rank_test <- survdiff(Surv(YearsAtCompany, Attrition=='Event') ~ JobLevel, 
                          data = ds[c(which(ds$JobLevel == 4),which(ds$JobLevel == 5)),])
log_rank_test
# p = 0.9, no substancial difference

data_joblevel = ds
data_joblevel$JobLevel = ifelse(ds$JobLevel == 1, "Group 1", ifelse(ds$JobLevel == 2, "Group 2/3", ifelse(ds$JobLevel == 3, "Group 2/3", "Group 4/5"))) 

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




ds$BusinessTravel <- factor(ds$BusinessTravel)
ds$Education <- factor(ds$Education)
ds$EducationField <- factor(ds$EducationField)
ds$Gender <- factor(ds$Gender)
ds$EnvironmentSatisfaction <- factor(ds$EnvironmentSatisfaction)
ds$JobInvolvement <- factor(ds$JobInvolvement)
ds$JobRole <- factor(ds$JobRole)
ds$JobSatisfaction <- factor(ds$JobSatisfaction)
ds$MaritalStatus <- factor(ds$MaritalStatus)
ds$OverTime <- factor(ds$OverTime)
ds$PerformanceRating <- factor(ds$PerformanceRating)
ds$RelationshipSatisfaction <- factor(ds$RelationshipSatisfaction)
ds$StockOptionLevel <- factor(ds$StockOptionLevel)
ds$WorkLifeBalance <- factor(ds$WorkLifeBalance)

glimpse(ds)

# Assumption : covariates time independent, so we removed YearsInCurrentRole, YearsSinceLastPromotion, 
#              YearsSinceLastPromotion, YearsWithCurrManager, TotalWorkingYears 

mod.cox <- coxph(Surv(YearsAtCompany, Attrition == 'Event') ~ 
                   Age + BusinessTravel + DistanceFromHome +
                   EnvironmentSatisfaction + Gender + JobInvolvement + I(JobLevel == 1) +
                   I(JobLevel == 2)  + I(JobLevel == 4)  +
                   I(JobRole == "Research Scientist") + I(JobRole == "Sales Executive") + 
                   I(JobRole == "Sales Representative") + I(JobSatisfaction == 4) + MonthlyIncome + 
                   NumCompaniesWorked + OverTime + RelationshipSatisfaction + 
                   I(StockOptionLevel == 0) + TrainingTimesLastYear + WorkLifeBalance , data = ds)
summary(mod.cox)

ggcoxdiagnostics(mod.cox, type = "martingale")

ggcoxdiagnostics(mod.cox, type = "deviance")

ggcoxdiagnostics(mod.cox, type = "schoenfeld")

ggcoxdiagnostics(mod.cox, type = "scaledsch")

test.ph <- cox.zph(mod.cox)
test.ph

for(i in 1:20){
  plot(test.ph[i])
  abline(h=0, col='red')
}

plot(test.ph[1], lwd = 2.5)
abline(h=0, col='red', lwd = 2.5)

plot(test.ph[14], lwd = 2.5)
abline(h=0, col='red', lwd = 2.5)


