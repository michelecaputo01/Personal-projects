library(fastDummies)
library(MASS)
library(rgl)
library(DepthProc)
library(readxl)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(corrplot)
library(glmperm)
library(lmtest)
library(robustbase)

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

