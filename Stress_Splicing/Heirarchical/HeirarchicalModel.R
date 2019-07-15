### HEIARCHICAL MODEL ###

library(stringr)
library(dplyr)
library(pracma)
library(MASS)


# MONTH 1 (2018_6 / JUNE) CALIBRATED DATA FRAME 
# <<<<<<< HEAD

# # Mac Directory
# setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_6")
# deriv_complete=read.csv(file = "2018_6_1_qPCR_Output.csv", header=FALSE)
# deriv = deriv_complete
# # Remove extra labels column 
# deriv = deriv[,-1]
# # Transpose derivatives to be in equivalent format as raw plate data
# deriv = as.data.frame(t(deriv), header=FALSE)
# # Rename columns
# colnames(deriv)=c("plateID", "reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
# ### Removing NTC and gblock-Minus values ###
# # Indicate if sample is NTC (negative control)
# deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# # Remove NTC samples, indicator (T/F) column, and cpD2 values
# ntc = which(deriv$sampleID_NTC)
# deriv = deriv[-ntc,]
# deriv = deriv[,-c(6,7)]
# # Indicate if sample is 'Plus' or 'Minus'
# deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# # Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
# minus = which(deriv$sampleID_Minus)
# deriv = deriv[-minus,]
# deriv = deriv[,-6]
# deriv$cpD1 = as.numeric(as.character(deriv$cpD1))
# # Remove unusual observations from initial data frame (CT value less than 10)
# deriv = deriv %>% filter(deriv$cpD1 >= 10)
# # Create data frame for Calibrated values
# calib_data = deriv %>% filter(str_detect(sampleID, "g"))
# # Sort by starting quantity
# calib_data = calib_data[order(calib_data$starting_quantity),]
# calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
# calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))
# # Create empty vectors for for-loop to input cpD1 values
# test1 = c()
# allP = c()
# startq = c()
# # For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
# for(i in 1:length(calib_data$starting_quantity)){
#   sq <- calib_data$starting_quantity[i]
#   if(i %% 6 == 1){
#     startq = c(startq,sq,sq,sq)
#   }
#   val <- toString(calib_data$reaction_type[i])
#   if(strcmp(val, "test1")){
#     test1 = c(test1, calib_data$cpD1[i])
#   }
#   if(strcmp(val, "all_products")){
#     allP = c(allP, calib_data$cpD1[i])
#   }
# }
# # Bind test1 and allProd cpD1 values by starting quantity
# calib_data = as.data.frame(cbind(startq, test1, allP))
# # Format starting quantity values as decimals, not scientific notation
# calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
# calib_data$startq=as.factor(calib_data$startq)
# # Append ratios to data set
# calib_data_6 =calib_data
# # Create month indicator column
# calib_data_6$month = strrep('june', length(calib_data_6))
# =======
calib_data_6 = read.csv("../2018_6/calib_2018_6.csv")[,-1]
calib_data_6$month ='june'
#>>>>>>> 5f33e18e3fce60a5c0d91821babf3d5b9f6982a3
calib_data_6

# MONTH 2 (2018_8 / AUGUST) CALIBRATED DATA FRAME
calib_data_8 = read.csv("../2018_8/calib_2018_8.csv")[,-1]
calib_data_8$month ='aug'
calib_data_8

# MONTH 3 (2018_11 / NOVEMBER) CALIBRATED DATA FRAME
calib_data_11 = read.csv("../2018_11/calib_2018_11.csv")[,-1]
calib_data_11$month ='nov'
calib_data_11


# Combined Calib d.f. for all months
calib_data = rbind(calib_data_6, calib_data_8, calib_data_11)

# Create dummy varible columns for each month
calib_data$june = ifelse(str_detect(calib_data[,4], "june"), 1, 0)
calib_data$aug = ifelse(str_detect(calib_data[,4], "aug"), 1, 0)
calib_data = calib_data[,-4]

# Drop rows containing NA
calib_data = na.omit(calib_data)

# Calculating test1 and allp zscores
calib_data$ztest1 = (calib_data$test1 - mean(calib_data$test1))/sd(calib_data$test1)
calib_data$zallP = (calib_data$allP - mean(calib_data$allP))/sd(calib_data$allP)
calib_subset = calib_data[,c(1, 4:7)]

# Ordinal Logistic Regression Model 
model = polr(as.factor(calib_subset$startq) ~ ., data=calib_subset, Hess = TRUE)
#(summary(model))
(ctable <- coef(summary(model)))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
options(scipen=999)
## combined table
(ctable <- cbind(ctable, "p value" = p))


# OLRM - SQ ~ Test1
model1 = polr(as.factor(calib_subset$startq) ~ ztest1 + june + aug, data = calib_subset, Hess = TRUE)
#(summary(model))
(ctable <- coef(summary(model1)))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
options(scipen=999)
## combined table
(ctable <- cbind(ctable, "p value" = p))

# OLRM - SQ ~ allP
model2 = polr(as.factor(calib_subset$startq) ~ zallP + june + aug, data = calib_subset, Hess = TRUE)
#(summary(model))
(ctable <- coef(summary(model2)))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
options(scipen=999)
## combined table
(ctable <- cbind(ctable, "p value" = p))



## Ordinal Net package ##
library("ordinalNet")

#define ordinal model startq~ztest1+month
ordmod1 = ordinalNet(as.matrix(calib_subset[,2:4]), calib_subset$startq)
summary(ordmod1)
coef(ordmod1, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit1 = ordinalNetTune(as.matrix(calib_subse[,c(2:4)]), calib_subset$startq,
                         family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit1$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit1$loglik))
head(coef(ordfit1$fit, matrix = TRUE, whichLambda = bestLambdaIndex))
#interpretation???

#define ordinal model starq~zallp+month
ordmod2 = ordinalNet(as.matrix(calib_subset[,c(2,3,5)]), calib_subset$startq)
summary(ordmod2)
coef(ordmod2, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit2 = ordinalNetTune(as.matrix(calib_subset[,c(2,3,5)]), calib_subset$startq,
                         family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit2$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit2$loglik))
head(coef(ordfit2$fit, matrix = TRUE, whichLambda = bestLambdaIndex))
#interpretation???

#define ordinal model starq~zallP+ztest1+month
ordmod3 = ordinalNet(as.matrix(calib_subset[,4:5]), calib_subset$startq)
summary(ordmod3)
coef(ordmod3, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit3 = ordinalNetTune(as.matrix(calib_subset[,4:5]), calib_subset$startq, family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit3$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit3$loglik))
head(coef(ordfit3$fit, matrix = TRUE, whichLambda = bestLambdaIndex))
#interpretation???



# #### Month 1 (2018_6) CT ####
# setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_6")
# m11 = read.csv(file = "2018_6_1_qPCR_Output.csv", header=FALSE)
# #Initial framing
# m11 = m11[,-1]
# # Transpose derivatives to be in equivalent format as raw plate data
# m11 = as.data.frame(t(m11), header=FALSE)
# m11 = m11[,-6]
# # Rename columns
# colnames(m11)=c("plateID", "reaction_type", "sampleID", "starting_quantity", "cpD1")
# # Remove NTC, gblock_minus values
# m11['sampleID_NTC'] = grepl('NTC', m11$sampleID)
# ntc = which(m11$sampleID_NTC)
# m11 = m11[-ntc,]
# m11 = m11[,-6]
# m11['sampleID_Minus'] = grepl('minus', m11$sampleID)
# minus = which(m11$sampleID_Minus)
# m11 = m11[-minus,]
# m11 = m11[,-6]
# m11$cpD1 = as.numeric(as.character(m11$cpD1))
# # Remove CT vals <10
# m11 = m11 %>% filter(m11$cpD1 >= 10)
# # Create calibrated data frame
# calib_data = m11 %>% filter(str_detect(sampleID, "g"))

# 
# #Month 2 (2018_8) CT
# setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_8")
# m21 = read.csv(file = "2018_8_1_qPCR_Output.csv", header=FALSE)
# m21 = m21[,-1]
# m22 = read.csv(file = "2018_8_2_qPCR_Output.csv", header=FALSE)
# m22 = m22[,-1]
# m23 = read.csv(file = "2018_8_3_qPCR_Output.csv", header=FALSE)
# m23 = m23[,-1]
# 
# #Month 3 (2018_11) CT
# setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_11")
# m31 = read.csv(file = "2018_11_1_qPCR_Output.csv", header=FALSE)
# m31 = m31[,-1]
# m32 = read.csv(file = "2018_11_2_qPCR_Output.csv", header=FALSE)
# m32 = m32[,-1]
# 
# #Combine months
# month = data.frame(m11, m21, m22, m23, m31, m32)
# month = as.data.frame(t(month))
# month$june = ifelse(str_detect(month[,1], "2018_6"), 1, 0)
# month$V1 = str_replace(month$V1, "2018_9", "2018_8")
# month$aug = ifelse(str_detect(month[,1], "2018_8"), 1, 0)
# 
# #Remove CT<10, NTC, gblock_minus, V6 (CT2)
# month$V5 = as.numeric(as.character(month$V5))
# month$V6 = as.numeric(as.character(month$V6))
# month = month %>% filter(month$V5 >= 10)
# month = month %>% filter(month$V6 >= 10)
# month['NTC'] = grepl('NTC', month$V3)
# ntc = which(month$NTC)
# month = month[-ntc,]
# month = month[,-9]
# month['Minus'] = grepl('minus', month$V3)
# minus = which(month$Minus)
# month = month[-minus,]
# month = month[,-c(6,9)]

