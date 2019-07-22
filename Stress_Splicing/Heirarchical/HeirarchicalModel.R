### HEIARCHICAL MODEL ###

library(MASS)
library(stringr)
## Ordinal Net package ##
library("ordinalNet")

# MONTH 1 (2018_6 / JUNE) CALIBRATED DATA FRAME 
calib_data_6 = na.omit(read.csv("../2018_6/calib_2018_6.csv")[,-1])
calib_data_6$ztest1 = (calib_data_6$test1 - mean(calib_data_6$test1))/sd(calib_data_6$test1)
calib_data_6$zallP = (calib_data_6$allP - mean(calib_data_6$allP))/sd(calib_data_6$allP)
calib_data_6$month ='june'
#calib_data_6

# MONTH 2 (2018_8 / AUGUST) CALIBRATED DATA FRAME
calib_data_8 = na.omit(read.csv("../2018_8/calib_2018_8.csv")[,-1])
calib_data_8$ztest1 = (calib_data_8$test1 - mean(calib_data_8$test1))/sd(calib_data_8$test1)
calib_data_8$zallP = (calib_data_8$allP - mean(calib_data_8$allP))/sd(calib_data_8$allP)
calib_data_8$month ='aug'
#calib_data_8

# MONTH 3 (2018_11 / NOVEMBER) CALIBRATED DATA FRAME
calib_data_11 = na.omit(read.csv("../2018_11/calib_2018_11.csv")[,-1])
calib_data_11$ztest1 = (calib_data_11$test1 - mean(calib_data_11$test1))/sd(calib_data_11$test1)
calib_data_11$zallP = (calib_data_11$allP - mean(calib_data_11$allP))/sd(calib_data_11$allP)
calib_data_11$month ='nov'
#calib_data_11

# Combined Calib d.f. for all months
calib_data = rbind(calib_data_6, calib_data_8, calib_data_11)

# Create dummy varible columns for each month
calib_data$june = ifelse(str_detect(calib_data[,6], "june"), 1, 0)
calib_data$aug = ifelse(str_detect(calib_data[,6], "aug"), 1, 0)
calib_data$nov = ifelse(str_detect(calib_data[,6], "nov"), 1, 0)
calib_data = calib_data[,-6]

# Drop rows containing NA
calib_data = na.omit(calib_data)
calib_subset = calib_data[,c(1, 4:8)]

#graphing log starting quantity to cp values
plot(calib_data$allP,log(calib_data$startq), col = 'red', main = "Hierarchical Log Starting Quanitty vs. Cp Values", 
      ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data$test1,log(calib_data$startq), col = 'blue')
abline(lm(log(calib_data$startq)~calib_data$allP), col = 'red')
abline(lm(log(calib_data$startq)~calib_data$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)


#graphing log starting quantity to cp values by month
#June
plot(calib_data_6$allP,log(calib_data_6$startq), col = 'red', main = "June Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data_6$test1,log(calib_data_6$startq), col = 'blue')
abline(lm(log(calib_data_6$startq)~calib_data_6$allP), col = 'red')
abline(lm(log(calib_data_6$startq)~calib_data_6$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)

#August
plot(calib_data_8$allP,log(calib_data_8$startq), col = 'red', main = "August Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data_8$test1,log(calib_data_8$startq), col = 'blue')
abline(lm(log(calib_data_8$startq)~calib_data_8$allP), col = 'red')
abline(lm(log(calib_data_8$startq)~calib_data_8$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)

#November
plot(calib_data_11$allP,log(calib_data_11$startq), col = 'red', main = "November Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data_11$test1,log(calib_data_11$startq), col = 'blue')
abline(lm(log(calib_data_11$startq)~calib_data_11$allP), col = 'red')
abline(lm(log(calib_data_11$startq)~calib_data_11$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)


###### POLR models ######
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


###### Ordinal Net models #######


### model 1 ###
#define ordinal model startq~ztest1+month
ordmod1 = ordinalNet(as.matrix(calib_subset[,2:4]), as.factor(calib_subset$startq))
summary(ordmod1)
coef(ordmod1, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit1 = ordinalNetTune(as.matrix(calib_subset[,c(2:4)]), as.factor(calib_subset$startq),
                         family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit1$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit1$loglik))
head(coef(ordfit1$fit, matrix = TRUE, whichLambda = bestLambdaIndex))

### model 2 ###
#define ordinal model starq~zallp+month
ordmod2 = ordinalNet(as.matrix(calib_subset[,c(2,3,5)]), as.factor(calib_subset$startq))
summary(ordmod2)
coef(ordmod2, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit2 = ordinalNetTune(as.matrix(calib_subset[,c(2,3,5)]), as.factor(calib_subset$startq),
                         family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit2$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit2$loglik))
head(coef(ordfit2$fit, matrix = TRUE, whichLambda = bestLambdaIndex))

### model 3 ###
#define ordinal model starq~zallP+ztest1+month
ordmod3 = ordinalNet(as.matrix(calib_subset[,2:5]), as.factor(calib_subset$startq))
summary(ordmod3)
coef(ordmod3, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit3 = ordinalNetTune(as.matrix(calib_subset[,2:5]), as.factor(calib_subset$startq), family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit3$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit3$loglik))
head(coef(ordfit3$fit, matrix = TRUE, whichLambda = bestLambdaIndex))






##### building separate month models into one #####
ordmod6 = ordinalNet(as.matrix(calib_data_6[,4:5]), as.factor(calib_data_6$startq))
summary(ordmod3)
coef(ordmod3, matrix=TRUE)

ordmod8 = ordinalNet(as.matrix(calib_data_8[,4:5]), as.factor(calib_data_6$startq))
summary(ordmod3)
coef(ordmod3, matrix=TRUE)

ordmod11 = ordinalNet(as.matrix(calib_data_11[,4:5]), as.factor(calib_data_6$startq))
summary(ordmod3)
coef(ordmod3, matrix=TRUE)










#######################################
####### k-Fold Cross Validation ####### ordinalNet model
#######################################
library(caret)


n = length(calib_subset[,1])
set.seed(1)
f <- 5
folds <- rep_len(1:f, length.out = dim(calib_subset)[1])
folds <- sample(folds, size = dim(calib_subset)[1], replace = F)
table(folds)


preds = NULL
truth = NULL
for (k in 1:f){
  test.ID = which(folds == k)
  train_y = calib_subset[-test.ID, "startq"]
  train_x = calib_subset[-test.ID, 2:5]
  test_y = calib_subset[test.ID, "startq"]
  test_x = calib_subset[test.ID, 2:5]
  truth = c(truth, test_y)
  model.fit = ordinalNet(as.matrix(train_x), as.factor(train_y)) 
  preds = c(preds, predict(model.fit, test_x)) #need to fix s.t. results in 9 values rather than 45
}

confusionMatrix(preds,truth)
ordinalNet(as.matrix(calib_subset[,2:5]), as.factor(calib_subset$startq))

