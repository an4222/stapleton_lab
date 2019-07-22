### generated test data for qpcr adjustment ###
library('dplyr')
library('ordinalNet')

### simulation n = 21 ###
#create the dataset with starting quantites, all products and test 1 cp values
sd.val = 0.2
sd.val.test = 0.1
eps = rnorm(21, 0, sd.val)
eps.test= rnorm(21, 0, sd.val.test)
stq =sort(rep(c(0.010, 0.050, 0.100, 0.500, 1.000, 0.005, 5.000), 3))
allp = (-2.5)*log(stq)+25+eps
test1 = -2.7*log(stq)+23+eps.test

#plot the two cp types against log starting quantity
plot( allp,log(stq), col = 'red', main = "Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(15, 45), ylim = c(-6, 2))
points(test1,log(stq), col = 'blue')
abline(lm(log(stq)~allp), col = 'red')
abline(lm(log(stq)~test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)


# create the data frame and include the z scores of test 1 and all p
data = as.data.frame(cbind(stq, allp, test1))
data$ztest1 = (data$test1 - mean(data$test1))/sd(data$test1)
data$zallP = (data$allp - mean(data$allp))/sd(data$allp)


#### Ordinal Net Models ####

## Model 1 ##
#define ordinal model starq~zallP+ztest1
covmat = as.matrix(data[, c(4,5)])
ordmod = ordinalNet(covmat, as.factor(data$stq))
summary(ordmod)
coef(ordmod, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit = ordinalNetTune(covmat, as.factor(data$stq), family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit$loglik))
head(coef(ordfit$fit, matrix = TRUE, whichLambda = bestLambdaIndex))


# ##finding ajustment value##
# group = split.data.frame(data, data$stq)
# 
# adj <- function(AllP, Test1){
#   adjust = ave(AllP)-ave(Test1)
#   return(adjust)
# }
# 
# adjval = NULL
# for (k in group){
#   adjval = c(adjval,adj(k$allp, k$test1))
# }
# 
# data$adjval = adjval
# data$adjusted_test1 = data$test1 + adjval
# 
# 
# plot( allp,log(stq), col = 'red', main = "Log Starting Quanitty vs. Cp Values", 
#      ylab = "log(Starting Quantity)", xlab = "Cp Values", xlim= c(10, 45), ylim = c(-6, 2))
# points(data$adjusted_test1,log(stq), col = 'blue')
# abline(lm(log(stq)~allp), col = "red")
# abline(lm(log(stq)~data$adjusted_test1), col = "blue")
# legend('topright', legend=c("Adjusted Test 1", "All Products"),
#        col=c("blue", "red"), lty = 1, cex=0.8)

##########################################################################################

### simulation n = 99

#create the dataset with starting quantites, all products and test 1 cp values
sd.val = 1
sd.val.test = 0.5
eps = rnorm(21, 0, sd.val)
eps.test= rnorm(21, 0, sd.val.test)
stq =sort(rep(c(0.010, 0.050, 0.100, 0.500, 1.000, 0.005, 5.000), 33))
allp = (-2.5)*log(stq)+25+eps
test1 = -2.7*log(stq)+23+eps.test

#plot the two cp types against log starting quantity
plot( allp,log(stq), col = 'red', main = "Log Starting Quanitty vs. Cp Values", 
      ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(15, 45), ylim = c(-6, 2))
points(test1,log(stq), col = 'blue')
abline(lm(log(stq)~allp), col = 'red')
abline(lm(log(stq)~test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)



# create the data frame and include the z scores of test 1 and all p
data = as.data.frame(cbind(stq, allp, test1))
data$ztest1 = (data$test1 - mean(data$test1))/sd(data$test1)
data$zallP = (data$allp - mean(data$allp))/sd(data$allp)


## Model 2 ##
#define ordinal model starq~zallP+ztest1
covmat = as.matrix(data[, c(4,5)])
ordmod = ordinalNet(covmat, as.factor(data$stq))
summary(ordmod)
coef(ordmod, matrix=TRUE)
#kfold cv
set.seed(123)
ordfit = ordinalNetTune(covmat, as.factor(data$stq), family = "cumulative",
                        link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                        warn = FALSE, printProgress = FALSE)
head(ordfit$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit$loglik))
head(coef(ordfit$fit, matrix = TRUE, whichLambda = bestLambdaIndex))

