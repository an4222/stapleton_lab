########################################################## 
############## QPCR PLATE & ADJUSTMENT MODEL #############
########################################################## 

library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
library(MASS)
library(glm.predict)
library(Stack)

# Mac Directory
setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_6")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")
# PC Directory
#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing/qPCR")

### READ IN DERIVATIVE DATA ###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
#deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
#deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
#deriv=cbind(deriv.1, deriv.2)

# In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)
deriv_complete=read.csv(file = "2018_6_1_qPCR_Output.csv", header=FALSE)

########################################################## 
################### Initial Data Framing #################
########################################################## 

deriv = deriv_complete
# Remove extra labels column 
deriv = deriv[,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=FALSE)
# Remove blank column (4th)
#deriv = deriv[,-5]
# Rename columns
colnames(deriv)=c("plateID", "reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
### Removing NTC and gblock-Minus values ###
# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv = deriv[,-c(6,7)]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
deriv = deriv[-minus,]
deriv = deriv[,-6]
deriv$cpD1 = as.numeric(as.character(deriv$cpD1))

### COMPLETED INITIAL DATA FRAMING ###

########################################################## 
############ Removing Ununsual Observations ##############
##########################################################

# Remove unusual observations from initial data frame (CT value less than 10)
deriv = deriv %>% filter(deriv$cpD1 >= 10)
# Read in raw cycle data - may need to combine multiple files
cycle1 = read.csv(file = "2018_6_1_plate.csv", header = FALSE)
# Create complete set of reaction data (derivative and cycle)
reaction = Stack(deriv_complete, cycle1)
# Remove repeat labeling
replace = reaction[7:10,]
reaction = reaction[-c(1:4, 7:10),]
reaction = Stack(replace, reaction)
# Transpose so column headers at top
reaction = as.data.frame(t(reaction))
reaction = reaction[,-c(6:7)]
# Replace column names with first row
colnames(reaction) <- as.character(unlist(reaction[1,]))
reaction = reaction[-1,]
colnames(reaction)[5] = "cpD1"
reaction$cpD1 = as.numeric(as.character(reaction$cpD1))
# Filter unusual observations (CT value less than 10)
unusual_obs_2018_6 = reaction %>% filter(reaction$cpD1 < 10)
# Write CSV file 
#write.csv(unusual_obs_2018_6, file="Unusual_Obs_2018_6.csv")

# ### COMPLETED UNUSUAL OBSERVATIONS REMOVAL/REPORTING ###


########################################################## 
################# Calibrated Data Framing ################
########################################################## 
library("rowr")
# Create/Write data frame for Calibrated values
calib_data = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_data = calib_data[order(calib_data$starting_quantity),]

calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))


test1 = filter(calib_data, reaction_type=="test1")[,5]
allP = filter(calib_data, reaction_type=="all_products")[,4:5]

#Combine test1 and allP obs, with NA in blank cells
calib_data = as.data.frame(cbind.fill(allP, test1, fill = NA))
colnames(calib_data) = c("startq", 'allP', "test1" )
# Format starting quantity values as decimals, not scientific notation
calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq=as.factor(calib_data$startq)

#Apply log scale to test1 and allP CT values
calib_data$allPln = log(calib_data$allP)
calib_data$test1ln = log(calib_data$test1)

write.csv(calib_data, file = "calib_2018_6.csv")


### COMPLETED CALIBRATED DATA FRAME ###

########################################################## 
############### Experimental Data Framing ################
########################################################## 

# Create/Write data frame for Experimental values
exp_data = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_data = exp_data[order(exp_data$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
# exp_data = exp_data[-1,]
# exp_data = exp_data[-nrow(exp_data),]
exp_data$cpD1 = as.numeric(as.character(exp_data$cpD1))
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find counts of each unique sampleID; for sample with a count not equal to 2, remove from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
# Remove invalid observations from data set
exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp = toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID, convert to data frame
exp_data = as.data.frame(cbind(sampleID.exp, test1.exp, allP.exp))
exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))

# Apply log scale to test1 and allP CT values
exp_data$test1.exp.ln = log(exp_data$test1.exp)
exp_data$allP.exp.ln = log(exp_data$allP.exp)

write.csv(exp_data, file = "exp_2018_6.csv")
### COMPLETED EXPERIMENTAL DATA FRAME ###



### CONFUSTION MATRIX ###
library(caret)
numLlvs <- 4
confusionMatrix(
  factor(sample(rep(letters[1:numLlvs], 200), 50)),
  factor(sample(rep(letters[1:numLlvs], 200), 50))) 


##########################################################
########## PROBABILITY MODEL - Calibrated Data ###########
##########################################################

##finding ajustment value##
group = split.data.frame(calib_data, calib_data$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(k$allP, k$test1))
}

calib_data$adjval = adjval
calib_data$adjusted_test1 = calib_data$test1 + adjval

##### Finding the average adjusted test 1 #########

#average function takes the mean of each 
average <- function(col1){
  avg = NULL;
  for (i in col1){
    avg = c(avg,mean(col1))
  }
  return(avg)
}
##Subset data by starting quantity
group = split.data.frame(calib_data, calib_data$startq)

adj.test1.avg = NULL;
for (k in group){
  adj.test1.avg = c(adj.test1.avg, average(k$adjusted_test1))
}
print(adj.test1.avg)

calib_adj = unique(as.data.frame(cbind(as.numeric(as.character(calib_data$startq)), adj.test1.avg)))

# Rename columns
colnames(calib_adj)=c("startq", "adj.test1.avg")
#############################################

##### Creating the components of the ORLM ######
# Calculate z-score for calibrated data
zscore = (newratiosvector - mean(newratiosvector))/sd(newratiosvector)
newratios.calib$zscore = zscore
# Predict calibrated data ratios using experimental data
pred.ratio = zscore*sd(ratio.exp)+mean(ratio.exp) #why do we need this???
newratios.calib$pred.ratio = pred.ratio
# Append y (predicted calibrated ratios) to calibrated data frame -- CALIBRATED RATIOS IN TERMS OF EXPERIMENTAL PARAMETERS
#no way of knowing if the pred.ratio and the correct starting quantities are lining up correctly #
#calib_data = cbind(calib_data,newratiosvector, pred.ratio) 

# Ordinal Logistic Regression Model - starting quantity as response to calibrated z-score
model = polr(as.factor(newratios.calib$startqvector) ~ newratios.calib$pred.ratio, Hess = TRUE)
# should the model be based off of the relation between startq and pred.ratio or zscore?
model = polr(as.factor(newratios.calib$startqvector) ~ newratios.calib$zscore, Hess = TRUE)
summary(model)
(ctable <- coef(summary(model)))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
options(scipen=999)
## combined table
(ctable <- cbind(ctable, "p value" = p))
#############################################

# Calculate experimental data z-score
zscore = (exp_data$ratio.exp - mean(exp_data$ratio.exp))/sd(exp_data$ratio.exp)
prob.matrix = predict(model, zscore, type='p')

# Apply probability matrix to the adjusted test 1 averages
apply(prob.matrix, 1, function(x) x*calib_adj$adj.test1.avg)
exp_data$exp.adjust = colSums(apply(prob.matrix, 1, function(x) x*calib_adj$adj.test1.avg))

# Create new column with stress product (VQTL input)
exp_data$stress = exp_data$allP.exp - exp_data$exp.adjust


###PLOTS###
# Calibrated data - s.q. vs. ratio
plot(as.factor(newratios.calib$startqvector), as.numeric(newratios.calib$newratiosvector), xlab='Starting Quantity', ylab='Ratio', 
     main='2018_6 Calibrated Data - Starting Quantities vs. Ratios')


# plot(as.factor(newratios.calib$startqvector), as.numeric(newratios.calib$zscore), xlab='Starting Quantity', ylab='Ratio', 
#      main='2018_6 Calibrated Data - Starting Quantities vs. Ratios')
# Histogram - calib allP vs. exp allP
hist(calib_data$allP, xlim=c(0,40), ylim=c(0,100), col=rgb(1,0,0,0.5), main='2018_6 Histogram of All Products', xlab='All Products Derivative')
hist(exp_data$allP.exp, xlim=c(0,40), ylim=c(0,100), add=T, col=rgb(0,0,1,0.5))
legend("topleft",
       c("Calibration", "Experimental"),
       fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
#dev.off()
# Histogram - calib test1 vs. exp test1
hist(calib_data$test1, xlim=c(0,30), ylim=c(0,30), col=rgb(1,0,0,0.5), main='2018_6 Histogram of Test 1', xlab='Test 1 Derivative')
hist(exp_data$test1.exp, xlim=c(0,30), ylim=c(0,30), add=T, col=rgb(0,0,1,0.5))
legend("topleft",
       c("Calibration", "Experimental"),
       fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
dev.off()


###UPDATE NONPAIRWISE RATIO PLOTS###
#Ratios - Calibrated#
hist(newratios.calib$combratio, xlim=c(0,3), ylim=c(0,70), col=rgb(1,0,0,0.5), main='Histogram of Non-Pairwise Ratios', xlab='Ratio')
#Ratios - Experimental#
# Values excluded from histogram that will be further investigated later (their index)
x=c(8,82,141,148,149,153,161,170,172,175,180,188)
exp_data2=exp_data[-x,]
hist(exp_data2$ratio.exp, xlim=c(0,3), ylim=c(0,70), col=rgb(0,0,1,0.5), add=T)
legend("topleft",
       c("Calibration", "Experimental"),
       fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")


#######################################
####### k-Fold Cross Validation #######
#######################################
library(caret)


n = length(newratios.calib[,1])
set.seed(1)
f <- 6
folds <- rep_len(1:f, length.out = dim(newratios.calib)[1])
folds <- sample(folds, size = dim(newratios.calib)[1], replace = F)
table(folds)


preds = NULL
truth = NULL
for (k in 1:f){
  test.ID = which(folds == k)
  train_y = newratios.calib[-test.ID, "startqvector"]
  train_x = newratios.calib[-test.ID, "zscore"]
  test_y = newratios.calib[test.ID, "startqvector"]
  test_x = newratios.calib[test.ID, "zscore"]
  truth = c(truth, test_y)
  model.fit = polr(as.factor(train_y) ~ train_x, Hess = TRUE) 
  preds = c(preds, predict(model.fit, test_x)) #need to fix s.t. results in 9 values rather than 45
}

confusionMatrix(preds,truth)


###### OLD CODE #######


# # Adjustment: allP - test1
# calib_data$diff = calib_data$allP - calib_data$adjusted_test1
# 
# calib_data = cbind(calib_data, calib.zscore)
# plot(as.factor(calib_data$startq), calib_data$calib.zscore)
# 
# # Ordinal Logistic Regression Model - starting quantity as response to calibrated z-score
# model = polr(as.factor(startqvalues) ~ calib.zscore, Hess = TRUE)
# summary(model)
# # Calculate experimental data z-score
# exp_data$exp.zscore = (exp_data$ratio.exp - mean(exp_data$ratio.exp))/sd(exp_data$ratio.exp)
#prob.matrix = predict(model, exp_data$exp.zscore, type='p') 
# apply(prob.matrix, 1, function(x) x*calib_data$diff)
# exp_data$VQTL = colSums(apply(prob.matrix, 1, function(x) x*calib_data$diff))



#Predict s.q. using ordinal logistic mode
#startq~Zcalb
#pred(model~zexp)
### just using z scores to predict

###### use the probability starting quantity matrix from the predict function
#to find the weighted average 

# ########################################################## 
# #### PREDICT STARTING QUANTITIES - Experimental Data ####
# ########################################################## 
# 
# # Using the adjustment model on the expiremental data
# new = data.frame((ratio = exp_data$allP/exp_data$test1), sampleID.exp)
# exp_predict_sq = as.data.frame(predict(OLR, new, interval = "confidence"))
# # Append sample ID's and corresponding starting quantity predictions
# exp_predict_sq = cbind(exp_predict_sq, exp_data$sampleID.exp)
# # Rename columns, re-order
# colnames(exp_predict_sq)=c("predicted_sq", "sampleID.exp")
# exp_predict_sq = exp_predict_sq[c(2,1)]
# #exp_predict_sq$predicted_sq = as.numeric(as.character(exp_predict_sq$predicted_sq))
# #exp_predict_sq$sampleID.exp = as.numeric(as.character(exp_predict_sq$sampleID.exp))
# # Merge complete experimental data frame with predicted starting quantities data frame by sample ID
# exp_data_predict = merge.data.frame(exp_data, exp_predict_sq, by="sampleID.exp")
# 
# ### COMPLETED ADJUSTMENT MODEL - EXPERIMENTAL DATA ###
# 
# ########################################################## 
# #### ADJUSTMENT VALUE BY STARTING QUANTITY ####
# ########################################################## 
# 
# # Create new data frame containing only predicted s.q. and adj_val
# sq.adjval = as.data.frame(cbind(as.numeric(as.character(calib_data$startq)), calib_data$adj_val))
# sq.adjval$V1=as.factor(format(sq.adjval$V1, scientific=FALSE))
# colnames(sq.adjval)=c("predicted_sq", "adj_val")
# # Remove repeat rows
# sq.adjval = sq.adjval[!duplicated(sq.adjval), ]
# # Match adjustment value for each s.q. with corresponding predicted s.q. in experimental data frame
# exp_data_predict = merge(exp_data_predict, sq.adjval, by='predicted_sq')
# exp_data_predict = exp_data_predict[order(exp_data_predict$sampleID.exp),]
# exp_data_predict = exp_data_predict[c(2,3,4,5,6,1)]
# # 
# #Calculate z-score for calibrated data
# calib.zscore = (calib_data$ratio - mean(calib_data$ratio))/sd(calib_data$ratio)
# #Predict calibrated data ratios using experimental data
# pred.ratio = calib.zscore*sd(ratio.exp)+mean(ratio.exp)
# #Append y (predicted calibrated ratios) to calibrated data frame
# calib_data = cbind(calib_data, pred.ratio)
# # Create empty vectors for for-loop input
# calib_data = as.data.frame(calib_data)
# calib_data$test1 = as.numeric(as.character(calib_data$test1))
# calib_data$allP = as.numeric(as.character(calib_data$allP))
# adj_val = c()
# allP = c()
# startq = c()
# ratio =calib_data$allP/calib_data$test1
# # Itterating through each set of (3) observations performing U-Stats on each set of inputs
# for (i in 1:(nrow(calib_data)/3)){
#   t_x <- c(calib_data$allP[3*i - 2], calib_data$allP[3*i - 1], calib_data$allP[3*i])
#   t_y <- c(calib_data$test1[3*i - 2], calib_data$test1[3*i - 1], calib_data$test1[3*i])
#   adj <- mean(outer(t_x, t_y, "-"))
#   adj_val <- c(adj_val, adj, adj, adj)
# }
# adjusted_test1 <- test1 + adj_val
# # Append adjusted test1 values and adjustment value to data set
# calib_data=cbind(calib_data,adjusted_test1,adj_val)
# # Write Calibrated Data CSV --> Used in "qPCR_Plotting" code for visuals
# #write.csv(file="YEAR_MONTH_Calibrated_DF", calib_data)





