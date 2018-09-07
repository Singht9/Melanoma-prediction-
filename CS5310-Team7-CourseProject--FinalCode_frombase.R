
# Identifying genetic markers in skin cancer cells #
# ================================================= #
# Team 7 - Emre Karakus | Tripti Singh
# Received 23689 X 4645 data. Identified two sets of 33 genes each : data_set1 and data_set2
# Classifiers developed on Genetic markers in data_set1 are 98-99 % accurate with similar kappa values
"planning on publishing the results to attract more scrutiny from subject matter experts"
"planning on publishing a paper on extracting genetic markers to classify Malignant & benign cells, based on our learnings during this project"


# STEP 1] READ DATA IN R -> basedata :==========================================================================
# Read basedata in csv file = GSE72056.csv =====================================================================

basedata <- read.table("C:/Users/Amit/Downloads/Documents/Final Project Files/GSE72056.csv", header=TRUE, sep=",",stringsAsFactors = FALSE)
str(basedata)

# STEP II] transpose to convert in to Row X Column dataframe -> basedata_t :====================================
#===============================================================================================================
basedata_t <- as.data.frame(t(basedata[-1]),stringsAsFactors=FALSE) 
str(basedata_t)

head(basedata_t[1:5,1:5]) # crosscheck

names(basedata_t) <- t(basedata[1]) # add column names = Genes
str(basedata_t)  

#write.csv(basedata_t,file="transposedBaseData.csv") # 4645 X 23689

# STEP III] remove '0' value columns : ==========================================================================
# basedata_t_V01 4645 X 22848 ===================================================================================

## checking '0' values using column sum = col. average ( as shown below)----
sum_vector <- 0  

for (k in 1:23689){
  sum_basedata_t <- sum(basedata_t[k])
  sum_vector <- c(sum_basedata_t==0,sum_vector)
}
table(sum_vector)
# Find columns with only 0 values using average ----- crosscheck -----------
sum_vector <- 0  

for (k in 1:23689){
  sum_basedata_t <- mean(basedata_t[,k])
  sum_vector <- c(sum_basedata_t==0,sum_vector)
}
table(sum_vector)


basedata_t_V01 <- basedata_t[,-(which(colSums(basedata_t) == 0))] 
str(basedata_t_V01)

# STEP IV] remove 3rd column listing non -malignant cell types ==================================================
# basedata_t_V02 4645 X 22847 ===================================================================================

basedata_t_V02 <- basedata_t_V01[-3]
str(basedata_t_V02)

# quick check to see understand data ##############
range(sapply(basedata_t_V02[,3:50], mean)) 

range(sapply(basedata_t_V02[,3:50], sd) )

table(basedata_t_V02[2])

Melanoma_data <- basedata_t_V02 ########################## henceforth will be using 4645 X 22847 <- Melanoma_data

table(Melanoma_data[2])

Melanoma_data$`malignant(1=no,2=yes,0=unresolved)`<- factor(Melanoma_data$malignant)

# STEP V] remove 'unresolved' malignant cells= keep only Benign and Malignant ===================================
# remove_unresolved 4513 X 22846 ===================================================================================

remove_unresolved <- subset(Melanoma_data,Melanoma_data$`malignant(1=no,2=yes,0=unresolved)` != 0) # only 0 and 1

table(remove_unresolved[2])

remove_unresolved$`malignant(1=no,2=yes,0=unresolved)` <- factor(remove_unresolved$`malignant(1=no,2=yes,0=unresolved)`,
                                                                 levels = c(1,2), labels = c("Benign","Malignant")) 

str(remove_unresolved)

## name numeric and class variables for convenience's sake ##
factor_vector <- remove_unresolved$`malignant(1=no,2=yes,0=unresolved)` ####################
num_df <- remove_unresolved[-2] ############################################################

# STEP VI] remove 'nearZeroVar' attributes =========================================================================
# cleanBase 4513 X 5626 ===================================================================================

library(ggplot2)
library(caret)

var_col <- nearZeroVar(num_df,saveMetrics = FALSE)
str(var_col)

Var_filtered <- num_df[,-var_col]   # entire dataset,NOTincluding factor variable, with nrZeroVar removed
str(Var_filtered)

cleanBase <- cbind(factor_vector,Var_filtered)
str(cleanBase)

# write.csv(cleanBase,file="cleanBase.csv")

factor_cleanBase <- cleanBase$factor_vector ##################################################
num_cleanBase <- cleanBase[-1]################################################################
str(num_cleanBase)############################################################################

# STEP VI] correlated predictors removed =========================================================================
# reduced_data 4513 X 5564 =======================================================================================
# reading data in to a csv for furher modelling - DATA PREPRATION STEPS CONCLUDED ================================

cor_filtered <- cor(num_cleanBase)
summary(cor_filtered[upper.tri(cor_filtered)])
remove_cor <- findCorrelation(cor_filtered,cutoff = 0.75)

str(remove_cor)

filtered_cor <- num_cleanBase[,-remove_cor]
str(filtered_cor)

reduced_data <- cbind(factor_vector,filtered_cor)
str(reduced_data)

# write.csv(reduced_data,file="reduced_data.csv")

####################  reduced_data used for further modelling #####################################################
# ref. PCA : https://little-book-of-r-for-multivariate-analysis.readthedocs.io/en/latest/src/multivariateanalysis.html
# ref. Caret: https://topepo.github.io/caret/parallel-processing.html
##################################################################################################################

### MODELLING :
# [A] reading data prepared for modelling and furher processing,as per model requirement

melanoma <- read.table("reduced_data.csv", header=TRUE, sep=",",stringsAsFactors = FALSE)
str(melanoma) # base data

factor_redData <- melanoma$factor_vector # 4513 X 1 class variable
num_redData <- melanoma[-c(1:2)] # 4513 X 5564 numeric attributes only(-cell name and class variable)
str(num_redData) 

baseData <- cbind(factor_redData,num_redData) # df we are using for further analysis = baseData
str(baseData)

# [B] CULLING OUT data_set1 from baseData using rpart -------------------------------------------------------------
############################################ Let's run rpart ------------------------------------------------------
library(rpart)
library(gmodels)
rpart_melanoma <- rpart(factor_redData~., data=baseData)
pred_melanoma <- predict(rpart_melanoma,baseData[1:100,],type = "class")
CrossTable(factor_redData[1:100],pred_melanoma) # 98% accurate with 2 False positives on first 100 datapoints

library(rpart.plot)
rpart.plot(rpart_melanoma, digits = 2)

rpart_melanoma$variable.importance # rpart is running on 5565 variables and getting a list of 33 imp. variables

Var_rpart <- rpart_melanoma$variable.importance

rpart_melanoma$cptable

reduced_rpart <- baseData[,names(baseData) %in% c(strsplit(names(Var_rpart), " "))] # df with 33 variables

##############################################################################################################
##############################################################################################################
data_set1 <- cbind(factor_redData,reduced_rpart) # data frame 4513 X 34 ==== will use for further analysis ===
##############################################################################################################
##############################################################################################################

# ==========================================================================================
# Developing a blind-holdout for testing :
# ==========================================================================================
Index <- createDataPartition(data_set1$factor_redData, p=.75,list = FALSE)

Train_set1 <- data_set1[Index,]
str(Train_set1)
Test_set1 <- data_set1[-Index,]
str(Test_set1)
# ==========================================================================================
# ==========================================================================================

# [C] DEVELOPING MODELS ON data_set1 ------------- -------------------------------------------------------------
# [C.1] ------------- extending rpart to cross validation with parameter tunning -------------------------------

ctrl <- trainControl(method = "cv", number = 10,selectionFunction = "oneSE")

grid <- expand.grid(cp=c(.010,.014,.018,.021,.039,.102,.718))
#modelLookup("rpart")

rpart_train <- train(factor_redData~., Train_set1,method="rpart",
                     trControl = ctrl,
                     tuneGrid = grid,
                     metric = "kappa")

rpart_train

rpart_train$results

rpart_train$finalModel # subject matter expert can help analyze TMEM98, DUSP4, ARHGDIB and CD63

dev.new()

rpart.plot(rpart_train$finalModel, digits = 2)

predict_rpart <- predict(rpart_train,Test_set1)
CrossTable(Test_set1$factor_redData,predict_rpart)

# [C] DEVELOPING MODELS ON data_set1 ------------- -------------------------------------------------------------
# [C.2] ------------- running knn on data_set1 --------------------------------- -------------------------------

library(ggplot2)
library(class)
library(gmodels)
library(caret)

set.seed(123)
ctrl <- trainControl(method = "repeatedcv", number = 10,repeats = 10,
                     selectionFunction = "oneSE")
# modelLookup("knn")
grid <- expand.grid(k=c(3,5,7,9,11))

knn_train <- train(factor_redData~., data=Train_set1,method ="knn",
                   metric = "Accuracy",
                   trControl=ctrl,
                   tuneGrid = grid
)

knn_train

knn_train$results

knn_train$finalModel 

predict_knn <- predict(knn_train,Test_set1)
CrossTable(Test_set1$factor_redData,predict_knn)

# [C] DEVELOPING MODELS ON data_set1 ------------- -------------------------------------------------------------
# [C.3] ------------- running naive on data_set1 --------------------------------- -----------------------------

library(MASS)
library(e1071)
library(klaR)


ctrl <- trainControl(method = "repeatedcv", number=10,repeats = 10,selectionFunction = "oneSE")

grid <- expand.grid(fL=1, usekernel=TRUE, adjust=1)

naive_train <- train(factor_redData~., data=Train_set1,method = "nb",
                     metric = "kappa",
                     trControl = ctrl,
                     tuneGrid = grid)
naive_train

naive_train$results

predict_naive <- predict(naive_train,Test_set1)
CrossTable(Test_set1$factor_redData,predict_naive)

#  Let's check on how data_set1 is clustered for k=2 ########################################################
###################################################### Let's check if data clustering improves on data_set1 :
#############################################################################################################

# standardize data
library(ggplot2)
library(caret)
rpartRed_z <- data.frame(lapply(reduced_rpart,scale)) # standardized numerical data of baseData
str(rpartRed_z)
summary(rpartRed_z[,1:5])
boxplot(rpartRed_z[,])

# RUNNING KMEANS :

# perform cluster analysis -------- ITERATION 2 on data_set1 ------------------------------------------------------
set.seed(76964058) #Set the seed for reproducibility
library(pracma)
library(cluster)

Kmeans_plus <- kmeans(rpartRed_z,2,iter.max = 100)
Kmeans_plus$size
table(factor_redData) 

dev.new()
clusplot(rpartRed_z, Kmeans_plus$cluster,color = TRUE,labels = 5) # clusters are overlapping but less than for data_set2 


###############################################################################################################
# ENSEMBLE MODELS ON data_set1 ------------- -----------------------------------------------------------------
# [C.4] ------------- Bagging on rpart --------------------------------- -------------------------------------
###############################################################################################################

library(ipred)
library(plyr)
library(caret)

set.seed(300)
ctrl <- trainControl(method = "cv", number = 10)
#grid <- expand.grid(cp=.01)

bagged_rpart <- train(factor_redData~., data=Train_set1,
                      method = "treebag",
                      #tuneGrid = grid,
                      trControl = ctrl)
#metric = "kappa")

bagged_rpart
bagged_rpart$results


predict_bagged <- predict(bagged_rpart,Test_set1)
CrossTable(Test_set1$factor_redData,predict_bagged)

###############################################################################################################
# ENSEMBLE MODELS ON data_set1 ------------- -----------------------------------------------------------------
# [C.5] ------------- random Forest on data_set1 ----------------------- -------------------------------------
###############################################################################################################

library(randomForest)

ctrl <- trainControl(method = "cv",number = 10)
grid_rf <- expand.grid(.mtry = c(3,6))

set.seed(321)
rf_rpart <- train(factor_redData~., data=data_set1,
                  method = "rf",
                  #metric = "kappa",
                  trControl = ctrl,
                  tuneGrid = grid_rf)

rf_rpart

rf_rpart$results


predict_rf <- predict(rf_rpart,Test_set1)
CrossTable(Test_set1$factor_redData,predict_rf)

##############################################################################################################
##############################################################################################################
# data_set2 extracted using ANOVA ---- dataset of 33 variables ,resulting in max separation ------------------
##############################################################################################################
##############################################################################################################

calcWithinGroupsVariance <- function(variable,groupvariable){
  # find out how many values the group variable can take
  groupvariable2 <- as.factor(groupvariable[[1]])
  levels <- levels(groupvariable2)
  numlevels <- length(levels)
  # get the mean and standard deviation for each group:
  numtotal <- 0
  denomtotal <- 0
  for (i in 1:numlevels)
  {
    leveli <- levels[i]
    levelidata <- variable[groupvariable==leveli,]
    levelilength <- length(levelidata)
    # get the standard deviation for group i:
    sdi <- sd(levelidata)
    numi <- (levelilength - 1)*(sdi * sdi)
    denomi <- levelilength
    numtotal <- numtotal + numi
    denomtotal <- denomtotal + denomi
  }
  # calculate the within-groups variance
  Vw <- numtotal / (denomtotal - numlevels)
  return(Vw)
}



my.calcBetween <- function(vec.dat, vec.grp){
  # vec.dat = data vector
  # vec.grp = group vector, same length as data vector
  groups <- unique(vec.grp)
  ngroups <- length(groups)
  grandmean <- mean(vec.dat)
  numtotal <- 0
  #	denomtotal <- 0
  for(i in groups){
    groupidata <- vec.dat[vec.grp == i]
    n.groupi <- length(groupidata)
    xbari <- mean(groupidata)
    sdi <- sd(groupidata)
    numi <- n.groupi * ((xbari - grandmean)^2)
    numtotal <- numtotal + numi
    # # 		denomi <- n.groupi
    # # 		denomtotal <- denomtotal + denomi
  }
  Vb <- numtotal/(ngroups-1)
  return(Vb)
}


separations <- NULL
calcSeparations <- function(variables,groupvariable){
  # find out how many variables we have
  variables <- as.data.frame(variables)
  numvariables <- length(variables)
  # find the variable names
  variablenames <- colnames(variables)
  # calculate the separation for each variable
  for (i in 1:numvariables)
  {
    variablei <- variables[i]
    variablename <- variablenames[i]
    Vw <- calcWithinGroupsVariance(variablei, groupvariable)
    #Vb <- calcBetweenGroupsVariance(variablei, groupvariable)
    Vb <- my.calcBetween(unlist(variablei), unlist(groupvariable))
    sep <- Vb/Vw
    separations <- rbind(c(variablename,sep),separations)
    #print(paste("variable",variablename,"Vw=",Vw,"Vb=",Vb,"separation=",sep))
  }
  return(separations)
}	

# [D] CULLING OUT data_set2 from baseData using ANOVA -------------------------------------------------------------
############################################ ----------------------------------------------------------------------

# dataframe for running calcSeparations
str(melanoma)
factor_redData <- melanoma$factor_vector # 4513 X 1
num_redData <- melanoma[-c(1:2)] # 4513 X 5564
str(num_redData) 

baseData <- cbind(factor_redData,num_redData) # df we are using for further analysis

data_z <- data.frame(lapply(num_redData,scale)) # standardized numerical data of baseData
pca_df <- cbind(factor_redData,data_z)
str(pca_df)

separations_pca <- calcSeparations(pca_df[2:5565],pca_df[1])

sorted_sep <- sort(separations_pca,decreasing = TRUE)

var_MaxSep <- sorted_sep[1:33]
var_MaxSep
Var_rpart <- names(Var_rpart)
Var_rpart

length(which(var_MaxSep==Var_rpart)) # 0 

# create new data frame including var_MaxSep
reduced_MaxSep <- baseData[names(baseData) %in% c(var_MaxSep)] # df with 33 variables
str(reduced_MaxSep)

#########========================================================================
#########========================================================================
data_set2 <- cbind(factor_redData,reduced_MaxSep)
#########========================================================================
#########========================================================================
str(data_set2)


# [E] DEVELOPING MODELS ON data_set2 ; not as good as on data_set1 so not perfecting these models  -------------
# [E.1] ------------- rpart on data_set2;cross validation with parameter tunning -------------------------------

library(ggplot2)
library(caret)
library(rpart.plot)
library(gmodels)

ctrl <- trainControl(method = "cv", number = 10,selectionFunction = "oneSE")

grid <- expand.grid(cp=c(.010,.014,.018,.021,.039,.102,.718))
#modelLookup("rpart")

rpart2_train <- train(factor_redData~., data_set2,method="rpart",
                      trControl = ctrl,
                      tuneGrid = grid,
                      metric = "kappa")

rpart2_train

rpart2_train$results

rpart2_train$finalModel # subject matter expert can help analyze ZNF664,ZNF714,ZNF706,ZZZ3,ZNHIT1,ZNF850,ZRANB2 and ZYX

dev.new()

rpart.plot(rpart2_train$finalModel, digits = 2)

predict_rpart2 <- predict(rpart2_train,data_set2)
CrossTable(factor_redData,predict_rpart2)


# [E] DEVELOPING MODELS ON data_set2 ; not as good as on data_set1 so not perfecting these models  -------------
# [E.2] ------------- knn on data_set2;cross validation with parameter tunning ---------------------------------

library(class)


set.seed(123)
ctrl <- trainControl(method = "repeatedcv", number = 10,repeats = 10,
                     selectionFunction = "oneSE")
# modelLookup("knn")
grid <- expand.grid(k=c(3,5,7,9,11))

knn_train <- train(factor_redData~., data=data_set2,method ="knn",
                   metric = "kappa",
                   trControl=ctrl,
                   tuneGrid = grid
)

knn_train

knn_train$results

knn_train$finalModel

predict_knn <- predict(knn_train,data_set2)
CrossTable(factor_redData,predict_knn)

# [E] DEVELOPING MODELS ON data_set2 ; not as good as on data_set1 so not perfecting these models  -------------
# [E.3] ------------- naive on data_set2;cross validation with parameter tunning -------------------------------

library(e1071)
library(klaR)
library(MASS)

ctrl <- trainControl(method = "repeatedcv", number=10,repeats = 10,selectionFunction = "oneSE")

grid <- expand.grid(fL=1, usekernel=TRUE, adjust=1)

naive_train <- train(factor_redData~., data=data_set2,method = "nb",
                     metric = "kappa",
                     trControl = ctrl,
                     tuneGrid = grid)
naive_train

naive_train$results

predict_naive <- predict(naive_train,data_set2)
CrossTable(factor_redData,predict_naive)

#  Let's check on how data_set2 is clustered for k=2 ########################################################
########################################################  as expected clustering is worse than on data_set1 :
#############################################################################################################

# standardize data_reduced_MaxSep
library(ggplot2)
library(caret)
MaxSepRed_z <- data.frame(lapply(reduced_MaxSep,scale)) # standardized numerical data of baseData
str(MaxSepRed_z)
summary(MaxSepRed_z[,1:5])
boxplot(MaxSepRed_z[,])

# RUNNING KMEANS :

# perform cluster analysis -------- ITERATION 2 on data_set2 ------------------------------------------------
set.seed(76964059) #Set the seed for reproducibility
library(pracma)
library(cluster)

Kmeans_plus2 <- kmeans(MaxSepRed_z,2,iter.max = 100)
Kmeans_plus2$size
table(factor_redData) # grown worse

dev.new()
clusplot(MaxSepRed_z, Kmeans_plus2$cluster,color = TRUE,labels = 5) # data_set2 isn't optimum for our analysis 

"We are presenting data_set1 variables as the Gene Markers that can be instrumental in identifying malignant cells"
################## END OF THE CODE ###################
##########################################
################################
######################
##############
#######
###
#













