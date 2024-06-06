
# Code used to support dissertation analyses
#
#  Challenges in Exposure Characterization for Healthcare-Associated Infections: 
#  Applications in the Clinical Epidemiology of Clostridioides difficile Infection
#
# Jessica Lynn Webster, PhD, MPH
# Drexel University, Dornsife School of Public Health
# May 2024


##-------------------------------------##
# Aim 1                                 #
# 1.1. factor analysis of mixed data    #
# 1.2. logistic regression              #
# 1.3. LASSO regression                 #
#                                       #
# Updated 07.10.2023                    #
##-------------------------------------##

#### 1.1. Factor analysis of mixed data ####
# Calling in libraries
library(tidyverse)
library(glmnet)
library(corrplot)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)

# Pulling in FULL wide dataset
wide <- read.csv("HUH_clean_wide_full.csv", header=T, na.strings = c("","NA"))

#### Data Preparation ####
#define response variable
y <- wide$culture_cdiff

# only exposure variables (binary and continuous)
exp_vars_b <- c("currentabx_any","aminoglycosides_any","monobactam_any","penicillin_any","rifamycins_any",
                "sulfonamides_any","b_lactamase_inhib_any","carbapenem_any","cephalosporin_any","fluoroquinolone_any",
                "glycolipopeptides_lipopeptides_any","macrolides_lincosamide_any","misc_any")
exp_vars_c <- c("dose_gl","currentabx_days","prop_abx_days","currentabx_total","classes")

# only covariates (binary and continuous)
co_vars_b <- c("sex","race_eth","insur_cat","refer_bin","icu","culture_other_bin","surg","currentproton","currentsteroid",
               "currentchemo","prior_hosp","priorabx_bin","prior_proton","prior_steroid","prior_chemo")
co_vars_c <- c("age","bmi","los","cci")

# interactions
## abx and ppi, current and prior
wide$curr_abx_ppi <- wide$currentabx_any*wide$currentproton
wide$prior_abx_ppi <- wide$priorabx_bin*wide$prior_proton

int_vars_b <- c("curr_abx_ppi","prior_abx_ppi")

## define matrix of predictor variables: exposure variables
# only exposure variables (binary and continuous)
exp_vars <- wide[, c(exp_vars_b,exp_vars_c)]
#scaling continuous exposure vars
exp_vars[, c(exp_vars_c)] <- scale(exp_vars[, c(exp_vars_c)])

# only covariates (binary and continuous)
co_vars <- wide[,c(co_vars_b, co_vars_c)]
co_vars[,c(co_vars_c)] <- scale(co_vars[,c(co_vars_c)])#scaling continuous vars

# only interaction variables
int_vars <- wide[, int_vars_b]

# converting binary and categorical variables to factors
exp_vars <- exp_vars %>% mutate(across(exp_vars_b, as.factor))
co_vars <- co_vars %>% mutate(across(co_vars_b, as.factor))
int_vars <- int_vars %>% mutate(across(int_vars_b, as.factor))

# combining exposure vars and covariates
x_df <- cbind(exp_vars,co_vars,int_vars)

## X with all exposure vars AND covariates, imputing missing ##
x <- makeX(x_df, na.impute = T)

#subsetting data to only outcome and exposure vars
wide_sub <- cbind(exp_vars,y)
wide_sub$y <- as.factor(wide_sub$y)
levels(wide_sub$y) <- c("No","Yes")

#subsetting data to outcome, exposure vars, and covariates
wide_sub_all <- as.data.frame(cbind(x_df,y))

#--- Unsupervised feature selection ---#
#### factor analysis of mixed data ####

# correlation matrices
model.matrix(~0+., data=wide_sub) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag=T, type="lower", lab=TRUE, lab_size=2, 
             colors = c("#00AFBB", "white", "#FC4E07"))

corr <- model.matrix(~0+., data=wide_sub) %>% 
  cor(use="pairwise.complete.obs")

# factor analysis of mixed data (FAMD) using FactoMineR
# only with exposure variables
res.famd <- FAMD(wide_sub, 
                 ncp=6,
                 sup.var = 19,  ## Set the outcome variable "y" as a supplementary variable
                 graph = FALSE)

## Inspect principal components
get_eigenvalue(res.famd)
summary(res.famd)

## Set figure size
options(repr.plot.width = 14, repr.plot.height = 12)

summary(res.famd)
fviz_eig(res.famd, addlabels = TRUE) # skree plot
fviz_famd_var(res.famd, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE) #contributions
fviz_famd_var(res.famd, col.var = "cos2",gradient.cols=c("#5c3391"),repel = TRUE) #cos2

# contribution plot for all dimensions
var <- get_famd_var(res.famd)
corrplot(var$contrib, is.corr=FALSE, tl.cex = 1)

# https://f0nzie.github.io/machine_learning_compilation/detailed-study-of-principal-component-analysis.html
# var$coord: coordinates of variables to create a scatter plot
# var$cos2: represents the quality of representation for variables on the factor map. 
#   It's calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# var$contrib: contains the contributions (in percentage) of the variables to the principal components. 
#   The contribution of a variable (var) to a given principal component is (in percentage): 
#   (var.cos2 * 100) / (total cos2 of the component).

# Contribution to all six dimensions
# dashed line is the expected average value if all contributions were uniform
fviz_contrib(res.famd, "var", axes = 1:6)
fviz_cos2(res.famd, "var", axes = 1:2)


#### 1.2. Logistic regression ####
# Calling in libraries
library(caret)

# subsetting wide data
vars <- c("culture_cdiff","currentabx_any","currentabx_total","classes","dose_gl","currentabx_days","prop_abx_days"
          ,"aminoglycosides_any","b_lactamase_inhib_any","carbapenem_any","cephalosporin_any"
          ,"fluoroquinolone_any","glycolipopeptides_lipopeptides_any","macrolides_lincosamide_any","monobactam_any"
          ,"penicillin_any","rifamycins_any","sulfonamides_any","misc_any"
          ,"insur_cat","los","cci","icu","culture_other_bin","currentchemo","prior_hosp","priorabx_bin" 
          ,"prior_proton","prior_steroid","currentproton")

wide_vars <- wide[,vars]

### Logistic regression models for each exposure measurement ###
# all adjusted for insurance, LOS, CCI, ICU, prior abx, prior hosp, non-c.diff infection
# class variables adjusted for all classes

## Model with EVERYTHING ! ##
m <- glm(culture_cdiff ~ currentabx_any + currentabx_total + classes + dose_gl + currentabx_days + prop_abx_days
         + aminoglycosides_any+b_lactamase_inhib_any+carbapenem_any+cephalosporin_any
         + fluoroquinolone_any+glycolipopeptides_lipopeptides_any+macrolides_lincosamide_any+monobactam_any
         + penicillin_any+rifamycins_any+sulfonamides_any+misc_any
         + insur_cat + los + cci + icu + culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
         + prior_proton + prior_steroid + (currentabx_any*currentproton)
         + (priorabx_bin*prior_proton), data=wide, family="binomial")
summary(m);confint(m)


## Model 1: Ever received abx ##
m1 <- glm(culture_cdiff ~ currentabx_any 
          + insur_cat + los + cci + icu + culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
          + prior_proton + prior_steroid + currentproton + (currentabx_any*currentproton), data=wide_vars, family="binomial")
summary(m1);confint(m1)
exp(m1$coefficients)

## Model 2: cumulative number of unique abx treatments ##
m2 <- glm(culture_cdiff ~ currentabx_total 
          + insur_cat + los + cci + icu + culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
          + prior_proton + prior_steroid + (currentabx_any*currentproton), data=wide, family="binomial")
summary(m2);confint(m2)

## Model 3a: class ##
m3 <- glm(culture_cdiff ~ aminoglycosides_any+b_lactamase_inhib_any+carbapenem_any+cephalosporin_any
          +fluoroquinolone_any+glycolipopeptides_lipopeptides_any+macrolides_lincosamide_any+monobactam_any
          +penicillin_any+rifamycins_any+sulfonamides_any+misc_any
          + insur_cat + los + cci + icu + culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
          + prior_proton + prior_steroid + (currentabx_any*currentproton), data=wide, family="binomial")
summary(m3);confint(m3)

## Model 4: cumulative number of unique classes ##
m4 <- glm(culture_cdiff ~ classes 
          + insur_cat + los + cci + icu + culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
          + prior_proton + prior_steroid + (currentabx_any*currentproton), data=wide, family="binomial")
summary(m4);confint(m4)

## Model 5: dose ##
m5 <- glm(culture_cdiff ~ dose_gl 
          + insur_cat + los + cci + icu + culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
          + prior_proton + prior_steroid + (currentabx_any*currentproton), data=wide, family="binomial")
summary(m5);confint(m5)

## Model 6: cumulative number of days on abx ##
m6 <- glm(culture_cdiff ~ currentabx_days
          + insur_cat + los + cci + icu + culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
          + prior_proton + prior_steroid + (currentabx_any*currentproton), data=wide, family="binomial")
summary(m6);confint(m6)

## Model 7: proportion of stay on abx ##
# removed LOS since exposure variable is calculated based on LOS
m7 <- glm(culture_cdiff ~ prop_abx_days  
          + insur_cat + los + cci + #icu + 
            culture_other_bin + currentchemo + prior_hosp + priorabx_bin 
          + prior_proton + prior_steroid + (currentabx_any*currentproton), data=wide, family="binomial")
summary(m7);confint(m7)



#### 1.3. LASSO regression ####
library(caret)

#### with GLMNET package ####
set.seed(1234)

# cross validation to find lambda
glmnet.model<- cv.glmnet(x=x,y=y,
                         family = "binomial",
                         alpha=1) #alpha=1 is lasso

plot(glmnet.model)
l.min <- glmnet.model$lambda.min
l.1se <-glmnet.model$lambda.1se

# running lasso model with min lambda
lasso.model_l.min <- glmnet( x=x,y=y, 
                             family = "binomial", 
                             alpha=1, 
                             lambda = l.min)
# sparse matrix
round(lasso.model_l.min$beta,4)

# assessing lasso model
assess.glmnet(lasso.model_l.min, # in this case, we are evaluating the model
              newx = x, # in the same data used to fit the model
              newy = y)  # so newx=X and newY=Y

plot(roc.glmnet(lasso.model_l.min,newx = x, newy = y), type="l") #produces the ROC plot

# importance
coefList <- coef(lasso.model_l.min, s='lambda.1se')
coefList <- data.frame(coefList@Dimnames[[1]][coefList@i+1],coefList@x)
names(coefList) <- c('var','val')

coefList %>%
  arrange(-abs(val)) %>%
  print(.)

varImp <- function(object, lambda = NULL, ...) {
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out)
  } else out <- data.frame(Overall = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out <- out/max(out)
  out[order(out$Overall, decreasing = TRUE),,drop=FALSE]
}

lass0 <- round(varImp(lasso.model_l.min, lambda = lasso.model_l.min$lambda.1se),2);lass0$vars <- row.names(lass0) 
lass0_vars <- lass0[lass0$Overall>0,]
lass0_vars$vars

# confusion matrix
itrain <- sample(c(TRUE,FALSE), nrow(x),rep=TRUE, prob=c(.67, .33))
cfit <- cv.glmnet(x[itrain, ], y[itrain], family = "binomial", alpha=1)
cnf <- confusion.glmnet(cfit, newx = x[-itrain, ], newy = y[-itrain])
print(cnf)
summary(cnf)
# assessing model
assess.glmnet(lasso.model_l.min,newx = x, newy = y )

## running lasso model with 1 se lambda
lasso.model_l.1se <- glmnet( x=x,y=y, 
                             family = "binomial", 
                             alpha=1, 
                             lambda = l.1se)
# sparse matrix
round(lasso.model_l.1se$beta,4)

# assessing lasso model
assess.glmnet(lasso.model_l.1se,
              newx = x, 
              newy = y)

plot(roc.glmnet(lasso.model_l.1se,newx = x, newy = y), type="l")

# importance
coefList <- coef(lasso.model_l.1se, s='lambda.1se')
coefList <- data.frame(coefList@Dimnames[[1]][coefList@i+1],coefList@x)
names(coefList) <- c('var','val')

coefList %>%
  arrange(-abs(val)) %>%
  print(.)

varImp <- function(object, lambda = NULL, ...) {
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out)
  } else out <- data.frame(Overall = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out <- out/max(out)
  out[order(out$Overall, decreasing = TRUE),,drop=FALSE]
}

lass0 <- round(varImp(lasso.model_l.1se, lambda = lasso.model_l.1se$lambda.1se), 2)
lass0$vars <- row.names(lass0) 
lass0_vars <- lass0[lass0$Overall>0,]
lass0_vars$vars

# confusion matrix
itrain <- sample(c(TRUE,FALSE), nrow(x),rep=TRUE, prob=c(.67, .33))
cfit <- cv.glmnet(x[itrain, ], y[itrain], family = "binomial", alpha=1)
cnf <- confusion.glmnet(cfit, newx = x[-itrain, ], newy = y[-itrain])
print(cnf)
summary(cnf)
# assessing model
assess.glmnet(lasso.model_l.1se,newx = x,newy = y )


## running a separate predictive model with the selected variables
# creating list of variables to include
lasso_vars_l.min <- wide_sub_all[,c("icu","priorabx_bin","culture_other_bin","currentchemo",         
                                    "insur_cat","rifamycins_any","carbapenem_any","currentabx_any","misc_any","surg",
                                    "prior_steroid","fluoroquinolone_any","currentsteroid","sex", 
                                    "b_lactamase_inhib_any","prior_proton","sulfonamides_any","prop_abx_days",
                                    "cephalosporin_any","race_eth","monobactam_any","curr_abx_ppi",
                                    "y")]
lasso_vars_l.1se <- wide_sub_all[,c("priorabx_bin","currentabx_any","insur_cat","culture_other_bin","icu","misc_any",
                                    "prop_abx_days","carbapenem_any","prior_abx_ppi","cephalosporin_any",
                                    "currentchemo","rifamycins_any","prior_proton",
                                    "y")]
lasso_comp_l.min <- lasso_vars_l.min[complete.cases(lasso_vars_l.min),]
lasso_comp_l.1se <- lasso_vars_l.1se[complete.cases(lasso_vars_l.1se),]

## testing predictive performance of model with final variables selected ** LMIN
set.seed(1234)
# creating index
index <- createDataPartition(lasso_comp_l.min$y, p=.8, list=FALSE, times=1)
# Create test and train data frames
train_df <- lasso_comp_l.min[index,]
test_df <- lasso_comp_l.min[-index,]

# Verify number of rows (cases) in each data frame
nrow(train_df);nrow(test_df)

# Re-label values of outcome variable for train_df
train_df$y[train_df$y==0] <- "control"
train_df$y[train_df$y==1] <- "case"

# Re-label values of outcome variable for test_df
test_df$y[test_df$y==0] <- "control"
test_df$y[test_df$y==1] <- "case"

# Convert outcome variable to factor for each data frame
train_df$y <- as.factor(train_df$y);test_df$y <- as.factor(test_df$y)

# Specify type of training method used and the number of folds
ctrlspecs <- trainControl(method="repeatedcv", 
                          number=10, 
                          savePredictions="all",
                          repeats = 10,
                          classProbs=TRUE)

# Specify logistic regression model to be estimated using training data
# and k-fold cross-validation process
model_l.min <- train(y ~ ., data=train_df, 
                     method="glm", 
                     family=binomial(link="logit"),
                     trControl=ctrlspecs)

# Print information about model
print(model_l.min)
# Print results of final model estimated using training data
summary(model_l.min)
# Predict outcome using model from training data based on testing data
predictions <- predict(model_l.min, newdata=test_df)
# Create confusion matrix to assess model fit/performance on test data
confusionMatrix(data=predictions, test_df$y)


## testing predictive performance of model with final variables selected ** L1SE
set.seed(1234)
# creating index
index <- createDataPartition(lasso_comp_l.1se$y, p=.8, list=FALSE, times=1)
# Create test and train data frames
train_df <- lasso_comp_l.1se[index,]
test_df <- lasso_comp_l.1se[-index,]

# Verify number of rows (cases) in each data frame
nrow(train_df);nrow(test_df)

# Re-label values of outcome variable for train_df
train_df$y[train_df$y==0] <- "control"
train_df$y[train_df$y==1] <- "case"

# Re-label values of outcome variable for test_df
test_df$y[test_df$y==0] <- "control"
test_df$y[test_df$y==1] <- "case"

# Convert outcome variable to factor for each data frame
train_df$y <- as.factor(train_df$y);test_df$y <- as.factor(test_df$y)

# Specify type of training method used and the number of folds
ctrlspecs <- trainControl(method="repeatedcv", 
                          number=10, 
                          savePredictions="all",
                          repeats = 10,
                          classProbs=TRUE)

# Specify logistic regression model to be estimated using training data
# and k-fold cross-validation process
model_l.1se <- train(y ~ ., data=train_df, 
                     method="glm", 
                     family=binomial(link="logit"),
                     trControl=ctrlspecs)

# Print information about model
print(model_l.1se)
# Print results of final model estimated using training data
summary(model_l.1se)
# Predict outcome using model from training data based on testing data
predictions <- predict(model_l.1se, newdata=test_df)
# Create confusion matrix to assess model fit/performance on test data
confusionMatrix(data=predictions, test_df$y)





