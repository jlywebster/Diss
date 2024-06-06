
# Code used to support dissertation analyses
#
#  Challenges in Exposure Characterization for Healthcare-Associated Infections: 
#  Applications in the Clinical Epidemiology of Clostridioides difficile Infection
#
# Jessica Lynn Webster, PhD, MPH
# Drexel University, Dornsife School of Public Health
# May 2024


##-------------------------------------##
# Aim 2                                 #
# 2.1. inverse probability weighting    #
# 2.2. marginal structural models       #
#                                       #
# Updated 10.01.2023                    #
##-------------------------------------##

# Calling in libraries
library(tidyverse)
library(twang)
library(cobalt)
library(survey)
library(tableone)
library(WeightIt)
library(geepack)
library(glmnet)

# Pulling in full long dataset
long_a <- read.csv("CDI_long_analytic.csv", header=T, na.strings = c("","NA"))

#### 2.1. inverse probability of exposure weighting ####

## calculating denominator
# exposure = any abx (currentabx_bin)
# creating cumulative variables to represent days in hospital stay and abx exposure history up until the previous day
long_a$hosp_day_yest <- long_a$hosp_day - 1
long_a$cum_abx_days_yest <- ifelse(long_a$cum_abx_days==0, 0, long_a$cum_abx_days - 1)
long_a$cum_prop_abx_yest <- ifelse(long_a$cum_abx_days_yest==0 & long_a$hosp_day_yest==0, 0, long_a$cum_abx_days_yest/long_a$hosp_day_yest)

# creating a cumulative number of days on current abx
long_a <- long_a %>%
  group_by(study_redcap_id, grp = cumsum(currentabx_bin == 0)) %>%
  mutate(cum_course_days_yest = (row_number() - 1)*currentabx_bin) %>%
  ungroup()

long_a <- as.data.frame(long_a)

#### Using LASSO to select variables that predict exposure ####
# define response variable
y <- long_a$currentabx_bin

# covariates (binary and continuous)
co_vars_b <- c("demog_gender","race_eth","insur_cat","refer_bin","course_icu","culture_other_bin","course_surg",
               "currentproton","currentsteroid",
               "currentchemo","hosp_any","priorabx_bin","prior_proton","prior_steroid","prior_chemo")
co_vars_c <- c("demog_age","demog_bmi","cci_total","hosp_day_yest","cum_abx_days_yest")

## define matrix of predictor variables: covariates
co_vars <- long_a[,c(co_vars_b, co_vars_c)]
co_vars[,c(co_vars_c)] <- scale(co_vars[,c(co_vars_c)])#scaling continuous vars

# converting binary and categorical variables to factors
co_vars <- co_vars %>% mutate(across(co_vars_b, as.factor))

## covariates, imputing missing ##
x <- makeX(co_vars, na.impute = T)

### with GLMNET package ###
set.seed(1234)

# cross validation to find lambda
glmnet.model<- cv.glmnet(x=x,y=y,
                         family = "binomial",
                         alpha=1)

plot(glmnet.model)
l.min <- glmnet.model$lambda.min
l.1se <-glmnet.model$lambda.1se

# running lasso model with min lambda
lasso.model <- glmnet(x=x,y=y,
                      family = "binomial",
                      alpha=1,
                      lambda = l.1se)
# sparse matrix
round(lasso.model$beta,4)

# assessing lasso model
assess.glmnet(lasso.model,
              newx = x,
              newy = y)

plot(roc.glmnet(lasso.model,newx = x, newy = y), type="l")

# importance
coefList <- coef(lasso.model, s='lambda.1se')
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

lass0 <- round(varImp(lasso.model, lambda = lasso.model$lambda.1se),2);lass0$vars <- row.names(lass0)
lass0_vars <- lass0[lass0$Overall>0,]
lass0_vars$vars

# LASSO model selected
# demog_gender,race_eth,course_surg,refer_bin,culture_other_bin,cci_total,course_icu,
# currentproton, currentsteroid, currentchemo
# hosp_any, priorabx_bin, 
# cum_abx_days_yest, hosp_day_yest

## some cleaning 
long_a$cci_cat <- ifelse(is.na(long_a$cci_total), "Missing", 
                         ifelse(long_a$cci_total<=2, "Mild",
                                ifelse(long_a$cci_total==3 | long_a$cci_total==4, "Mod", "Severe")))
long_a$bmi_cat <- ifelse(is.na(long_a$demog_bmi), "Missing", 
                         ifelse(long_a$demog_bmi<=18.5, "UW",
                                ifelse(long_a$demog_bmi<25, "H",
                                       ifelse(long_a$demog_bmi<30, "OW",
                                              ifelse(long_a$demog_bmi<40, "O","SO")))))

# converting binary/categorical variables to factor
long_a$demog_gender_f <-  as.factor(long_a$demog_gender)
long_a$race_eth_f <-  as.factor(long_a$race_eth)
long_a$insur_cat_f <-  as.factor(long_a$insur_cat)
long_a$insur_cat_f <-  as.factor(long_a$insur_cat)
long_a$course_surg_f <-  as.factor(long_a$course_surg)
long_a$course_icu_f <-  as.factor(long_a$course_icu)
long_a$refer_bin_f <-  as.factor(ifelse(is.na(long_a$refer_bin), "Missing", long_a$refer_bin))
long_a$culture_other_bin_f <-  as.factor(long_a$culture_other_bin)
long_a$currentproton_f <-  as.factor(long_a$currentproton)
long_a$currentsteroid_f <-  as.factor(long_a$currentsteroid)
long_a$currentchemo_f <-  as.factor(long_a$currentchemo)
long_a$hosp_any_f <-  as.factor(ifelse(is.na(long_a$hosp_any), "Missing", long_a$hosp_any))
long_a$prior_proton_f <-  as.factor(long_a$prior_proton)
long_a$prior_steroid_f <-  as.factor(long_a$prior_steroid)
long_a$prior_chemo_f <-  as.factor(long_a$prior_chemo)
long_a$priorabx_bin_f <-  as.factor(long_a$priorabx_bin)
long_a$bmi_cat_f <-  as.factor(long_a$bmi_cat)
long_a$cci_cat_f <-  as.factor(long_a$cci_cat)


#### Calculating IPW using variables selected by LASSO + treatment history ####
## using three methods: regression, CART, weight it

## variables included:
# demog_gender,race_eth,course_surg,refer_bin,culture_other_bin,cci_total,course_icu,
# currentproton, currentsteroid, currentchemo
# hosp_any, priorabx_bin, 
# cum_abx_days_yest, hosp_day_yest

### Denominator
## Regression
ps.model <- glm(currentabx_bin ~ cum_abx_days_yest+hosp_day_yest+
                  demog_gender_f+race_eth_f+course_surg_f+refer_bin_f+course_icu_f+
                  culture_other_bin_f+cci_cat_f+currentproton_f+currentsteroid_f+currentchemo_f+
                  hosp_any_f+priorabx_bin_f,
                data=long_a, family=binomial(link="logit"))
long_a$PS <- predict(ps.model, long_a, type="response")

## CART
ps.model.CART <- ps(currentabx_bin ~ cum_abx_days_yest+hosp_day_yest+
                      demog_gender_f+race_eth_f+course_surg_f+refer_bin_f+course_icu_f+
                      culture_other_bin_f+cci_cat_f+currentproton_f+currentsteroid_f+currentchemo_f+
                      hosp_any_f+priorabx_bin_f,
                    data=long_a)
long_a$PS_CART <- ps.model.CART$ps[, 1]

### Numerator
# treatment history and baseline covariates (stabilized)
## Regression
stab.model <- glm(currentabx_bin ~ cum_abx_days_yest+hosp_day_yest,
                  data=long_a, family=binomial(link="logit"))
long_a$stab_PS <- predict(stab.model, long_a, type="response")

## CART
set.seed(1234)
stab.model.CART <- ps(currentabx_bin ~ cum_abx_days_yest+hosp_day_yest,
                      data=long_a)
long_a$stab_PS_CART <- stab.model.CART$ps[,1]


### Calculating stabilized inverse probability weights
## Regression
long_a$PS_weights_reg <- ifelse(long_a$currentabx_bin==1, long_a$stab_PS/long_a$PS, long_a$stab_PS/(1-long_a$PS))
## CART
long_a$PS_weights_cart <- ifelse(long_a$currentabx_bin==1, long_a$stab_PS_CART/long_a$PS_CART, long_a$stab_PS_CART/(1-long_a$PS_CART))

# visualizing
boxplot(long_a$PS_weights_reg~long_a$currentabx_bin)
boxplot(long_a$PS_weights_cart~long_a$currentabx_bin)

# TRIMMING the stabilized IPW weights
long_a$PS_weights_reg <- ifelse(long_a$PS_weights_reg < quantile(long_a$PS_weights_reg, 0.01, na.rm=T), quantile(long_a$PS_weights_reg, 0.01, na.rm=T), long_a$PS_weights_reg)
long_a$PS_weights_reg <- ifelse(long_a$PS_weights_reg > quantile(long_a$PS_weights_reg, 0.99, na.rm=T), quantile(long_a$PS_weights_reg, 0.99, na.rm=T), long_a$PS_weights_reg)

long_a$PS_weights_cart <- ifelse(long_a$PS_weights_cart < quantile(long_a$PS_weights_cart, 0.01, na.rm=T), quantile(long_a$PS_weights_cart, 0.01, na.rm=T), long_a$PS_weights_cart)
long_a$PS_weights_cart <- ifelse(long_a$PS_weights_cart > quantile(long_a$PS_weights_cart, 0.99, na.rm=T), quantile(long_a$PS_weights_cart, 0.99, na.rm=T), long_a$PS_weights_cart)

# visualizing again
boxplot(long_a$PS_weights_reg~long_a$currentabx_bin)
boxplot(long_a$PS_weights_cart~long_a$currentabx_bin)

## Weight diagnostics
# variables included
variables <- c("cum_abx_days_yest","hosp_day_yest",
               "demog_gender_f","race_eth_f","course_surg_f","refer_bin_f",
               "course_icu_f","culture_other_bin_f","cci_cat_f",
               "currentproton_f","currentsteroid_f","currentchemo_f",
               "hosp_any_f","priorabx_bin_f")
covs <- subset(long_a, select = variables)

# balance table
balance_reg <- bal.tab(covs, treat = long_a$currentabx_bin, weights = long_a$PS_weights_reg, thresholds = c(m = .1, v = 2));balance_reg
balance_cart <- bal.tab(covs, treat = long_a$currentabx_bin, weights = long_a$PS_weights_cart, thresholds = c(m = .1, v = 2));balance_cart

## Creating weights using the "weight it" package
iptw <- weightit(
  currentabx_bin ~ cum_abx_days_yest+hosp_day_yest+
    demog_gender_f+race_eth_f+course_surg_f+refer_bin_f+course_icu_f+
    culture_other_bin_f+cci_cat_f+currentproton_f+currentsteroid_f+currentchemo_f+
    hosp_any_f+priorabx_bin_f,
  data = long_a, estimand = "ATE", method = "ps", stabilize = T, include.obj = T)

# trimming
iptw.t <- trim(iptw, at = .99)
set.cobalt.options(binary = "std")

# balance table
table_iptw <- bal.tab(iptw.t, stats = c("m"), thresholds = c(m = 0.1));table_iptw

# creating a weight variable
long_a$weightit_weights <- iptw.t$weights

## love.plot(w.out1.t)
new.names <- c(hosp_day_yest = "Cumulative length of stay",
               cum_abx_days_yest = "Cumulative number of days on antibiotics",
               demog_gender_f_2 = "Gender (F)",
               race_eth_f_Hisp = "Race: Hispanic",
               race_eth_f_NHA = "Race: Non-Hispanic Asian",
               race_eth_f_NHB = "Race: Non-Hispanic Black",
               race_eth_f_NHW = "Race: Non-Hispanic White",
               'race_eth_f_Other/Unknown' = "Race: Other/Unknown",
               insur_cat_f_Missing = "Insurance: Missing",
               insur_cat_f_Other = "Insurance: Other",
               insur_cat_f_Private = "Insurance: Private",
               insur_cat_f_Public = "Insurance: Public",
               course_surg_f_0 = "Surgery during stay: No",
               course_surg_f_1 = "Surgery during stay: Yes",
               course_surg_f_Missing = "Surgery during stay: Missing",
               priorabx = "Number of prior antibiotic treatments",
               priorabx_bin_f = "Prior antibiotic treatment: Yes",
               prior_proton_f = "Prior PPI administration",
               prior_chemo_f = "Prior chemotherapy",
               prior_steroid_f = "Prior steroid administration",
               currentsteroid_f = "Current steroid administration",
               culture_other_bin_f = "Non-C.difficile infection",
               currentproton_f = "Current PPI administration",
               currentchemo_f = "Current chemotherapy",
               cci_cat_f_Mild = "CCI: Mild",
               cci_cat_f_Missing = "CCI: Missing",
               cci_cat_f_Severe = "CCI: Severe",
               cci_cat_f_Mod = "CCI: Moderate",
               course_icu_f_1 = "ICU: Yes",
               course_icu_f_0 = "ICU: No",
               course_icu_f_Missing = "ICU: Missing",
               hosp_any_f_1 = "Prior hospitalization: Yes",
               hosp_any_f_0 = "Prior hospitalization: No",
               hosp_any_f_Missing = "Prior hospitalization: Missing",
               bmi_cat_f_H = "BMI: Healthy",
               bmi_cat_f_UW = "BMI: Underweight",
               bmi_cat_f_OW = "BMI: Overweight",
               bmi_cat_f_O = "BMI: Obese",
               bmi_cat_f_SO = "BMI: Severely obese",
               bmi_cat_f_Missing = "BMI: Missing",
               refer_bin_f_0 = "Referral: Other healthcare facility",
               refer_bin_f_1 = "Referral: Home",
               refer_bin_f_Missing = "Referral: Missing"
)

love.plot(currentabx_bin ~ cum_abx_days_yest+hosp_day_yest+
            demog_gender_f+race_eth_f+course_surg_f+refer_bin_f+course_icu_f+
            culture_other_bin_f+cci_cat_f+currentproton_f+currentsteroid_f+currentchemo_f+
            hosp_any_f+priorabx_bin_f, 
          data = long_a, weights = long_a$PS_weights_reg,
          drop.distance = TRUE, 
          var.order = "unadjusted",
          abs = TRUE,
          line = TRUE, 
          thresholds = c(m = .1),
          var.names = new.names,
          colors = c("#d1482a", "#2896b5"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted"),
          limits = c(0, .5),
          position = c(.75, .25)) +
  theme(legend.box.background = element_rect(), 
        legend.box.margin = margin(1, 1, 1, 1))


## Weighted data
longsvy_reg <- svydesign(ids = ~ 1, data = long_a, weights = ~ PS_weights_reg)
longsvy_cart <- svydesign(ids = ~ 1, data = long_a, weights = ~ PS_weights_cart)
longsvy_weightit <- svydesign(ids = ~ 1, data = long_a, weights = ~ weightit_weights)

## Construct a table
tabWeighted_reg <- svyCreateTableOne(vars = variables, strata = "currentabx_bin", data = longsvy_reg, test = FALSE)
tabWeighted_cart <- svyCreateTableOne(vars = variables, strata = "currentabx_bin", data = longsvy_cart, test = FALSE)
tabWeighted_weightit <- svyCreateTableOne(vars = variables, strata = "currentabx_bin", data = longsvy_weightit, test = FALSE)

## Show table with SMD
tabreg <- print(tabWeighted_reg, smd = TRUE)
tabcart <- print(tabWeighted_cart, smd = TRUE)
tabwi <- print(tabWeighted_weightit, smd = TRUE)

# exporting to CSV
write.csv(tabreg, "tabWeighted_reg.csv")
write.csv(tabcart, "tabWeighted_cart.csv")
write.csv(tabwi, "tabWeighted_weightit.csv")

# exporting full dataset with weights
long_w <- long_a
write.csv(long_w, "CDI_long_analytic_weight.csv", row.names = F)


#### 2.2. marginal structural models ####

## Pulling in weighted dataset 
# long_w <- read.csv("CDI_long_analytic_weight.csv", header=T, na.strings = c("","NA"))

# function to calculate confidence intervals in GEE
confint.geeglm <- function(object, parm, level = 0.95, ...) {
  cc <- coef(summary(object))
  mult <- qnorm((1+level)/2)
  citab <- with(as.data.frame(cc),
                cbind(lwr=Estimate-mult*Std.err,
                      upr=Estimate+mult*Std.err))
  rownames(citab) <- rownames(cc)
  citab[parm,]
}

## Running generalized estimating equations non-weighted
gee0 <- geeglm(cdiff ~  currentabx_bin
               + cum_abx_days_yest+hosp_day_yest+hosp_any_f+priorabx_bin_f, #unbalanced 
               data=long_w, id=study_redcap_id, 
               family=binomial(link="logit"), corstr="ar1")
summary(gee0);QIC(gee0);
confint(gee0, level = 0.95)

### Running GEE with weights
## using Regression weights because they were the most balanced
## including variables that remained unbalanced after weighting
# Exposure = any antibiotic exposure
gee1 <- geeglm(cdiff ~  currentabx_bin
               + cum_abx_days_yest+hosp_day_yest+hosp_any_f+priorabx_bin_f, #unbalanced 
               data=long_w, id=study_redcap_id, weights=PS_weights_reg, 
               family=binomial(link="logit"), corstr="ar1")
summary(gee1);QIC(gee1);
confint(gee1, level = 0.95)

# Exposure = specific medication
gee2 <- geeglm(cdiff ~  b_lactamase_inhib_any + carbapenem_any + cephalosporin_any 
               + fluoroquinolone_any + monobactam_any + rifamycins_any + sulfonamides_any
               + misc_any 
               + cum_abx_days_yest+hosp_day_yest+hosp_any_f+priorabx_bin_f, #unbalanced 
               data=long_w, id=study_redcap_id, weights=PS_weights_reg, 
               family=binomial(link="logit"), corstr="ar1")
summary(gee2);QIC(gee2)
exp(confint(gee2, level = 0.95))
exp(gee2$coefficients)

# Exposure = proportion of hospitalization on antibiotics
gee3 <- geeglm(cdiff ~  cum_prop_abx_yest
               + hosp_any_f+priorabx_bin_f, #unbalanced
               # + (currentabx_bin*currentproton_f) + (priorabx_bin_f*prior_proton_f), #interactions  
               data=long_w, id=study_redcap_id, weights=PS_weights_reg, 
               family=binomial(link="logit"), corstr="ar1")
summary(gee3);QIC(gee3)
confint(gee3, level = 0.95)



