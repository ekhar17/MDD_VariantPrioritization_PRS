# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(purrr)
library(bigsnpr)
library(caret)
library(data.table)
library(tidyverse)
library(glmnet)
library(glmnetUtils)

## Arguments Loaded in
df_total = data.frame()
LD_ref =  "1kg"
i=commandArgs(trailingOnly=TRUE)[1] #1,2,3
type=commandArgs(trailingOnly=TRUE)[2]

set.seed(1)
print(i)
print(LD_ref)
print(type)

out_dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/single_trait/")


## Train Test file
testphenotype_file = paste(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/splits/split"),i,"_QC_MDD_testing.txt", sep = "")
trainphenotype_file = paste(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/splits/split"),i,"_QC_MDD_alltraining.txt", sep = "")
filename_covar = mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/QC_MDD_covariate_info.txt")


## Load
trainpheno = fread(trainphenotype_file)
testpheno = fread(testphenotype_file)
covar = fread(filename_covar)

## Obtain Unique ID
train_id = paste(trainpheno$FID, trainpheno$IID, sep = "_")
test_id = paste(testpheno$FID, testpheno$IID, sep = "_")
covar$covar_id = paste(covar$FID, covar$IID, sep = "_")


## Obtain y train and y test
y01.train = trainpheno$MDD
y01.test = testpheno$MDD

## Obtain indices of covar to match train and test
ind.train_cov =  match(train_id, covar$covar_id)
ind.test_cov =  match(test_id, covar$covar_id)

## Check that matched
if(!all.equal(covar$covar_id[ind.train_cov], train_id) | !all.equal(covar$covar_id[ind.test_cov], test_id) ){
  stop("something went wrong - train ids not match ind.train for covariance or test ids not match ind.test for covariance")
}

## Covar for train and covar for test
covar.train = as.data.frame(covar[ind.train_cov,-c("FID", "IID", "covar_id")])
covar.test = as.data.frame(covar[ind.test_cov, -c("FID", "IID", "covar_id")])

## Change assessment center to factor
covar.train$assessment_center = as.factor(covar.train$assessment_center)
covar.test$assessment_center = as.factor(covar.test$assessment_center)


## Match Id's
sample = fread(mdd_path("MDD_SAMPLE", ""))
head(sample,3)

## remove first individual
sample = sample[-1,]
print("remove placeholder")
head(sample,3)


## Object ID
object_id = paste(sample$ID_1, sample$ID_2, sep = "_")
ind.test =  match(test_id, object_id)
ind.train =  match(train_id, object_id)

## Get Ids of train and test
if(!all.equal(object_id[ind.test], test_id) | !all.equal(object_id[ind.train], train_id) ){
  print("something went wrong - test ids not match ind.test or train ids not match ind.train")
  stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
}


## Load MDD PRS CS files
prscs_file_twassig = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_twassig_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
prs_cs_twassig = read.table(prscs_file_twassig)

prscs_file_cr95set = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_cr95set_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
prs_cs_cr95set = read.table(prscs_file_cr95set)

prscs_file_logbf_over4 = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
prs_cs_logbf_over4 = read.table(prscs_file_logbf_over4)

if(!all.equal(nrow(sample), nrow(prs_cs_twassig)) | !all.equal(nrow(sample), nrow(prs_cs_cr95set))| !all.equal(nrow(sample), nrow(prs_cs_logbf_over4))){
  print("something went wrong - diff number of samples")
  stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
}


## Obtain train and test
prs_cs_twassig_test = prs_cs_twassig$x[ind.test]
prs_cs_twassig_train = prs_cs_twassig$x[ind.train]

prs_cs_cr95set_test = prs_cs_cr95set$x[ind.test]
prs_cs_cr95set_train = prs_cs_cr95set$x[ind.train]

prs_cs_logbf_over4_test = prs_cs_logbf_over4$x[ind.test]
prs_cs_logbf_over4_train = prs_cs_logbf_over4$x[ind.train]


## Create dataframe and scale
df_train_unscaled = data.frame(prs_cs_cr95set_train, prs_cs_logbf_over4_train, prs_cs_twassig_train, covar.train)
df_train_scale_full = df_train_unscaled
df_train_scale_full[,c(1:11,13)] = (scale(df_train_scale_full[,c(1:11,13)]))

df_test_unscaled = data.frame(prs_cs_cr95set_test, prs_cs_logbf_over4_test, prs_cs_twassig_test, covar.test)
df_test_scale_full = df_test_unscaled
df_test_scale_full[,c(1:11,13)] = (scale(df_test_scale_full[,c(1:11,13)]))


## Create character train
y01.train_character = y01.train
y01.train_character[as.logical(y01.train)] = "MDD"
y01.train_character[!as.logical(y01.train)] = "control"


## Tuning grid, and functions for caret and glmnet
tg <- expand.grid(lambda=c(.01,.05,.1,.25,.5,.75,0.9,0.95,1), alpha=c(0, 0.05, 0.1,.25,.5,.75,0.9,0.95,1))
get_alpha <- function(fit) {
  alpha <- fit$alpha
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  alpha[which.min(error)]
}
get_model_params <- function(fit) {
  alpha <- fit$alpha
  lambdaMin <- sapply(fit$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(fit$modlist, `[[`, "lambda.1se")
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  best <- which.min(error)
  data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
             lambdaSE = lambdaSE[best], eror = error[best])
}
fitControl <- trainControl(## 10-fold CV)
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  ## repeated ten times
  classProbs = T)
dev_summary <- function(data, lev = NULL, model = NULL) {
  is_class1 <- ifelse(data$obs == lev[1], 1, 0)
  prob_class1 <- data[, lev[1]]

  c(deviance = -2*sum(is_class1*log(prob_class1) +
                        ((1-is_class1)*log(1-prob_class1))),
    twoClassSummary(data, lev = lev))
}
ctrl <- trainControl(summaryFunction = dev_summary,
                     classProbs = TRUE,
                     method = "repeatedcv",
                     number = 10,
                     repeats = 10)


## Check that class factor
if(class(df_train_scale_full$assessment_center) != "factor"){
  stop("ASSESSMENT CENTER NEEDS TO BE FACTOR")
}

## Loop over types
#for( type in c("ALLCovar_wTWAS" , "ALLCovar_NoTWAS",  "NoCovar_wTWAS" , "NoCovar_NoTWAS") ){

  ## Columns keep for each type
      if(type == "ALLCovar_wTWAS"){
      print(type)
      col_index = 1:15
      print(head(df_train_scale_full[,col_index]))
      print(head(df_test_scale_full[,col_index]))
    }
    if(type == "ALLCovar_NoTWAS"){
  print(type)
  col_index = c(1:2,4:15)
  print(head(df_train_scale_full[,col_index]))
  print(head(df_test_scale_full[,col_index]))
}
if(type == "NoCovar_wTWAS"){
  print(type)
  col_index = c(1:3)
  print(head(df_train_scale_full[,col_index]))
  print(head(df_test_scale_full[,col_index]))
}
if(type == "NoCovar_NoTWAS"){
  print(type)
  col_index = c(1:2)
  print(head(df_train_scale_full[,col_index]))
  print(head(df_test_scale_full[,col_index]))

  }


## Train and test column keep
    df_train_scale = df_train_scale_full[,col_index]
  df_test_scale = df_test_scale_full[,col_index]

  ## Glmnet models
  lasso_scale = cv.glmnet(x = as.matrix(df_train_scale),
                    y = y01.train, family = "binomial")
  print("cv glmnet lasso scale")
  print(lasso_scale)

  glmnet_scale_alphagrid = cva.glmnet(x = as.matrix(df_train_scale),
                y = y01.train, family = "binomial")
  print("cv glmnet alpha & lasso scale")
  print(glmnet_scale_alphagrid)

  glmnet_alpha0.5_scale = cv.glmnet(x = as.matrix(df_train_scale),
                       y = y01.train, family = "binomial",alpha = 0.5)
  print("results alpha - 0.5 enet scale")
  print(glmnet_alpha0.5_scale)

  ridge_scale = cv.glmnet(x = as.matrix(df_train_scale),
                     y = y01.train, family = "binomial",alpha = 0)
  print("results ridge scale")
  print(ridge_scale)


  ## Obtain coefficients
  coef_df = cbind(coef(lasso_scale),
      coef(ridge_scale),
      coef(glmnet_scale_alphagrid, alpha = get_alpha(glmnet_scale_alphagrid)),
      coef(glmnet_alpha0.5_scale))
  colnames(coef_df) = c("lasso", "ridge", "alphagrid", "alpha0.5")

  ## Obtain test set predictions
  p_lasso_scale = predict(lasso_scale, newx = as.matrix(df_test_scale))
  p_enet_scale = predict(glmnet_scale_alphagrid, newx = as.matrix(df_test_scale), alpha = get_alpha(glmnet_scale_alphagrid))
  p_enet0.5_scale = predict(glmnet_alpha0.5_scale, newx = as.matrix(df_test_scale))
  p_ridge_scale = predict(ridge_scale, newx = as.matrix(df_test_scale))


  ## Caret models
  caret_model_Kappa <- train(
    x = df_train_scale,
    y = y01.train_character,
    method = "glmnet",
    trControl = fitControl,
    tuneGrid=tg,
    verbose = T,metric = "Kappa"
  )
  best_lambda_Kappa <- caret_model_Kappa$bestTune$lambda
  print("results - train glmnet kappa center scale")
  print(caret_model_Kappa)

  deviance_model <- train(
    x = df_train_scale,
    y = y01.train_character,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid=tg,
    verbose = T
  )
  best_lambda_deviance <- deviance_model$bestTune$lambda
  print("results - train glmnet deviance center scale")
  print(deviance_model)

  roc_model <- train(
    x = df_train_scale,
    y = y01.train_character,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid=tg,
    verbose = T, metric = "ROC"
  )
  best_lambda_roc <- roc_model$bestTune$lambda
  print("results - train glmnet ROC center scale")
  print(roc_model)


  ## Coefficients
  coef_df_caret = cbind(coef(caret_model_Kappa$finalModel, s = best_lambda_Kappa)[, 1],
                        coef(deviance_model$finalModel, s = best_lambda_deviance)[, 1],
                        coef(roc_model$finalModel, s = best_lambda_roc)[, 1])
  colnames(coef_df_caret) = c("caret_kappa", "caret_deviance", "caret_roc")

  ## Test set predictions
  colnames(df_test_scale) = colnames(df_train_scale)
  p_caret_kappa = predict(caret_model_Kappa, newdata = df_test_scale, type = "prob")
  p_caret_deviance = predict(deviance_model, newdata = df_test_scale, type = "prob")
  p_caret_roc = predict(roc_model, newdata = df_test_scale, type = "prob")


## combine results
res = data.frame(p_lasso_scale, p_enet_scale ,
                p_enet0.5_scale, p_ridge_scale,
                p_caret_kappa, p_caret_deviance, p_caret_roc)

## Save
out_file = paste(out_dir, "/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_res.txt", sep = "")
print(out_file)
write.table(res, out_file)

coef_file = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_coefglmnet.txt", sep = "")
print(coef_file)
write.table(as.matrix(coef_df), coef_file)

coef_file_caret = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_coefcaret.txt", sep = "")
print(coef_file_caret)
write.table(as.matrix(coef_df_caret), coef_file_caret)


## Save models
all_models = list(lasso_scale, ridge_scale, glmnet_scale_alphagrid, glmnet_alpha0.5_scale,
     caret_model_Kappa, deviance_model,roc_model)
file_models = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_models.rds", sep = "")
print(file_models )
saveRDS(all_models, file_models )


#}
