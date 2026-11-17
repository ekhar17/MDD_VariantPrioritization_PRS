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

LD_ref = "1kg"
i=commandArgs(trailingOnly=TRUE)[1] #1,2,3
type=commandArgs(trailingOnly=TRUE)[2]


set.seed(1)
print(i)
print(LD_ref)
print(type)


out_dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/multi_trait/")


##  Test file
testphenotype_file = paste(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/splits/split"),i,"_QC_MDD_testing.txt", sep = "")
filename_covar = mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/QC_MDD_covariate_info.txt")


## Load
testpheno = fread(testphenotype_file)
covar = fread(filename_covar)

## Obtain Unique ID
test_id = paste(testpheno$FID, testpheno$IID, sep = "_")
covar$covar_id = paste(covar$FID, covar$IID, sep = "_")


## Obtain indices of covar to match train and test
ind.test_cov =  match(test_id, covar$covar_id)

## Check that matched
if( !all.equal(covar$covar_id[ind.test_cov], test_id) ){
  stop("something went wrong - train ids not match ind.train for covariance or test ids not match ind.test for covariance")
}

## Covar for train and covar for test
covar.test = as.data.frame(covar[ind.test_cov, -c("FID", "IID", "covar_id")])

## Change assessment center to factor
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

## Get Ids of train and test
if(!all.equal(object_id[ind.test], test_id) ){
  print("something went wrong - test ids not match ind.test or train ids not match ind.train")
  stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
}


## Load MDD PRS CS files
prscs_file_cr95set = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_cr95set_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
prs_cs_cr95set = read.table(prscs_file_cr95set)

prscs_file_logbf_over4 = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
prs_cs_logbf_over4 = read.table(prscs_file_logbf_over4)

if( !all.equal(nrow(sample), nrow(prs_cs_cr95set))| !all.equal(nrow(sample), nrow(prs_cs_logbf_over4))){
  print("something went wrong - diff number of samples")
  stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
}


########### Load BIP, SCZ, ADHD prediction #########################


prscs_file_cr95set_BIP = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/BIP_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_cr95set_BIP = read.table(prscs_file_cr95set_BIP)
prscs_file_logbf_over4_BIP = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/BIP_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_logbf_over4_BIP = read.table(prscs_file_logbf_over4_BIP)


prscs_file_cr95set_SCZ = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/SCZ_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_cr95set_SCZ = read.table(prscs_file_cr95set_SCZ)
prscs_file_logbf_over4_SCZ = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/SCZ_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_logbf_over4_SCZ = read.table(prscs_file_logbf_over4_SCZ )


prscs_file_cr95set_ADHD = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ADHD_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_cr95set_ADHD = read.table(prscs_file_cr95set_ADHD)
prscs_file_logbf_over4_ADHD = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ADHD_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_logbf_over4_ADHD = read.table(prscs_file_logbf_over4_ADHD)


prscs_file_cr95set_ANX = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ANXmetasubtract_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_cr95set_ANX = read.table(prscs_file_cr95set_ANX)
prscs_file_logbf_over4_ANX = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ANXmetasubtract_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
prs_cs_logbf_over4_ANX = read.table(prscs_file_logbf_over4_ANX)

if(!all.equal(nrow(sample), nrow(prs_cs_cr95set_BIP))| !all.equal(nrow(sample), nrow(prs_cs_logbf_over4_BIP))|!all.equal(nrow(sample), nrow(prs_cs_cr95set_SCZ))| !all.equal(nrow(sample), nrow(prs_cs_logbf_over4_SCZ))| !all.equal(nrow(sample), nrow(prs_cs_cr95set_ADHD))| !all.equal(nrow(sample), nrow(prs_cs_logbf_over4_ADHD)) |
   !all.equal(nrow(sample), nrow(prs_cs_cr95set_ANX))| !all.equal(nrow(sample), nrow(prs_cs_logbf_over4_ANX))
){
  print("something went wrong - diff number of samples")
  stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
}


########################### Obtain train and test #########################
prs_cs_cr95set_test = prs_cs_cr95set$x[ind.test]
prs_cs_logbf_over4_test = prs_cs_logbf_over4$x[ind.test]


########################### Repeat for other traits #########################
prs_cs_cr95set_BIP_test = prs_cs_cr95set_BIP$x[ind.test]
prs_cs_logbf_over4_BIP_test = prs_cs_logbf_over4_BIP$x[ind.test]

prs_cs_cr95set_SCZ_test = prs_cs_cr95set_SCZ$x[ind.test]
prs_cs_logbf_over4_SCZ_test = prs_cs_logbf_over4_SCZ$x[ind.test]

prs_cs_cr95set_ADHD_test = prs_cs_cr95set_ADHD$x[ind.test]
prs_cs_logbf_over4_ADHD_test = prs_cs_logbf_over4_ADHD$x[ind.test]


prs_cs_cr95set_ANX_test = prs_cs_cr95set_ANX$x[ind.test]
prs_cs_logbf_over4_ANX_test = prs_cs_logbf_over4_ANX$x[ind.test]


## Create dataframe and scale
df_test_unscaled = data.frame(prs_cs_cr95set_test, prs_cs_logbf_over4_test,
                              prs_cs_cr95set_BIP_test, prs_cs_logbf_over4_BIP_test,
                              prs_cs_cr95set_SCZ_test, prs_cs_logbf_over4_SCZ_test,
                              prs_cs_cr95set_ADHD_test, prs_cs_logbf_over4_ADHD_test,
                              prs_cs_cr95set_ANX_test, prs_cs_logbf_over4_ANX_test,
                              covar.test)
df_test_scale_full = df_test_unscaled
df_test_scale_full[,c(1:18,20)] = (scale(df_test_scale_full[,c(1:18,20)]))


## Make all covar0
df_test_scale_full_nocovar = df_test_scale_full
df_test_scale_full_nocovar[,11:22] = 0
df_test_scale_full_nocovar_caret = df_test_scale_full_nocovar
colnames(df_test_scale_full_nocovar_caret)[1:10] = c("prs_cs_cr95set_train",
                                                    "prs_cs_logbf_over4_train",
                                                    "prs_cs_cr95set_BIP_train" ,
                                                    "prs_cs_logbf_over4_BIP_train",
                                                    "prs_cs_cr95set_SCZ_train",
                                                    "prs_cs_logbf_over4_SCZ_train",
                                                    "prs_cs_cr95set_ADHD_train",
                                                    "prs_cs_logbf_over4_ADHD_train",
                                                    "prs_cs_cr95set_ANX_train" ,
                                                    "prs_cs_logbf_over4_ANX_train" )

## Columns keep for each type
if(type == "cr95setlogbf4_ALLCovar"){
  print(type)
  col_index = 1:22
  print(head(df_test_scale_full[,col_index]))
  col_apply = 1:10
  index_coef = 2:11
}
if(type == "cr95setonly_ALLCovar"){
  print(type)
  col_index = c(1,3,5,7,9,11:22)
  print(head(df_test_scale_full[,col_index]))
  col_apply = 1:5
  index_coef = 2:6


}
if(type ==  "cr95setlogbf4_NoCovar"){
  print(type)
  col_index = c(1:10)
  print(head(df_test_scale_full[,col_index]))
  col_apply = 1:10
  index_coef = 2:11


}
if(type == "cr95setonly_NoCovar"){
  print(type)
  col_index = c(1,3,5,7,9)
  print(head(df_test_scale_full[,col_index]))
  col_apply = 1:5
  index_coef = 2:6

}

df_test_scale_full_nocovar = df_test_scale_full_nocovar[,col_index]
df_test_scale_full_nocovar_caret = df_test_scale_full_nocovar_caret[,col_index]
head(df_test_scale_full_nocovar)
head(df_test_scale_full_nocovar_caret)
head(df_test_scale_full[,col_apply])

## Load models
file_models = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_models_multitrait.rds", sep = "")
models = readRDS(file_models)

## Load coefficients Glmnet
coef_file = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_coefglmnet_multitrait.txt", sep = "")
coef_df = fread(coef_file) %>% as.data.frame()
rownames(coef_df) = coef_df$V1
coef_df = coef_df[,-1]

## Load coefficients caret
coef_file_caret = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_coefcaret_multitrait.txt", sep = "")
coef_df_caret = fread(coef_file_caret) %>% as.data.frame()
rownames(coef_df_caret) = coef_df_caret$V1
coef_df_caret = coef_df_caret[,-1]

## Combine
coef_full = cbind(coef_df, coef_df_caret)
coef_full[index_coef,]

## Helper function
get_alpha <- function(fit) {
  alpha <- fit$alpha
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  alpha[which.min(error)]
}

pred_test_covar0 = matrix(NA, ncol = 7, nrow = nrow(df_test_scale_full_nocovar))

## Apply to nocovar test set
for(k in 1:7){

  print(k)
  ## Apply via predict
  if(k %in% c(1,2,4)){
    predict_funct_mat = predict(models[[k]], newx = as.matrix(df_test_scale_full_nocovar))
    predict_funct = predict_funct_mat[,1]
  }
  if(k == 3){
    predict_funct_mat = predict(models[[k]], newx = as.matrix(df_test_scale_full_nocovar), get_alpha(models[[k]]))
    predict_funct = predict_funct_mat[,1]

  }
  if(k %in% c(5:7)){
    predict_funct_both = predict(models[[k]], newdata = df_test_scale_full_nocovar_caret, type = "prob")
    predict_funct = predict_funct_both[,2]
  }

  ## Multiply saved coefficients
  if( k %in% c(1:4)){
    multiply_coef = as.matrix(df_test_scale_full_nocovar[,col_apply]) %*% matrix(coef_full[index_coef,k] ) + coef_full[1,k]
  }
  if( k %in% c(5:7)){
    logit = as.matrix(df_test_scale_full_nocovar[,col_apply]) %*% matrix(coef_full[index_coef,k] ) + coef_full[1,k]
    multiply_coef = exp(logit )/(1+exp(logit ))
  }

  if(!all.equal(predict_funct, c(multiply_coef))){
    stop("Not all equal - check")
  }

  pred_test_covar0[,k] = predict_funct
}


colnames(pred_test_covar0) = c("lasso", "ridge", "enet_alphagrid", "enet0.5",
                               "caret_kappa", "caret_deviance", "caret_roc")
write.csv(pred_test_covar0,
          paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_pred_covar0_multitrait.csv", sep = ""),
          row.names = F
)
print( paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_pred_covar0multitrait.csv", sep = ""))
