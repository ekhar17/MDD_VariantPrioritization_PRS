# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(data.table)
library(tidyverse)
library(glmnet)
library(glmnetUtils)
library(caret)

for(i in 1:10){
  print(i)
  for(type in c( "ALLCovar_wTWAS","ALLCovar_NoTWAS", "NoCovar_wTWAS","NoCovar_NoTWAS")){
## Load sample ID
sample = fread(mdd_path("MDD_SAMPLE", ""))
head(sample,3)
## remove first individual
sample = sample[-1,]
print("remove placeholder")
head(sample,3)

## Load PRS results
LD_ref  = "1kg"


  ## Load predictions
  prscs_file_hapmap = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_pst_eff_a1_b0.5_phiauto_update_prediction.txt", sep = "")
  prs_cs_hapmap = read.table(prscs_file_hapmap)

  prscs_file_twassig = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_twassig_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
  prs_cs_twassig = read.table(prscs_file_twassig)

  prscs_file_logbf_over4 = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
  prs_cs_logbf_over4 = read.table(prscs_file_logbf_over4)

  prscs_file_cr95set = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_cr95set_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
  prs_cs_cr95set = read.table(prscs_file_cr95set)

  ## Create dataframe
  df_all = data.frame(FID = sample$ID_1, IID = sample$ID_2,
                      cr95set_prs = -prs_cs_cr95set$x,
                      logbf_over4_prs = -prs_cs_logbf_over4$x,
                      twassig_prs = -prs_cs_twassig$x)


  ## Create object ID
  object_id = paste(sample$ID_1, sample$ID_2, sep = "_")


  ## Load other ancestries
  SAS_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/SAS_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  SAS_qc_unique_id = paste(SAS_QC_all$FID, SAS_QC_all$IID, sep = "_")
  AFR_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/AFR_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  AFR_qc_unique_id = paste(AFR_QC_all$FID, AFR_QC_all$IID, sep = "_")

  ## Add predictions to all
  if(!
     all.equal(object_id[match(SAS_qc_unique_id ,object_id)], SAS_qc_unique_id ) &
     all.equal(object_id[match(AFR_qc_unique_id ,object_id)], AFR_qc_unique_id ) ){
    stop("Not all ID equal - check")
  }

  ## Add info
  ## Combine prediction info, MDD status, and covariate info
  df_AFR = cbind(df_all[match(AFR_qc_unique_id,object_id),], AFR_QC_all[,3:15])
  df_SAS = cbind(df_all[match(SAS_qc_unique_id,object_id),], SAS_QC_all[,3:15])


  ## Create dataframe and scale
  df_test_unscaled_AFR = df_AFR[,c(3:5,7:18)]
  df_test_unscaled_SAS = df_SAS[,c(3:5,7:18)]
  df_test_unscaled_AFR$assessment_center = as.factor(df_test_unscaled_AFR$assessment_center)
  df_test_unscaled_SAS$assessment_center = as.factor(df_test_unscaled_SAS$assessment_center)

  df_test_scale_full_AFR = df_test_unscaled_AFR
  df_test_scale_full_AFR[,c(1:11,13)] = (scale(df_test_scale_full_AFR[,c(1:11,13)]))


  df_test_scale_full_SAS = df_test_unscaled_SAS
  df_test_scale_full_SAS[,c(1:11,13)] = (scale(df_test_scale_full_SAS[,c(1:11,13)]))


  ## Make all covar0
  df_test_scale_full_nocovar_AFR = df_test_scale_full_AFR
  df_test_scale_full_nocovar_AFR[,4:15] = 0
  df_test_scale_full_nocovar_caret_AFR = df_test_scale_full_nocovar_AFR
  colnames(df_test_scale_full_nocovar_caret_AFR)[1:3] = c("prs_cs_cr95set_train",
                                                      "prs_cs_logbf_over4_train",
                                                      "prs_cs_twassig_train")

  df_test_scale_full_nocovar_SAS = df_test_scale_full_SAS
  df_test_scale_full_nocovar_SAS[,4:15] = 0
  df_test_scale_full_nocovar_caret_SAS = df_test_scale_full_nocovar_SAS
  colnames(df_test_scale_full_nocovar_caret_SAS)[1:3] = c("prs_cs_cr95set_train",
                                                          "prs_cs_logbf_over4_train",
                                                          "prs_cs_twassig_train")

  df_test_scale_full_wcovar_caret_SAS = df_test_scale_full_SAS
  colnames(df_test_scale_full_wcovar_caret_SAS)[1:3] = c("prs_cs_cr95set_train",
                                                          "prs_cs_logbf_over4_train",
                                                          "prs_cs_twassig_train")

  df_test_scale_full_wcovar_caret_AFR = df_test_scale_full_AFR
  colnames(df_test_scale_full_wcovar_caret_AFR)[1:3] = c("prs_cs_cr95set_train",
                                                         "prs_cs_logbf_over4_train",
                                                         "prs_cs_twassig_train")


  out_dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/single_trait/")


  ## Columns keep for each type
  if(type == "ALLCovar_wTWAS"){
    print(type)
    col_index = 1:15
   # print(head(df_test_scale_full[,col_index]))
    col_apply = 1:3
    index_coef = 2:4


  }
  if(type == "ALLCovar_NoTWAS"){
    print(type)
    col_index = c(1:2,4:15)
   # print(head(df_test_scale_full[,col_index]))
    col_apply = 1:2
    index_coef = 2:3


  }
  if(type == "NoCovar_wTWAS"){
    print(type)
    col_index = c(1:3)
    #print(head(df_test_scale_full[,col_index]))
    col_apply = 1:3
    index_coef = 2:4


  }
  if(type == "NoCovar_NoTWAS"){
    print(type)
    col_index = c(1:2)
   # print(head(df_test_scale_full[,col_index]))
    col_apply = 1:2
    index_coef = 2:3


  }

  df_test_scale_full_nocovar_AFR = df_test_scale_full_nocovar_AFR[,col_index]
  df_test_scale_full_nocovar_caret_AFR = df_test_scale_full_nocovar_caret_AFR[,col_index]
  head(df_test_scale_full_nocovar_AFR)
  head(df_test_scale_full_nocovar_caret_AFR)
  head(df_test_scale_full_AFR[,col_apply])

  df_test_scale_full_nocovar_SAS = df_test_scale_full_nocovar_SAS[,col_index]
  df_test_scale_full_nocovar_caret_SAS = df_test_scale_full_nocovar_caret_SAS[,col_index]
  head(df_test_scale_full_nocovar_SAS)
  head(df_test_scale_full_nocovar_caret_SAS)
  head(df_test_scale_full_SAS[,col_apply])

  df_test_scale_full_wcovar_AFR = df_test_scale_full_AFR[,col_index]
  df_test_scale_full_wcovar_caret_AFR = df_test_scale_full_wcovar_caret_AFR[,col_index]
  head(df_test_scale_full_wcovar_AFR)
  head(df_test_scale_full_wcovar_caret_AFR)
  df_test_scale_full_wcovar_SAS = df_test_scale_full_SAS[,col_index]
  df_test_scale_full_wcovar_caret_SAS = df_test_scale_full_wcovar_caret_SAS[,col_index]
  head(df_test_scale_full_wcovar_SAS)
  head(df_test_scale_full_wcovar_caret_SAS)

  ## Load models
  file_models = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_models.rds", sep = "")
  models = readRDS(file_models)

  ## Load coefficients Glmnet
  coef_file = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_coefglmnet.txt", sep = "")
  coef_df = fread(coef_file) %>% as.data.frame()
  rownames(coef_df) = coef_df$V1
  coef_df = coef_df[,-1]

  ## Load coefficients caret
  coef_file_caret = paste(out_dir,"/results/LDref", LD_ref, "_split", i, "_predictortype-",type,"_coefcaret.txt", sep = "")
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

  pred_test_covar0_AFR = matrix(NA, ncol = 7, nrow = nrow(df_test_scale_full_nocovar_AFR))
  pred_test_covar0_SAS = matrix(NA, ncol = 7, nrow = nrow(df_test_scale_full_nocovar_SAS))

  pred_test_fullcovar_AFR = matrix(NA, ncol = 7, nrow = nrow(df_test_scale_full_nocovar_AFR))
  pred_test_fullcovar_SAS = matrix(NA, ncol = 7, nrow = nrow(df_test_scale_full_nocovar_SAS))

  ## Apply to nocovar test set
  for(k in 1:7){

    print(k)
    ## Apply via predict
    if(k %in% c(1,2,4)){
      predict_funct_mat_AFR = predict(models[[k]], newx = as.matrix(df_test_scale_full_nocovar_AFR))
      predict_funct_AFR = predict_funct_mat_AFR[,1]

      predict_funct_mat_SAS = predict(models[[k]], newx = as.matrix(df_test_scale_full_nocovar_SAS))
      predict_funct_SAS = predict_funct_mat_SAS[,1]

      predict_funct_mat_AFR_full = predict(models[[k]], newx = as.matrix(df_test_scale_full_wcovar_AFR))
      predict_funct_AFR_full = predict_funct_mat_AFR_full[,1]

      predict_funct_mat_SAS_full = predict(models[[k]], newx = as.matrix(df_test_scale_full_wcovar_SAS))
      predict_funct_SAS_full = predict_funct_mat_SAS_full[,1]
    }
    if(k == 3){
      predict_funct_mat_AFR = predict(models[[k]], newx = as.matrix(df_test_scale_full_nocovar_AFR), get_alpha(models[[k]]))
      predict_funct_AFR = predict_funct_mat_AFR[,1]

      predict_funct_mat_SAS = predict(models[[k]], newx = as.matrix(df_test_scale_full_nocovar_SAS), get_alpha(models[[k]]))
      predict_funct_SAS = predict_funct_mat_SAS[,1]

      predict_funct_mat_AFR_full = predict(models[[k]], newx = as.matrix(df_test_scale_full_wcovar_AFR), get_alpha(models[[k]]))
      predict_funct_AFR_full = predict_funct_mat_AFR_full[,1]

      predict_funct_mat_SAS_full = predict(models[[k]], newx = as.matrix(df_test_scale_full_wcovar_SAS), get_alpha(models[[k]]))
      predict_funct_SAS_full = predict_funct_mat_SAS_full[,1]

    }
    if(k %in% c(5:7)){
      predict_funct_both_AFR = predict(models[[k]], newdata = df_test_scale_full_nocovar_caret_AFR, type = "prob")
      predict_funct_AFR = predict_funct_both_AFR[,2]

      predict_funct_both_SAS = predict(models[[k]], newdata = df_test_scale_full_nocovar_caret_SAS, type = "prob")
      predict_funct_SAS = predict_funct_both_SAS[,2]


      predict_funct_both_AFR_full = predict(models[[k]], newdata = df_test_scale_full_wcovar_caret_AFR, type = "prob")
      predict_funct_AFR_full = predict_funct_both_AFR_full[,2]

      predict_funct_both_SAS_full = predict(models[[k]], newdata = df_test_scale_full_wcovar_caret_SAS, type = "prob")
      predict_funct_SAS_full = predict_funct_both_SAS_full[,2]


    }


    pred_test_covar0_AFR[,k] = predict_funct_AFR
    pred_test_covar0_SAS[,k] = predict_funct_SAS


    pred_test_fullcovar_AFR[,k] = predict_funct_AFR_full
    pred_test_fullcovar_SAS[,k] = predict_funct_SAS_full
  }


  colnames(pred_test_covar0_AFR) =colnames(pred_test_covar0_SAS) =  colnames(pred_test_fullcovar_AFR) =colnames(pred_test_fullcovar_SAS) =
    c("lasso", "ridge", "enet_alphagrid", "enet0.5",
                                 "caret_kappa", "caret_deviance", "caret_roc")
  write.csv(pred_test_covar0_AFR,
            paste(out_dir,"/results_altancestry/LDref", LD_ref, "_split", i, "_predictortype-",type,"_pred_covar0_AFR.csv", sep = ""),
            row.names = F
  )
  write.csv(pred_test_covar0_SAS,
            paste(out_dir,"/results_altancestry/LDref", LD_ref, "_split", i, "_predictortype-",type,"_pred_covar0_SAS.csv", sep = ""),
            row.names = F
  )

  write.csv(pred_test_fullcovar_AFR,
            paste(out_dir,"/results_altancestry/LDref", LD_ref, "_split", i, "_predictortype-",type,"_pred_fullcovar0_AFR.csv", sep = ""),
            row.names = F
  )
  write.csv(pred_test_fullcovar_SAS,
            paste(out_dir,"/results_altancestry/LDref", LD_ref, "_split", i, "_predictortype-",type,"_pred_fullcovar_SAS.csv", sep = ""),
            row.names = F
  )

}

}
