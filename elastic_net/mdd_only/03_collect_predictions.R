# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(data.table)
library(tidyverse)


###############################################################################################

############################## COVAR 0 ######################################################

###############################################################################################


dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/single_trait/")

## Load covariates
filename_covar = mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/QC_MDD_covariate_info.txt")
covar = fread(filename_covar)
covar$covar_id = paste(covar$FID, covar$IID, sep = "_")

for(i in 1:10){
  print(i)


  ## Load Phenotype File to Determine Individuals
  testphenotype_file = paste(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/splits/split"),i,"_QC_MDD_testing.txt", sep = "")

  ## Load
  testpheno = fread(testphenotype_file)

  ## Obtain Unique ID
  test_id = paste(testpheno$FID, testpheno$IID, sep = "_")


  ## Obtain y train and y test
  y01.test = testpheno$MDD


  ## Obtain indices of covar to match train and test
  ind.test_cov =  match(test_id, covar$covar_id)

  ## Check that matched

  ## Covar for train and covar for test
  covar.test = as.data.frame(covar[ind.test_cov, -c("FID", "IID", "covar_id")])

  ## Change assessment center to factor
  covar.test$assessment_center = as.factor(covar.test$assessment_center)

  ## Load hapmap
  prscs_1kg_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,"_redo.txt", sep = "")
  prscs_1kg = fread(prscs_1kg_file)
  colnames(prscs_1kg) = "hapmap"

  ## Load logbf4
  prscs_1kg_file_logbf_over4 = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,".txt", sep = "")
  prscs_1kg_logbf_over4 = fread(prscs_1kg_file_logbf_over4)
  colnames(prscs_1kg_logbf_over4) = "logbf4"

  ## Load twassig
  prscs_1kg_file_twassig = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_twassig_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,".txt", sep = "")
  prscs_1kg_twassig = fread(prscs_1kg_file_twassig)
  colnames(prscs_1kg_twassig) = "TWAS"

  ## Load cr95set
  prscs_1kg_file_cr95set = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_cr95set_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,".txt", sep = "")
  prscs_1kg_cr95set = fread(prscs_1kg_file_cr95set)
  colnames(prscs_1kg_cr95set) = "cr95set"

  ## Load data
  for( type in c("ALLCovar_wTWAS" , "ALLCovar_NoTWAS",  "NoCovar_wTWAS" , "NoCovar_NoTWAS") ){
    file = paste(dir, "/results/LDref1kg_split", i, "_predictortype-",type,"_pred_covar0.csv", sep = "")
    pred = fread(file)
    colnames(pred) = paste(type, colnames(pred), sep = "_")
    assign(paste(type),pred)
  }


  ## Check same number of rows
  if(nrow(prscs_1kg) != length(y01.test) | nrow(prscs_1kg_cr95set)  != length(y01.test) |  nrow(prscs_1kg_logbf_over4)  != length(y01.test) |
     nrow(prscs_1kg_twassig)  != length(y01.test) |   nrow(ALLCovar_wTWAS) != length(y01.test) | nrow(ALLCovar_NoTWAS)  != length(y01.test) |
     nrow(NoCovar_wTWAS)  != length(y01.test) |  nrow(NoCovar_NoTWAS)  != length(y01.test) |  nrow(covar.test)  != length(y01.test) ){
    stop("Not all same rows")
  }

  ## Combine
  cov_res = cbind(y01.test,
                  prscs_1kg, prscs_1kg_cr95set, prscs_1kg_logbf_over4, prscs_1kg_twassig,
                  ALLCovar_wTWAS, ALLCovar_NoTWAS, NoCovar_wTWAS, NoCovar_NoTWAS, covar.test)


  ## Save
  fwrite(cov_res, paste(dir,"/results/split",i,"_testset_covar0_pred.csv" ,sep = ""))

}


###############################################################################################

############################## WITH COVAR ######################################################

###############################################################################################

dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/single_trait/")

## Load covariates
filename_covar = mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/QC_MDD_covariate_info.txt")
covar = fread(filename_covar)
covar$covar_id = paste(covar$FID, covar$IID, sep = "_")

for(i in 1:10){
  print(i)


  ## Load Phenotype File to Determine Individuals
  testphenotype_file = paste(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/splits/split"),i,"_QC_MDD_testing.txt", sep = "")

  ## Load
  testpheno = fread(testphenotype_file)

  ## Obtain Unique ID
  test_id = paste(testpheno$FID, testpheno$IID, sep = "_")


  ## Obtain y train and y test
  y01.test = testpheno$MDD


  ## Obtain indices of covar to match train and test
  ind.test_cov =  match(test_id, covar$covar_id)

  ## Check that matched

  ## Covar for train and covar for test
  covar.test = as.data.frame(covar[ind.test_cov, -c("FID", "IID", "covar_id")])

  ## Change assessment center to factor
  covar.test$assessment_center = as.factor(covar.test$assessment_center)

  ## Load hapmap
  prscs_1kg_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,"_redo.txt", sep = "")
  prscs_1kg = fread(prscs_1kg_file)
  colnames(prscs_1kg) = "hapmap"

  ## Load logbf4
  prscs_1kg_file_logbf_over4 = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,".txt", sep = "")
  prscs_1kg_logbf_over4 = fread(prscs_1kg_file_logbf_over4)
  colnames(prscs_1kg_logbf_over4) = "logbf4"

  ## Load twassig
  prscs_1kg_file_twassig = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_twassig_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,".txt", sep = "")
  prscs_1kg_twassig = fread(prscs_1kg_file_twassig)
  colnames(prscs_1kg_twassig) = "TWAS"

  ## Load cr95set
  prscs_1kg_file_cr95set = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_1kg_LD_ref_cr95set_pst_eff_a1_b0.5_phiauto_prediction_split_"),i,".txt", sep = "")
  prscs_1kg_cr95set = fread(prscs_1kg_file_cr95set)
  colnames(prscs_1kg_cr95set) = "cr95set"

  ## Load data
  for( type in c("ALLCovar_wTWAS" , "ALLCovar_NoTWAS",  "NoCovar_wTWAS" , "NoCovar_NoTWAS") ){
    file = paste(dir, "/results/LDref1kg_split", i, "_predictortype-",type,"_res.txt", sep = "")
    pred = fread(file)
    pred = pred[,c(2:5,7,9,11)]
    colnames(pred) = c("lasso", "enet_alphagrid", "enet_0.5", "ridge", "caret_kappa", "caret_deviance", "caret_roc")
    colnames(pred) = paste(type, colnames(pred), sep = "_")
    assign(paste(type),pred)
  }


  ## Check same number of rows
  if(nrow(prscs_1kg) != length(y01.test) | nrow(prscs_1kg_cr95set)  != length(y01.test) |  nrow(prscs_1kg_logbf_over4)  != length(y01.test) |
    nrow(prscs_1kg_twassig)  != length(y01.test) |   nrow(ALLCovar_wTWAS) != length(y01.test) | nrow(ALLCovar_NoTWAS)  != length(y01.test) |
    nrow(NoCovar_wTWAS)  != length(y01.test) |  nrow(NoCovar_NoTWAS)  != length(y01.test) |  nrow(covar.test)  != length(y01.test) ){
    stop("Not all same rows")
  }

  ## Combine
  cov_res = cbind(y01.test,
                  prscs_1kg, prscs_1kg_cr95set, prscs_1kg_logbf_over4, prscs_1kg_twassig,
                  ALLCovar_wTWAS, ALLCovar_NoTWAS, NoCovar_wTWAS, NoCovar_NoTWAS, covar.test)


  ## Save
  fwrite(cov_res, paste(dir,"/results/split",i,"_testset_wcovar_pred.csv" ,sep = ""))

}
