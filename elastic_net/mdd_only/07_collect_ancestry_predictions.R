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


  ## Load other ancestries
  SAS_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/SAS_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  SAS_qc_unique_id = paste(SAS_QC_all$FID, SAS_QC_all$IID, sep = "_")
  AFR_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/AFR_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  AFR_qc_unique_id = paste(AFR_QC_all$FID, AFR_QC_all$IID, sep = "_")


  ## Load data
  for( type in c("ALLCovar_wTWAS" , "ALLCovar_NoTWAS",  "NoCovar_wTWAS" , "NoCovar_NoTWAS") ){
    file_AFR = paste(dir, "/results_altancestry/LDref1kg_split", i, "_predictortype-",type,"_pred_covar0_AFR.csv", sep = "")
    pred_AFR = fread(file_AFR)
    colnames(pred_AFR) = paste(type, colnames(pred_AFR), sep = "_")
    file_SAS = paste(dir, "/results_altancestry/LDref1kg_split", i, "_predictortype-",type,"_pred_covar0_SAS.csv", sep = "")
    pred_SAS = fread(file_SAS)
    colnames(pred_SAS) = paste(type, colnames(pred_SAS), sep = "_")
    assign(paste(type,"AFR", sep = "_"),pred_AFR)
    assign(paste(type,"SAS", sep = "_"),pred_SAS)

    file_AFR_cov = paste(dir, "/results_altancestry/LDref1kg_split", i, "_predictortype-",type,"_pred_fullcovar0_AFR.csv", sep = "")
    pred_AFR_cov= fread(file_AFR_cov)
    colnames(pred_AFR_cov) = paste(type, colnames(pred_AFR_cov), sep = "_")
    file_SAS_cov = paste(dir, "/results_altancestry/LDref1kg_split", i, "_predictortype-",type,"_pred_fullcovar_SAS.csv", sep = "")
    pred_SAS_cov = fread(file_SAS_cov)
    colnames(pred_SAS_cov) = paste(type, colnames(pred_SAS_cov), sep = "_")
    assign(paste(type,"AFR_wcov", sep = "_"),pred_AFR_cov)
    assign(paste(type,"SAS_wcov", sep = "_"),pred_SAS_cov)
  }


  ## Combine
  cov_res_SAS = cbind(SAS_QC_all,
                  ALLCovar_wTWAS_SAS, ALLCovar_NoTWAS_SAS, NoCovar_wTWAS_SAS, NoCovar_NoTWAS_SAS)
  cov_res_AFR = cbind(AFR_QC_all,
                      ALLCovar_wTWAS_AFR, ALLCovar_NoTWAS_AFR, NoCovar_wTWAS_AFR, NoCovar_NoTWAS_AFR)

  cov_res_SAS_wcov = cbind(SAS_QC_all,
                      ALLCovar_wTWAS_SAS_wcov, ALLCovar_NoTWAS_SAS_wcov, NoCovar_wTWAS_SAS_wcov, NoCovar_NoTWAS_SAS_wcov)
  cov_res_AFR_wcov = cbind(AFR_QC_all,
                      ALLCovar_wTWAS_AFR_wcov, ALLCovar_NoTWAS_AFR_wcov, NoCovar_wTWAS_AFR_wcov, NoCovar_NoTWAS_AFR_wcov)


  ## Save
  fwrite(cov_res_SAS, paste(dir,"/results_altancestry/split",i,"_testset_covar0_pred_SAS.csv" ,sep = ""))
  fwrite(cov_res_AFR, paste(dir,"/results_altancestry/split",i,"_testset_covar0_pred_AFR.csv" ,sep = ""))


  fwrite(cov_res_SAS_wcov, paste(dir,"/results_altancestry/split",i,"_testset_fullcovar_pred_SAS.csv" ,sep = ""))
  fwrite(cov_res_AFR_wcov, paste(dir,"/results_altancestry/split",i,"_testset_fullcovar_pred_AFR.csv" ,sep = ""))

}
