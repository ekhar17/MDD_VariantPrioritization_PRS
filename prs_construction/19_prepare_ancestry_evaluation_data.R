# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

##
library(data.table)
library(tidyverse)

## Load sample ID
sample = fread(mdd_path("MDD_SAMPLE", ""))
head(sample,3)
## remove first individual
sample = sample[-1,]
print("remove placeholder")
head(sample,3)

## Load PRS results
for(LD_ref in  c("1kg")){


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
                      hapmap_prs = -prs_cs_hapmap$x,
                      cr95set_prs = -prs_cs_cr95set$x,
                      logbf_over4_prs = -prs_cs_logbf_over4$x,
                      twassig_prs = -prs_cs_twassig$x)


  ## Create object ID
  object_id = paste(sample$ID_1, sample$ID_2, sep = "_")


  ## Load EUR
  EUR_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/EUR_Single_Het_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  EUR_qc_unique_id = paste(EUR_QC_all$FID, EUR_QC_all$IID, sep = "_")

  ## Add covar info
  EUR_covar = fread( mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/QC_MDD_covariate_info.txt"))
  EUR_covar_unique_id = paste(EUR_covar$FID, EUR_covar$IID, sep = "_")
  if(!all.equal(EUR_qc_unique_id, EUR_covar_unique_id)){
    stop("ID not match")
  }
  EUR_QC_all = cbind(EUR_QC_all, EUR_covar[,3:14])

  ## Load other ancestries
  SAS_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/SAS_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  SAS_qc_unique_id = paste(SAS_QC_all$FID, SAS_QC_all$IID, sep = "_")
  AFR_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/AFR_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  AFR_qc_unique_id = paste(AFR_QC_all$FID, AFR_QC_all$IID, sep = "_")

  ## Add predictions to all
  if(!all.equal(object_id[match(EUR_qc_unique_id,object_id)], EUR_qc_unique_id) &
  all.equal(object_id[match(SAS_qc_unique_id ,object_id)], SAS_qc_unique_id ) &
  all.equal(object_id[match(AFR_qc_unique_id ,object_id)], AFR_qc_unique_id ) ){
    stop("Not all ID equal - check")
  }

  ## Add info
  ## Combine prediction info, MDD status, and covariate info
  df_EUR = cbind(df_all[match(EUR_qc_unique_id,object_id),], EUR_QC_all[,3:15])
  df_AFR = cbind(df_all[match(AFR_qc_unique_id,object_id),], AFR_QC_all[,3:15])
  df_SAS = cbind(df_all[match(SAS_qc_unique_id,object_id),], SAS_QC_all[,3:15])

  ## Add covar info
  all.equal(df_EUR$FID, EUR_QC_all$FID)
  all.equal(df_AFR$FID, AFR_QC_all$FID)
  all.equal(df_SAS$FID, SAS_QC_all$FID)


  ## Check AUC
  print("EUR - Hapmap, cr95set, logbf4, twassig")
  print(bigstatsr::AUC(df_EUR$hapmap_prs, df_EUR$MDD))
  print(bigstatsr::AUC(df_EUR$cr95set_prs, df_EUR$MDD))
  print(bigstatsr::AUC(df_EUR$logbf_over4_prs, df_EUR$MDD))
  print(bigstatsr::AUC(df_EUR$twassig_prs, df_EUR$MDD))

  ## Check AUC
  print("AFR - Hapmap, cr95set, logbf4, twassig")
  print(bigstatsr::AUC(df_AFR$hapmap_prs, df_AFR$MDD))
  print(bigstatsr::AUC(df_AFR$cr95set_prs, df_AFR$MDD))
  print(bigstatsr::AUC(df_AFR$logbf_over4_prs, df_AFR$MDD))
  print(bigstatsr::AUC(df_AFR$twassig_prs, df_AFR$MDD))


  ## Check AUC
  print("SAS - Hapmap, cr95set, logbf4, twassig")
  print(bigstatsr::AUC(df_SAS$hapmap_prs, df_SAS$MDD))
  print(bigstatsr::AUC(df_SAS$cr95set_prs, df_SAS$MDD))
  print(bigstatsr::AUC(df_SAS$logbf_over4_prs, df_SAS$MDD))
  print(bigstatsr::AUC(df_SAS$twassig_prs, df_SAS$MDD))

  ## Save data
  fwrite(df_EUR, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_EUR_clean_prediction.txt"))
  fwrite(df_SAS, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_SAS_clean_prediction.txt"))
  fwrite(df_AFR, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_AFR_clean_prediction.txt"))


}
