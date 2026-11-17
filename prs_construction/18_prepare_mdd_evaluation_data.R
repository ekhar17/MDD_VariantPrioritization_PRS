# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))


## Load sample ID
sample = fread(mdd_path("MDD_SAMPLE", ""))
head(sample,3)
## remove first individual
sample = sample[-1,]
print("remove placeholder")
head(sample,3)

## Load PRS results
for(LD_ref in  c("1kg")){

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


  ## Load clean data
  QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/EUR_Single_Het_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
  qc_unique_id = paste(QC_all$FID, QC_all$IID, sep = "_")
  covar = fread( mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/QC_MDD_covariate_info.txt"))


  ## Keep info regarding only clean data
  df_keep =df_all[match(qc_unique_id,object_id),]
  all.equal(df_keep$FID, QC_all$FID)
  all.equal(df_keep$IID, QC_all$IID)

  ## Add MDD info
  df_keep$MDD = QC_all$MDD


  ## Add covar info
  all.equal(df_keep$FID, covar$FID)
  all.equal(df_keep$IID, covar$IID)
  df_keep_wcovar = cbind(df_keep, covar[,3:14])

  ## Check AUC
  print(bigstatsr::AUC(df_keep$hapmap_prs, df_keep$MDD))
  print(bigstatsr::AUC(df_keep$cr95set_prs, df_keep$MDD))
  print(bigstatsr::AUC(df_keep$logbf_over4_prs, df_keep$MDD))
  print(bigstatsr::AUC(df_keep$twassig_prs, df_keep$MDD))


  ## Save data
  fwrite(df_all, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_allprediction.txt"))
  fwrite(df_keep, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_QCset_prediction_wMDDstatus.txt"))
  fwrite(df_keep_wcovar, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_QCset_prediction_wMDDstatus_wcovar.txt"))


}
