# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

##########################################################################################################################################

#################################################### MULTITRAIT ##########################################################################

##########################################################################################################################################

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

  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/BIP_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_cr95set_prs_cs = read.table(BIP_prscs_file)
  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/BIP_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_logbf4_prs_cs = read.table(BIP_prscs_file)
  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/BIP_"),LD_ref,"_LD_ref_hapamap_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_hapmap_prs_cs = read.table(BIP_prscs_file)

  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/SCZ_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_cr95set_prs_cs = read.table(SCZ_prscs_file)
  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/SCZ_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_logbf4_prs_cs = read.table(SCZ_prscs_file)
  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/SCZ_"),LD_ref,"_LD_ref_hapamap_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_hapmap_prs_cs = read.table(SCZ_prscs_file)

  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ADHD_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_cr95set_prs_cs = read.table(ADHD_prscs_file)
  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ADHD_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_logbf4_prs_cs = read.table(ADHD_prscs_file)
  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ADHD_"),LD_ref,"_LD_ref_hapamap_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_hapmap_prs_cs = read.table(ADHD_prscs_file)


  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ANXmetasubtract_"),LD_ref,"_LD_ref_hapmap_cr95set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_cr95set_prs_cs = read.table(ANX_prscs_file)
  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ANXmetasubtract_"),LD_ref,"_LD_ref_hapmap_logbf_over4_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_logbf4_prs_cs = read.table(ANX_prscs_file)
  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ANXmetasubtract_"),LD_ref,"_LD_ref_hapamap_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_hapmap_prs_cs = read.table(ANX_prscs_file)


  ## Create dataframe
  df_all = data.frame(FID = sample$ID_1, IID = sample$ID_2,
                      BIP_cr95set_prs = BIP_cr95set_prs_cs$x,
                      BIP_logbf4_prs = BIP_logbf4_prs_cs$x,
                      BIP_hapmap_prs = BIP_hapmap_prs_cs$x,


                      ANX_cr95set_prs = ANX_cr95set_prs_cs$x,
                      ANX_logbf4_prs = ANX_logbf4_prs_cs$x,
                      ANX_hapmap_prs = ANX_hapmap_prs_cs$x,


                      SCZ_cr95set_prs = SCZ_cr95set_prs_cs$x,
                      SCZ_logbf4_prs = SCZ_logbf4_prs_cs$x,
                      SCZ_hapmap_prs = SCZ_hapmap_prs_cs$x,


                      ADHD_cr95set_prs = ADHD_cr95set_prs_cs$x,
                      ADHD_logbf4_prs = ADHD_logbf4_prs_cs$x,
                      ADHD_hapmap_prs = ADHD_hapmap_prs_cs$x)


  ## Create object ID
  object_id = paste(sample$ID_1, sample$ID_2, sep = "_")


  ## Load other ancestries

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

  ## Add predictions to all
  if(!
     all.equal(object_id[match(EUR_qc_unique_id ,object_id)], EUR_qc_unique_id )  ){
    stop("Not all ID equal - check")
  }

  ## Add info
  ## Combine prediction info, MDD status, and covariate info
  df_EUR = cbind(df_all[match(EUR_qc_unique_id,object_id),], EUR_QC_all[,3:15])
  ## Add covar info
  all.equal(df_EUR$FID, EUR_QC_all$FID)


  ## Save data
  fwrite(df_EUR, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_EUR_clean_prediction_multitrait.txt"))

}


##########################################################################################################################################

#################################################### AltTrait Set #####################################################################

##########################################################################################################################################

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

  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_BIP_cr95set_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_cr95set_prs_cs = read.table(BIP_prscs_file)
  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_BIP_logbfover4_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_logbf4_prs_cs = read.table(BIP_prscs_file)

  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_ADHD_cr95set_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_cr95set_prs_cs = read.table(ADHD_prscs_file)
  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_ADHD_logbfover4_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_logbf4_prs_cs = read.table(ADHD_prscs_file)

  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_SCZ_cr95set_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_cr95set_prs_cs = read.table(SCZ_prscs_file)
  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_SCZ_logbfover4_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_logbf4_prs_cs = read.table(SCZ_prscs_file)


  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_ANXmetasubtract_cr95set_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_cr95set_prs_cs = read.table(ANX_prscs_file)
  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_ANXmetasubtract_logbfover4_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_logbf4_prs_cs = read.table(ANX_prscs_file)

  ## Create dataframe
  df_all = data.frame(FID = sample$ID_1, IID = sample$ID_2,
                      BIP_cr95set_prs = BIP_cr95set_prs_cs$x,
                      BIP_logbf4_prs = BIP_logbf4_prs_cs$x,


                      ANX_cr95set_prs = ANX_cr95set_prs_cs$x,
                      ANX_logbf4_prs = ANX_logbf4_prs_cs$x,


                      SCZ_cr95set_prs = SCZ_cr95set_prs_cs$x,
                      SCZ_logbf4_prs = SCZ_logbf4_prs_cs$x,


                      ADHD_cr95set_prs = ADHD_cr95set_prs_cs$x,
                      ADHD_logbf4_prs = ADHD_logbf4_prs_cs$x)


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

  ## Add predictions to all
  if(!
     all.equal(object_id[match(EUR_qc_unique_id ,object_id)], EUR_qc_unique_id )  ){
    stop("Not all ID equal - check")
  }

  ## Add info
  ## Combine prediction info, MDD status, and covariate info
  df_EUR = cbind(df_all[match(EUR_qc_unique_id,object_id),], EUR_QC_all[,3:15])
  ## Add covar info
  all.equal(df_EUR$FID, EUR_QC_all$FID)


  ## Save data
  fwrite(df_EUR, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_EUR_clean_prediction_alttraitset.txt"))


}


##########################################################################################################################################

################################################# ALTTRAIT COEF ##########################################################################

##########################################################################################################################################

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

  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_BIPset_"),LD_ref,"_LD_ref_cr95set_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_cr95set_prs_cs = read.table(BIP_prscs_file)
  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_BIPset_"),LD_ref,"_LD_ref_logbfover4_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_logbf4_prs_cs = read.table(BIP_prscs_file)
  BIP_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_BIPset_"),LD_ref,"_LD_ref_twassig_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  BIP_twassig_prs_cs = read.table(BIP_prscs_file)


  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_ANXmetasubtractset_"),LD_ref,"_LD_ref_cr95set_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_cr95set_prs_cs = read.table(ANX_prscs_file)
  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_ANXmetasubtractset_"),LD_ref,"_LD_ref_logbfover4_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_logbf4_prs_cs = read.table(ANX_prscs_file)
  ANX_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_ANXmetasubtractset_"),LD_ref,"_LD_ref_twassig_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ANX_twassig_prs_cs = read.table(ANX_prscs_file)

  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_ADHDset_"),LD_ref,"_LD_ref_cr95set_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_cr95set_prs_cs = read.table(ADHD_prscs_file)
  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_ADHDset_"),LD_ref,"_LD_ref_logbfover4_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_logbf4_prs_cs = read.table(ADHD_prscs_file)
  ADHD_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_ADHDset_"),LD_ref,"_LD_ref_twassig_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  ADHD_twassig_prs_cs = read.table(ADHD_prscs_file)


  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_SCZset_"),LD_ref,"_LD_ref_cr95set_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_cr95set_prs_cs = read.table(SCZ_prscs_file)
  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_SCZset_"),LD_ref,"_LD_ref_logbfover4_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_logbf4_prs_cs = read.table(SCZ_prscs_file)
  SCZ_prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_SCZset_"),LD_ref,"_LD_ref_twassig_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
  SCZ_twassig_prs_cs = read.table(SCZ_prscs_file)

  ## Create dataframe
  df_all = data.frame(FID = sample$ID_1, IID = sample$ID_2,
                      BIP_cr95set_prs = BIP_cr95set_prs_cs$x,
                      BIP_logbf4_prs = BIP_logbf4_prs_cs$x,
                      BIP_twassig_prs = BIP_twassig_prs_cs$x,


                      ANX_cr95set_prs = ANX_cr95set_prs_cs$x,
                      ANX_logbf4_prs = ANX_logbf4_prs_cs$x,
                      ANX_twassig_prs = ANX_twassig_prs_cs$x,


                      SCZ_cr95set_prs = SCZ_cr95set_prs_cs$x,
                      SCZ_logbf4_prs = SCZ_logbf4_prs_cs$x,
                      SCZ_twassig_prs = SCZ_twassig_prs_cs$x,


                      ADHD_cr95set_prs = ADHD_cr95set_prs_cs$x,
                      ADHD_logbf4_prs = ADHD_logbf4_prs_cs$x,
                      ADHD_twassig_prs = ADHD_twassig_prs_cs$x)


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

  ## Add predictions to all
  if(!
     all.equal(object_id[match(EUR_qc_unique_id ,object_id)], EUR_qc_unique_id )  ){
    stop("Not all ID equal - check")
  }

  ## Add info
  ## Combine prediction info, MDD status, and covariate info
  df_EUR = cbind(df_all[match(EUR_qc_unique_id,object_id),], EUR_QC_all[,3:15])
  ## Add covar info
  all.equal(df_EUR$FID, EUR_QC_all$FID)


  ## Save data
  fwrite(df_EUR, paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_EUR_clean_prediction_alttraitcoef.txt"))


}
