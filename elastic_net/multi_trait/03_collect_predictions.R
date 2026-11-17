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


dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/multi_trait/")

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
  for(trait in c("BIP", "SCZ", "ADHD", "ANXmetasubtract")){
    for(group in c("hapmap_cr95set", "hapmap_logbf_over4", "hapamap")){
      file =  paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/"), trait,"_1kg_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_prediction_split_",i,"_updated.txt", sep = "")
      pred = fread(file)

      if(trait == "ANXmetasubtract"){
        trait_name = "ANX"
      } else {
        trait_name = trait
      }

      if(group == "hapmap_cr95set"){
        group_name = "cr95set"
      }
      if(group == "hapmap_logbf_over4"){
        group_name = "logbf4"
      }
      if(group == "hapamap"){
        group_name = "hapmap"
      }
      colnames(pred) = paste(trait_name, group_name, sep = "_")
      assign(paste(trait_name, group_name, sep = "_"),pred)
    }
  }


  ## Load data
  for( type in c("cr95setlogbf4_ALLCovar" ,  "cr95setonly_ALLCovar",  "cr95setlogbf4_NoCovar" ,  "cr95setonly_NoCovar") ){
    file = paste(dir, "/results/LDref1kg_split", i, "_predictortype-",type,"_pred_covar0_multitrait.csv", sep = "")
    pred = fread(file)
    colnames(pred) = paste(type, colnames(pred), sep = "_")
    assign(paste(type),pred)
  }


  ## Check same number of rows
  if(nrow(BIP_hapmap) != length(y01.test) | nrow(BIP_cr95set)  != length(y01.test) |  nrow(BIP_logbf4)  != length(y01.test) |
     nrow(SCZ_hapmap) != length(y01.test) | nrow(SCZ_cr95set)  != length(y01.test) |  nrow(SCZ_logbf4)  != length(y01.test) |
     nrow(ADHD_hapmap) != length(y01.test) | nrow(ADHD_cr95set)  != length(y01.test) |  nrow(ADHD_logbf4)  != length(y01.test) |
     nrow(ANX_hapmap) != length(y01.test) | nrow(ANX_cr95set)  != length(y01.test) |  nrow(ANX_logbf4)  != length(y01.test) |
      nrow(cr95setlogbf4_ALLCovar) != length(y01.test) | nrow(cr95setonly_ALLCovar)  != length(y01.test) |
     nrow(cr95setlogbf4_NoCovar)  != length(y01.test) |  nrow(cr95setonly_NoCovar)  != length(y01.test) |  nrow(covar.test)  != length(y01.test) ){
    stop("Not all same rows")
  }

  ## Combine
  cov_res = cbind(y01.test,
                  BIP_hapmap, BIP_cr95set, BIP_logbf4,
                  SCZ_hapmap, SCZ_cr95set, SCZ_logbf4,
                  ADHD_hapmap, ADHD_cr95set, ADHD_logbf4,
                  ANX_hapmap, ANX_cr95set, ANX_logbf4,
                  cr95setlogbf4_ALLCovar, cr95setonly_ALLCovar,
                  cr95setlogbf4_NoCovar, cr95setonly_NoCovar,
                  covar.test)


  ## Save
  fwrite(cov_res, paste(dir,"/results/split",i,"_testset_covar0_pred_multitrait.csv" ,sep = ""))

}


###############################################################################################

############################## WITH COVAR ######################################################

###############################################################################################

dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/multi_trait/")

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
  for(trait in c("BIP", "SCZ", "ADHD", "ANXmetasubtract")){
    for(group in c("hapmap_cr95set", "hapmap_logbf_over4", "hapamap")){
      file =  paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/"), trait,"_1kg_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_prediction_split_",i,"_updated.txt", sep = "")
      pred = fread(file)

      if(trait == "ANXmetasubtract"){
        trait_name = "ANX"
      } else {
        trait_name = trait
      }

      if(group == "hapmap_cr95set"){
        group_name = "cr95set"
      }
      if(group == "hapmap_logbf_over4"){
        group_name = "logbf4"
      }
      if(group == "hapamap"){
        group_name = "hapmap"
      }
      colnames(pred) = paste(trait_name, group_name, sep = "_")
      assign(paste(trait_name, group_name, sep = "_"),pred)
    }
  }


  ## Load data
  for( type in c("cr95setlogbf4_ALLCovar" ,  "cr95setonly_ALLCovar",  "cr95setlogbf4_NoCovar" ,  "cr95setonly_NoCovar") ){
    file = paste(dir, "/results/LDref1kg_split", i, "_predictortype-",type,"_res_multitrait.txt", sep = "")
    pred = fread(file)
    pred = pred[,c(2:5,7,9,11)]
    colnames(pred) = c("lasso", "enet_alphagrid", "enet_0.5", "ridge", "caret_kappa", "caret_deviance", "caret_roc")
    colnames(pred) = paste(type, colnames(pred), sep = "_")
    assign(paste(type),pred)
  }


  ## Check same number of rows
  if(nrow(BIP_hapmap) != length(y01.test) | nrow(BIP_cr95set)  != length(y01.test) |  nrow(BIP_logbf4)  != length(y01.test) |
     nrow(SCZ_hapmap) != length(y01.test) | nrow(SCZ_cr95set)  != length(y01.test) |  nrow(SCZ_logbf4)  != length(y01.test) |
     nrow(ADHD_hapmap) != length(y01.test) | nrow(ADHD_cr95set)  != length(y01.test) |  nrow(ADHD_logbf4)  != length(y01.test) |
     nrow(ANX_hapmap) != length(y01.test) | nrow(ANX_cr95set)  != length(y01.test) |  nrow(ANX_logbf4)  != length(y01.test) |
     nrow(cr95setlogbf4_ALLCovar) != length(y01.test) | nrow(cr95setonly_ALLCovar)  != length(y01.test) |
     nrow(cr95setlogbf4_NoCovar)  != length(y01.test) |  nrow(cr95setonly_NoCovar)  != length(y01.test) |  nrow(covar.test)  != length(y01.test) ){
    stop("Not all same rows")
  }

  ## Combine
  cov_res = cbind(y01.test,
                  BIP_hapmap, BIP_cr95set, BIP_logbf4,
                  SCZ_hapmap, SCZ_cr95set, SCZ_logbf4,
                  ADHD_hapmap, ADHD_cr95set, ADHD_logbf4,
                  ANX_hapmap, ANX_cr95set, ANX_logbf4,
                  cr95setlogbf4_ALLCovar, cr95setonly_ALLCovar,
                  cr95setlogbf4_NoCovar, cr95setonly_NoCovar,
                  covar.test)


  ## Save
  fwrite(cov_res, paste(dir,"/results/split",i,"_testset_wcovar_pred_multitrait.csv" ,sep = ""))

}

for( i in 1:10){

  f = paste(dir,"/results/split",i,"_testset_wcovar_pred_multitrait.csv" ,sep = "")
  f1 = paste(dir,"/results/split",i,"_testset_covar0_pred_multitrait.csv" ,sep = "")

  z = fread(f) %>% as.data.frame()
  z1 = fread(f1)  %>% as.data.frame()

  colnames(z) = sub("enet_0.5", "enet0.5", colnames(z))

  z2 = z1[,match(colnames(z), colnames(z1))]

  max_diff = apply(z2 - z, 2, function(x) max(abs(x)))
  print(max(max_diff[!grepl("ALLCovar",names(max_diff))]))

  }
