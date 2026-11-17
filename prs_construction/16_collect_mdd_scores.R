# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(purrr)
library(bigsnpr)
library(data.table)
library(tidyverse)


## Arguments Loaded in
df_total = data.frame()
for(LD_ref in c( "1kg")){
  for(i in 1:10){

    print(i)


    testphenotype_file = paste(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/splits/split"),i,"_QC_MDD_testing.txt", sep = "")


    ## Load
    testpheno = fread(testphenotype_file)

    ## Obtain Unique ID
    test_id = paste(testpheno$FID, testpheno$IID, sep = "_")


    ## Obtain y train and y test
    y01.test = testpheno$MDD


    ## Match Id's
    sample = fread(mdd_path("MDD_SAMPLE", ""))
    head(sample,3)

    ## remove first individual
    sample = sample[-1,]
    print("remove placeholder")
    head(sample,3)


    object_id = paste(sample$ID_1, sample$ID_2, sep = "_")

    ind.test =  match(test_id, object_id)

    if(!all.equal(object_id[ind.test], test_id)){
      print("something went wrong - test ids not match ind.test")
      stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
    }


    prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_pst_eff_a1_b0.5_phiauto_update_prediction.txt", sep = "")
    prs_cs = read.table(prscs_file)

    if(!all.equal(nrow(sample), nrow(prs_cs))){
      print("something went wrong - diff number of samples")
      stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
    }

    prs_cs_test = prs_cs$x[ind.test]

    testset_predfile = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_pst_eff_a1_b0.5_phiauto_prediction_split_",i,"_redo.txt", sep = "")
    write.table(prs_cs_test, testset_predfile, row.names = F)

    auc_val = bigstatsr::AUC(prs_cs_test, y01.test)
    neg_auc_val = 1-bigstatsr::AUC(prs_cs_test, y01.test)
    df_run = data.frame(split = i, LD_ref = LD_ref, test_auc = neg_auc_val, n_test_auc = auc_val)
    df_total = rbind(df_total, df_run)

  }
}


out_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs_cs_results_initial_redo.txt"), sep = "")
write.table(df_total, out_file)
df_total


################################################################# Repeat for Group by Group #####################################################

library(purrr)
library(bigsnpr)
library(data.table)
library(tidyverse)


## Arguments Loaded in
df_total = data.frame()
for(LD_ref in c("1kg")){
  for(group in c("cr95set", "logbf_over4", "twassig")){
    for(i in 1:10){

      print(i)


      testphenotype_file = paste(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/splits/split"),i,"_QC_MDD_testing.txt", sep = "")


      ## Load
      testpheno = fread(testphenotype_file)

      ## Obtain Unique ID
      test_id = paste(testpheno$FID, testpheno$IID, sep = "_")


      ## Obtain y train and y test
      y01.test = testpheno$MDD


      ## Match Id's
      sample = fread(mdd_path("MDD_SAMPLE", ""))
      head(sample,3)

      ## remove first individual
      sample = sample[-1,]
      print("remove placeholder")
      head(sample,3)


      object_id = paste(sample$ID_1, sample$ID_2, sep = "_")

      ind.test =  match(test_id, object_id)

      if(!all.equal(object_id[ind.test], test_id)){
        print("something went wrong - test ids not match ind.test")
        stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
      }


      prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
      prs_cs = read.table(prscs_file)

      if(!all.equal(nrow(sample), nrow(prs_cs))){
        print("something went wrong - diff number of samples")
        stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
      }

      prs_cs_test = prs_cs$x[ind.test]

      testset_predfile = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_prediction_split_",i,".txt", sep = "")
      write.table(prs_cs_test, testset_predfile, row.names = F)

      auc_val = bigstatsr::AUC(prs_cs_test, y01.test)
      neg_auc_val = 1-bigstatsr::AUC(prs_cs_test, y01.test)
      df_run = data.frame(split = i, LD_ref = LD_ref, test_auc = neg_auc_val, n_test_auc = auc_val, group = group)
      df_total = rbind(df_total, df_run)

    }
  }
}

out_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs_cs_results_initial_bygroup_redo.txt"), sep = "")
write.table(df_total, out_file)
df_total
