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
  for(trait in c("BIP", "SCZ", "ADHD",  "ANXmetasubtract")){
    for(group in c("hapmap_cr95set", "hapmap_logbf_over4", "hapamap")){
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


        prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/"), trait,"_",LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
        prs_cs = read.table(prscs_file)

        if(!all.equal(nrow(sample), nrow(prs_cs))){
          print("something went wrong - diff number of samples")
          stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
        }

        prs_cs_test = prs_cs$x[ind.test]

        testset_predfile = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/"), trait,"_",LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_prediction_split_",i,"_updated.txt", sep = "")
        write.table(prs_cs_test, testset_predfile, row.names = F)

        auc_val = AUC(prs_cs_test, y01.test)
        neg_auc_val = 1-AUC(prs_cs_test, y01.test)
        df_run = data.frame(split = i, LD_ref = LD_ref, test_auc = neg_auc_val, n_test_auc = auc_val, group = group, trait = trait)
        df_total = rbind(df_total, df_run)

      }
    }
  }
}

out_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs_cs_results_alttrait_altgwas_updated.txt"), sep = "")
write.table(df_total, out_file)
df_total

################################################### Alt Trait PRS Sets ######################################################

## Arguments Loaded in
df_total = data.frame()
for(LD_ref in c( "1kg")){
  for(trait in c("BIP", "SCZ", "ADHD",  "ANXmetasubtract")){
    for(group in c("cr95set", "logbfover4")){
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


        prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_", trait,"_",group,"_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")

        prs_cs = read.table(prscs_file)

        if(!all.equal(nrow(sample), nrow(prs_cs))){
          print("something went wrong - diff number of samples")
          stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
        }

        prs_cs_test = prs_cs$x[ind.test]

        testset_predfile = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/_MDD_GWAS_"),LD_ref,"_LD_ref_", trait, "_",group,"_set_pst_eff_a1_b0.5_phiauto_prediction_split_",i,"_updated.txt", sep = "")
        write.table(prs_cs_test, testset_predfile, row.names = F)

        auc_val = AUC(prs_cs_test, y01.test)
        neg_auc_val = 1-AUC(prs_cs_test, y01.test)
        df_run = data.frame(split = i, LD_ref = LD_ref, test_auc = neg_auc_val, n_test_auc = auc_val, group = group, trait = trait)
        df_total = rbind(df_total, df_run)

      }
    }
  }
}

out_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs_cs_results_alttraitsets_MDDgwas_updated.txt"), sep = "")
write.table(df_total, out_file)
df_total


################################################### Alt Trait PRS Coef MDD Sets ######################################################

## Arguments Loaded in
df_total = data.frame()
for(LD_ref in c( "1kg")){
  for(trait in c("BIP", "SCZ", "ADHD",  "ANXmetasubtract")){
    for(group in c("cr95set", "logbfover4", "twassig")){
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


        prscs_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),trait,"set_",LD_ref,"_LD_ref_",group,"_coef_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")

        prs_cs = read.table(prscs_file)

        if(!all.equal(nrow(sample), nrow(prs_cs))){
          print("something went wrong - diff number of samples")
          stop("Input identity or dimension check failed; inspect the preceding diagnostic.")
        }

        prs_cs_test = prs_cs$x[ind.test]

        testset_predfile = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),trait,"set_",LD_ref,"_LD_ref_",group,"_coef_pst_eff_a1_b0.5_phiauto_prediction_split_",i,"_updated.txt", sep = "")
        write.table(prs_cs_test, testset_predfile, row.names = F)

        auc_val = AUC(prs_cs_test, y01.test)
        neg_auc_val = 1-AUC(prs_cs_test, y01.test)
        df_run = data.frame(split = i, LD_ref = LD_ref, test_auc = neg_auc_val, n_test_auc = auc_val, group = group, trait = trait)
        df_total = rbind(df_total, df_run)

      }
    }
  }
}

out_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs_cs_results_alttraitcoef_MDDsets_updated.txt"), sep = "")
write.table(df_total, out_file)
df_total
