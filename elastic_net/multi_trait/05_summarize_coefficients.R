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
out_dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/multi_trait/")


df_all = data.frame()
for(i in 1:10){
  for(type in c( "cr95setlogbf4_ALLCovar",  "cr95setonly_ALLCovar",   "cr95setlogbf4_NoCovar",  "cr95setonly_NoCovar")){

    ## Columns keep for each type
    if(type == "cr95setlogbf4_ALLCovar"){
      print(type)
      index_coef = 2:11
    }
    if(type ==  "cr95setonly_ALLCovar"){
      print(type)
      index_coef = 2:6
    }
    if(type == "cr95setlogbf4_NoCovar"){
      print(type)
      index_coef =  2:11
    }
    if(type == "cr95setonly_NoCovar"){
      print(type)
      index_coef =  2:6
    }

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
    if(i == 1){
    print(coef_full[index_coef,])
    print(coef_full[grepl("BIP", rownames(coef_full))|
                      grepl("SCZ", rownames(coef_full))|
                      grepl("ADHD", rownames(coef_full))|
                      grepl("ANX", rownames(coef_full)),])
    }


    ## Calculate number nonzero
    number_nonzero_coef = colSums(coef_full[index_coef,]!=0)
    number_nonzero_alttraitcoef = colSums(coef_full[grepl("BIP", rownames(coef_full))|
                                                      grepl("SCZ", rownames(coef_full))|
                                                                  grepl("ADHD", rownames(coef_full))|
                                                                              grepl("ANX", rownames(coef_full)),]!=0)


    ## Combine into dataframe
    df_cur = data.frame(split= i, type = type, model_type = names(number_nonzero_coef), number_nonzero_coef, number_nonzero_alttraitcoef)
    df_all = rbind(df_all, df_cur)
  }
}


write.csv(df_all,
          paste(out_dir,"/results/ncoef_summary_multitrait.csv", sep = ""),
          row.names = F
)

ggplot(df_all, aes(x = model_type, y = number_nonzero_coef))+facet_wrap(~type)+geom_boxplot()
ggplot(df_all, aes(x = model_type, y = number_nonzero_alttraitcoef))+facet_wrap(~type)+geom_boxplot()
df_all[df_all$number_nonzero_coef > 1,]
