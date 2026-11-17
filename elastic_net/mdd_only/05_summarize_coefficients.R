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
out_dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/single_trait/")


df_all = data.frame()
for(i in 1:10){
  for(type in c( "ALLCovar_wTWAS",  "ALLCovar_NoTWAS",  "NoCovar_wTWAS",  "NoCovar_NoTWAS")){

    ## Columns keep for each type
    if(type == "ALLCovar_wTWAS"){
      print(type)
      index_coef = 2:4
    }
    if(type == "ALLCovar_NoTWAS"){
      print(type)
      index_coef = 2:3
    }
    if(type == "NoCovar_wTWAS"){
      print(type)
      index_coef = 2:4
    }
    if(type == "NoCovar_NoTWAS"){
      print(type)
      index_coef = 2:3
    }

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
    print(coef_full[index_coef,])

    ## Calculate number nonzero
    number_nonzero_coef = colSums(coef_full[index_coef,]!=0)
    number_nonzero_coef_logbf4_twassig = colSums(coef_full[grepl("logbf_over4", rownames(coef_full))|grepl("twassig", rownames(coef_full)),]!=0)
    number_nonzero_coef_logbf4 = colSums(coef_full[grepl("logbf_over4", rownames(coef_full)),]!=0)

    ## Combine into dataframe
    df_cur = data.frame(split= i, type = type, model_type = names(number_nonzero_coef), number_nonzero_coef_logbf4_twassig,
                        number_nonzero_coef_logbf4, number_nonzero_coef)
    df_all = rbind(df_all, df_cur)
  }
}


write.csv(df_all,
          paste(out_dir,"/results/ncoef_summary.csv", sep = ""),
          row.names = F
)

ggplot(df_all, aes(x = model_type, y = number_nonzero_coef))+facet_wrap(~type)+geom_boxplot()
df_all[df_all$number_nonzero_coef > 1,]
