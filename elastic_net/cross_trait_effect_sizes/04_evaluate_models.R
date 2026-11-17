# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(data.table)
library(tidyverse)
library(fmsb)
library(bigsnpr)


##########################################################################################################################################

############################################# PROCESS RESULTS COVAR 0 ####################################################################

##########################################################################################################################################

dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/muti_trait_MDDset_alttraitcoef/")

res = data.frame()
for(i in 1:10){
  print(i)

  ## Load in model
  pred_wcov = fread(paste(dir,"/results/split",i,"_testset_covar0_pred_alttraitcoef.csv" ,sep = "")) %>% as.data.frame()
  pred_wcov[,c(2:49,51)] = scale(pred_wcov[,c(2:49,51)])

  pred_wcov$assessment_center = as.factor(pred_wcov$assessment_center)

  ## Rename y01.test as MDD
  pred_wcov$MDD = pred_wcov$y01.test

  ## Fit null model
  null_mod = glm(MDD ~ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+
                   sex+age+as.factor(assessment_center)+genotype_axiom,family = "binomial",data = pred_wcov)

  ## Adjusted phenotype
  pheno_adj = null_mod$residuals
  pheno_adj_norm = qnorm((rank(pheno_adj ,na.last="keep")-0.5)/sum(!is.na(pheno_adj )))

  df = data.frame()
  for(k in 2:41){
    model_type = colnames(pred_wcov)[k]
    mod = try(glm(MDD ~ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+
                    sex+age+as.factor(assessment_center)+genotype_axiom+pred_wcov[,k],family = "binomial",data = pred_wcov))
    if(!"try-error" %in% class(mod)){
      coef_output <- summary(mod)$coefficients[rownames(summary(mod)$coefficients) == "pred_wcov[, k]", c(1:4)]
      AUC_fullmodel_wcovar = AUC(mod$fitted.values,pred_wcov$MDD)
      r2 =NagelkerkeR2(mod)$R2
      rr <- mod
      full_r2 <-  (1 - exp((rr$dev - rr$null)/nrow(rr$model)))/(1 - exp(-rr$null/nrow(rr$model)))
      AUC = AUC(pred_wcov[,k],pred_wcov$MDD)
      pred_p_quant_10 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.10))),
                                labels=c("0-10", "10-20","20-30",
                                         "30-40","40-50", "50-60",
                                         "60-70","70-80","80-90",
                                         "90-100"), include.lowest=TRUE))
          pred_p_quant_10_mod = try(cut(mod$fitted.values, breaks=c(quantile(mod$fitted.values, probs = seq(0, 1, by = 0.10))),
                                labels=c("0-10", "10-20","20-30",
                                         "30-40","40-50", "50-60",
                                         "60-70","70-80","80-90",
                                         "90-100"), include.lowest=TRUE))
      pred_p_quant_100_mod = try(cut(mod$fitted.values, breaks=c(quantile(mod$fitted.values, probs = seq(0, 1, by = 0.01))),
                                 labels=c(1:100), include.lowest=TRUE))
      pred_p_quant_20_mod = try(cut(mod$fitted.values, breaks=c(quantile(mod$fitted.values, probs = seq(0, 1, by = 0.05))),
                                labels=c(1:20), include.lowest=TRUE))

      pred_p_quant_50_mod = try(cut(mod$fitted.values, breaks=c(quantile(mod$fitted.values, probs = seq(0, 1, by = 0.02))),
                                labels=c(1:50), include.lowest=TRUE))

      pred_p_quant_10_table_mod = table(pred_p_quant_10_mod, pred_wcov$MDD)
      pred_p_quant_100_table_mod = table(pred_p_quant_100_mod, pred_wcov$MDD)
      pred_p_quant_20_table_mod = table(pred_p_quant_20_mod, pred_wcov$MDD)
      pred_p_quant_50_table_mod = table(pred_p_quant_50_mod, pred_wcov$MDD)
      pred_p_quant_10_table_top_odds_mod = pred_p_quant_10_table_mod[rownames(pred_p_quant_10_table_mod) == "90-100",][2]/pred_p_quant_10_table_mod[rownames(pred_p_quant_10_table_mod) == "90-100",][1]
      pred_p_quant_10_table_nottop_odds_mod = colSums(pred_p_quant_10_table_mod[rownames(pred_p_quant_10_table_mod) != "90-100",])[2]/colSums(pred_p_quant_10_table_mod[rownames(pred_p_quant_10_table_mod) != "90-100",])[1]

      pred_p_quant_100_table_top_odds_mod = pred_p_quant_100_table_mod[rownames(pred_p_quant_100_table_mod) == "100",][2]/pred_p_quant_100_table_mod[rownames(pred_p_quant_100_table_mod) == "100",][1]
      pred_p_quant_100_table_nottop_odds_mod = colSums(pred_p_quant_100_table_mod[rownames(pred_p_quant_100_table_mod) != "100",])[2]/colSums(pred_p_quant_100_table_mod[rownames(pred_p_quant_100_table_mod) != "100",])[1]

      pred_p_quant_20_table_top_odds_mod = pred_p_quant_20_table_mod[rownames(pred_p_quant_20_table_mod) == "20",][2]/pred_p_quant_20_table_mod[rownames(pred_p_quant_20_table_mod) == "20",][1]
      pred_p_quant_20_table_nottop_odds_mod = colSums(pred_p_quant_20_table_mod[rownames(pred_p_quant_20_table_mod) != "20",])[2]/colSums(pred_p_quant_20_table_mod[rownames(pred_p_quant_20_table_mod) != "20",])[1]

      pred_p_quant_50_table_top_odds_mod = pred_p_quant_50_table_mod[rownames(pred_p_quant_50_table_mod) == "50",][2]/pred_p_quant_50_table_mod[rownames(pred_p_quant_50_table_mod) == "50",][1]
      pred_p_quant_50_table_nottop_odds_mod = colSums(pred_p_quant_50_table_mod[rownames(pred_p_quant_50_table_mod) != "50",])[2]/colSums(pred_p_quant_50_table_mod[rownames(pred_p_quant_50_table_mod) != "50",])[1]

      top10_rest_mod = pred_p_quant_10_table_top_odds_mod/pred_p_quant_10_table_nottop_odds_mod
      top5_rest_mod = pred_p_quant_20_table_top_odds_mod/pred_p_quant_20_table_nottop_odds_mod
      top2_rest_mod = pred_p_quant_50_table_top_odds_mod/pred_p_quant_50_table_nottop_odds_mod
      top1_rest_mod = pred_p_quant_100_table_top_odds_mod/pred_p_quant_100_table_nottop_odds_mod

      if("try-error" %in% class(pred_p_quant_10)){

        if(AUC<0.5){
          AUC = 1 - AUC}
        cor_pheno_adj = cor(pheno_adj, pred_wcov[,k])
        decile_top_bottom_or =  NA
        decile_top_bottom_prev = NA
        decile_top_mid_or =  NA
        decile_top_mid_prev = NA
        quintile_top_bottom_or = NA
        quintile_top_bottom_prev =NA
        quintile_top_mid_or =  NA
        quintile_top_mid_prev = NA

        top10_rest = NA
        top5_rest = NA
        top2_rest = NA
        top1_rest = NA
      } else{

        pred_p_quant_100 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.01))),
                                   labels=c(1:100), include.lowest=TRUE))
        pred_p_quant_20 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.05))),
                                  labels=c(1:20), include.lowest=TRUE))

        pred_p_quant_50 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.02))),
                                  labels=c(1:50), include.lowest=TRUE))

        pred_p_quant_5 = cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.20))),
                             labels=c("0-20","20-40",
                                      "40-60",
                                      "60-80","80-100"), include.lowest=TRUE)

        pred_p_quant_10_table = table(pred_p_quant_10, pred_wcov$MDD)
        pred_p_quant_10_table_top_odds = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][2]/pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][1]
        pred_p_quant_10_table_top_prev = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][2]/(pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][1] + pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][2])
        pred_p_quant_10_table_bottom_odds = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][2]/pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][1]
        pred_p_quant_10_table_bottom_prev = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][2]/(pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][1] + pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][2])
        pred_p_quant_10_table_mid_odds = sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),][,2])/sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),][,1])
        pred_p_quant_10_table_mid_prev = sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),][,2])/sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),])
        pred_p_quant_10_table_nottop_odds = colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "90-100",])[2]/colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "90-100",])[1]
        pred_p_quant_10_table_notbottom_odds = colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "0-10",])[2]/colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "0-10",])[1]


        pred_p_quant_5_table = table(pred_p_quant_5, pred_wcov$MDD)
        pred_p_quant_5_table_top_odds = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][2]/pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][1]
        pred_p_quant_5_table_top_prev = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][2]/(pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][1] + pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][2])
        pred_p_quant_5_table_bottom_odds = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][2]/pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][1]
        pred_p_quant_5_table_bottom_prev = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][2]/(pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][1] + pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][2])
        pred_p_quant_5_table_mid_odds = pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][2]/pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][1]
        pred_p_quant_5_table_mid_prev = pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][2]/(pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][1] + pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][2])

        pred_p_quant_100_table = table(pred_p_quant_100, pred_wcov$MDD)
        pred_p_quant_100_table_top_odds = pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "100",][2]/ pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "100",][1]
        pred_p_quant_100_table_bottom_odds = pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "1",][2]/ pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "1",][1]
        pred_p_quant_100_table_nottop_odds = colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "100",])[2]/colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "100",])[1]
        pred_p_quant_100_table_notbottom_odds = colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "1",])[2]/colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "1",])[1]

        pred_p_quant_20_table = table(pred_p_quant_20, pred_wcov$MDD)
        pred_p_quant_20_table_top_odds = pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "20",][2]/ pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "20",][1]
        pred_p_quant_20_table_bottom_odds = pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "1",][2]/ pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "1",][1]
        pred_p_quant_20_table_nottop_odds = colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "20",])[2]/colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "20",])[1]
        pred_p_quant_20_table_notbottom_odds = colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "1",])[2]/colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "1",])[1]

        pred_p_quant_50_table = table(pred_p_quant_50, pred_wcov$MDD)
        pred_p_quant_50_table_top_odds = pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "50",][2]/ pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "50",][1]
        pred_p_quant_50_table_bottom_odds = pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "1",][2]/ pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "1",][1]
        pred_p_quant_50_table_nottop_odds = colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "50",])[2]/colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "50",])[1]
        pred_p_quant_50_table_notbottom_odds = colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "1",])[2]/colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "1",])[1]


        decile_top_bottom_or = pred_p_quant_10_table_top_odds/pred_p_quant_10_table_bottom_odds
        decile_top_bottom_prev = pred_p_quant_10_table_top_prev/pred_p_quant_10_table_bottom_prev

        decile_top_mid_or = pred_p_quant_10_table_top_odds/pred_p_quant_10_table_mid_odds
        decile_top_mid_prev = pred_p_quant_10_table_top_prev/pred_p_quant_10_table_mid_prev

        top10_rest = pred_p_quant_10_table_top_odds/pred_p_quant_10_table_nottop_odds
        top5_rest = pred_p_quant_20_table_top_odds/pred_p_quant_20_table_nottop_odds
        top2_rest = pred_p_quant_50_table_top_odds/pred_p_quant_50_table_nottop_odds
        top1_rest = pred_p_quant_100_table_top_odds/pred_p_quant_100_table_nottop_odds

        quintile_top_bottom_or = pred_p_quant_5_table_top_odds/pred_p_quant_5_table_bottom_odds
        quintile_top_bottom_prev = pred_p_quant_5_table_top_prev/pred_p_quant_5_table_bottom_prev

        quintile_top_mid_or = pred_p_quant_5_table_top_odds/pred_p_quant_5_table_mid_odds
        quintile_top_mid_prev = pred_p_quant_5_table_top_prev/pred_p_quant_5_table_mid_prev
        cor_pheno_adj = cor(pheno_adj, pred_wcov[,k])
        if(AUC<0.5){
          AUC = 1 - AUC

          decile_top_bottom_or = pred_p_quant_10_table_bottom_odds/pred_p_quant_10_table_top_odds
          decile_top_bottom_prev = pred_p_quant_10_table_bottom_prev/pred_p_quant_10_table_top_prev

          decile_top_mid_or = pred_p_quant_10_table_bottom_odds/pred_p_quant_10_table_mid_odds
          decile_top_mid_prev = pred_p_quant_10_table_bottom_prev/pred_p_quant_10_table_mid_prev

          quintile_top_bottom_or = pred_p_quant_5_table_bottom_odds/pred_p_quant_5_table_top_odds
          quintile_top_bottom_prev = pred_p_quant_5_table_bottom_prev/pred_p_quant_5_table_top_prev

          quintile_top_mid_or = pred_p_quant_5_table_bottom_odds/pred_p_quant_5_table_mid_odds
          quintile_top_mid_prev = pred_p_quant_5_table_bottom_prev/pred_p_quant_5_table_mid_prev

          cor_pheno_adj = -cor(pheno_adj, pred_wcov[,k])


          top10_rest = pred_p_quant_10_table_bottom_odds/pred_p_quant_10_table_notbottom_odds
          top5_rest = pred_p_quant_20_table_bottom_odds/pred_p_quant_20_table_notbottom_odds
          top2_rest = pred_p_quant_50_table_bottom_odds/pred_p_quant_50_table_notbottom_odds
          top1_rest = pred_p_quant_100_table_bottom_odds/pred_p_quant_100_table_notbottom_odds


        }

      }

      df_out = data.frame(split = i, model = model_type, nag_r2 = r2, test_auc = AUC,
                          cor_pheno_adj = cor_pheno_adj, decile_top_bottom_or =  decile_top_bottom_or,
                          decile_top_bottom_prev = decile_top_bottom_prev,

                          decile_top_mid_or =  decile_top_mid_or  ,
                          decile_top_mid_prev =  decile_top_mid_prev  ,

                          quintile_top_bottom_or = quintile_top_bottom_or,
                          quintile_top_bottom_prev = quintile_top_bottom_prev,

                          quintile_top_mid_or =  quintile_top_mid_or ,
                          quintile_top_mid_prev = quintile_top_mid_prev,

                          top10_rest,
                          top5_rest,
                          top2_rest,
                          top1_rest,
                          full_r2 = full_r2,
                          AUC_fullmodel_wcovar = AUC_fullmodel_wcovar,
                          beta = (coef_output)[1],
                          std_err = (coef_output)[2],
                          zscore = (coef_output)[3],
                          pval = (coef_output)[4],
                          OR=exp(as.numeric((coef_output)[1])),
                          lowerCI=exp(as.numeric((coef_output)[1])-1.96*as.numeric((coef_output)[2])),
                          upperCI=exp(as.numeric((coef_output)[1])+1.96*as.numeric((coef_output)[2])),
                          top10_rest_mod,
                          top5_rest_mod,
                          top2_rest_mod,
                          top1_rest_mod
      )

      df = rbind(df, df_out)

    }

  }

  res = rbind(res, df)

}


fwrite(res, paste0(dir,"/results/combined_results_covar0_alttraitcoef_071125.csv"))
head(res)


##########################################################################################################################################

############################################# PROCESS RESULTS W COVAR ####################################################################

##########################################################################################################################################


dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/muti_trait_MDDset_alttraitcoef/")


res = data.frame()
for(i in 1:10){
  print(i)

  ## Load in model
  pred_wcov = fread(paste(dir,"/results/split",i,"_testset_wcovar_pred_alttraitcoef.csv" ,sep = "")) %>% as.data.frame()
  pred_wcov[,c(2:49,51)] = scale(pred_wcov[,c(2:49,51)])

  pred_wcov$assessment_center = as.factor(pred_wcov$assessment_center)

  ## Rename y01.test as MDD
  pred_wcov$MDD = pred_wcov$y01.test
  df = data.frame()


  for(k in 2:41){
    model_type = colnames(pred_wcov)[k]
    mod = try(glm(MDD ~ pred_wcov[,k],family = "binomial",data = pred_wcov))
    if(!"try-error" %in% class(mod)){
      coef_output <- summary(mod)$coefficients[rownames(summary(mod)$coefficients) == "pred_wcov[, k]", c(1:4)]
      AUC_fullmodel_wcovar = AUC(mod$fitted.values,pred_wcov$MDD)
      rr <- mod
      full_r2 <-  (1 - exp((rr$dev - rr$null)/nrow(rr$model)))/(1 - exp(-rr$null/nrow(rr$model)))

      r2 =NagelkerkeR2(mod)$R2
      AUC = AUC(pred_wcov[,k],pred_wcov$MDD)
      pred_p_quant_10 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.10))),
                                labels=c("0-10", "10-20","20-30",
                                         "30-40","40-50", "50-60",
                                         "60-70","70-80","80-90",
                                         "90-100"), include.lowest=TRUE))

      if("try-error" %in% class(pred_p_quant_10)){

        if(AUC<0.5){
          AUC = 1 - AUC}
        decile_top_bottom_or =  NA
        decile_top_bottom_prev = NA
        decile_top_mid_or =  NA
        decile_top_mid_prev = NA
        quintile_top_bottom_or = NA
        quintile_top_bottom_prev =NA
        quintile_top_mid_or =  NA
        quintile_top_mid_prev = NA

        top10_rest = NA
        top5_rest = NA
        top2_rest = NA
        top1_rest = NA
      } else{

        pred_p_quant_100 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.01))),
                                   labels=c(1:100), include.lowest=TRUE))
        pred_p_quant_20 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.05))),
                                  labels=c(1:20), include.lowest=TRUE))

        pred_p_quant_50 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.02))),
                                  labels=c(1:50), include.lowest=TRUE))

        pred_p_quant_5 = cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.20))),
                             labels=c("0-20","20-40",
                                      "40-60",
                                      "60-80","80-100"), include.lowest=TRUE)

        pred_p_quant_10_table = table(pred_p_quant_10, pred_wcov$MDD)
        pred_p_quant_10_table_top_odds = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][2]/pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][1]
        pred_p_quant_10_table_top_prev = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][2]/(pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][1] + pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][2])
        pred_p_quant_10_table_bottom_odds = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][2]/pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][1]
        pred_p_quant_10_table_bottom_prev = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][2]/(pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][1] + pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "0-10",][2])
        pred_p_quant_10_table_mid_odds = sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),][,2])/sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),][,1])
        pred_p_quant_10_table_mid_prev = sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),][,2])/sum(pred_p_quant_10_table[rownames(pred_p_quant_10_table) %in% c("40-50", "50-60"),])
        pred_p_quant_10_table_nottop_odds = colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "90-100",])[2]/colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "90-100",])[1]
        pred_p_quant_10_table_notbottom_odds = colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "0-10",])[2]/colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "0-10",])[1]


        pred_p_quant_5_table = table(pred_p_quant_5, pred_wcov$MDD)
        pred_p_quant_5_table_top_odds = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][2]/pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][1]
        pred_p_quant_5_table_top_prev = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][2]/(pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][1] + pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "80-100",][2])
        pred_p_quant_5_table_bottom_odds = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][2]/pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][1]
        pred_p_quant_5_table_bottom_prev = pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][2]/(pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][1] + pred_p_quant_5_table[rownames(pred_p_quant_5_table) == "0-20",][2])
        pred_p_quant_5_table_mid_odds = pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][2]/pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][1]
        pred_p_quant_5_table_mid_prev = pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][2]/(pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][1] + pred_p_quant_5_table[rownames(pred_p_quant_5_table) %in% c("40-60"),][2])

        pred_p_quant_100_table = table(pred_p_quant_100, pred_wcov$MDD)
        pred_p_quant_100_table_top_odds = pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "100",][2]/ pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "100",][1]
        pred_p_quant_100_table_bottom_odds = pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "1",][2]/ pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "1",][1]
        pred_p_quant_100_table_nottop_odds = colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "100",])[2]/colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "100",])[1]
        pred_p_quant_100_table_notbottom_odds = colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "1",])[2]/colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "1",])[1]

        pred_p_quant_20_table = table(pred_p_quant_20, pred_wcov$MDD)
        pred_p_quant_20_table_top_odds = pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "20",][2]/ pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "20",][1]
        pred_p_quant_20_table_bottom_odds = pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "1",][2]/ pred_p_quant_20_table[rownames(pred_p_quant_20_table) == "1",][1]
        pred_p_quant_20_table_nottop_odds = colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "20",])[2]/colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "20",])[1]
        pred_p_quant_20_table_notbottom_odds = colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "1",])[2]/colSums(pred_p_quant_20_table[rownames(pred_p_quant_20_table) != "1",])[1]

        pred_p_quant_50_table = table(pred_p_quant_50, pred_wcov$MDD)
        pred_p_quant_50_table_top_odds = pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "50",][2]/ pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "50",][1]
        pred_p_quant_50_table_bottom_odds = pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "1",][2]/ pred_p_quant_50_table[rownames(pred_p_quant_50_table) == "1",][1]
        pred_p_quant_50_table_nottop_odds = colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "50",])[2]/colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "50",])[1]
        pred_p_quant_50_table_notbottom_odds = colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "1",])[2]/colSums(pred_p_quant_50_table[rownames(pred_p_quant_50_table) != "1",])[1]


        decile_top_bottom_or = pred_p_quant_10_table_top_odds/pred_p_quant_10_table_bottom_odds
        decile_top_bottom_prev = pred_p_quant_10_table_top_prev/pred_p_quant_10_table_bottom_prev

        decile_top_mid_or = pred_p_quant_10_table_top_odds/pred_p_quant_10_table_mid_odds
        decile_top_mid_prev = pred_p_quant_10_table_top_prev/pred_p_quant_10_table_mid_prev

        top10_rest = pred_p_quant_10_table_top_odds/pred_p_quant_10_table_nottop_odds
        top5_rest = pred_p_quant_20_table_top_odds/pred_p_quant_20_table_nottop_odds
        top2_rest = pred_p_quant_50_table_top_odds/pred_p_quant_50_table_nottop_odds
        top1_rest = pred_p_quant_100_table_top_odds/pred_p_quant_100_table_nottop_odds

        quintile_top_bottom_or = pred_p_quant_5_table_top_odds/pred_p_quant_5_table_bottom_odds
        quintile_top_bottom_prev = pred_p_quant_5_table_top_prev/pred_p_quant_5_table_bottom_prev

        quintile_top_mid_or = pred_p_quant_5_table_top_odds/pred_p_quant_5_table_mid_odds
        quintile_top_mid_prev = pred_p_quant_5_table_top_prev/pred_p_quant_5_table_mid_prev

        if(AUC<0.5){
          AUC = 1 - AUC

          decile_top_bottom_or = pred_p_quant_10_table_bottom_odds/pred_p_quant_10_table_top_odds
          decile_top_bottom_prev = pred_p_quant_10_table_bottom_prev/pred_p_quant_10_table_top_prev

          decile_top_mid_or = pred_p_quant_10_table_bottom_odds/pred_p_quant_10_table_mid_odds
          decile_top_mid_prev = pred_p_quant_10_table_bottom_prev/pred_p_quant_10_table_mid_prev

          quintile_top_bottom_or = pred_p_quant_5_table_bottom_odds/pred_p_quant_5_table_top_odds
          quintile_top_bottom_prev = pred_p_quant_5_table_bottom_prev/pred_p_quant_5_table_top_prev

          quintile_top_mid_or = pred_p_quant_5_table_bottom_odds/pred_p_quant_5_table_mid_odds
          quintile_top_mid_prev = pred_p_quant_5_table_bottom_prev/pred_p_quant_5_table_mid_prev


          top10_rest = pred_p_quant_10_table_bottom_odds/pred_p_quant_10_table_notbottom_odds
          top5_rest = pred_p_quant_20_table_bottom_odds/pred_p_quant_20_table_notbottom_odds
          top2_rest = pred_p_quant_50_table_bottom_odds/pred_p_quant_50_table_notbottom_odds
          top1_rest = pred_p_quant_100_table_bottom_odds/pred_p_quant_100_table_notbottom_odds


        }

      }

      df_out = data.frame(split = i, model = model_type, nag_r2 = r2, test_auc = AUC,
                          decile_top_bottom_or =  decile_top_bottom_or,
                          decile_top_bottom_prev = decile_top_bottom_prev,

                          decile_top_mid_or =  decile_top_mid_or  ,
                          decile_top_mid_prev =  decile_top_mid_prev  ,

                          quintile_top_bottom_or = quintile_top_bottom_or,
                          quintile_top_bottom_prev = quintile_top_bottom_prev,

                          quintile_top_mid_or =  quintile_top_mid_or ,
                          quintile_top_mid_prev = quintile_top_mid_prev,

                          top10_rest,
                          top5_rest,
                          top2_rest,
                          top1_rest,
                          full_r2 = full_r2,
                          AUC_fullmodel_wcovar = AUC_fullmodel_wcovar,
                          beta = (coef_output)[1],
                          std_err = (coef_output)[2],
                          zscore = (coef_output)[3],
                          pval = (coef_output)[4],
                          OR=exp(as.numeric((coef_output)[1])),
                          lowerCI=exp(as.numeric((coef_output)[1])-1.96*as.numeric((coef_output)[2])),
                          upperCI=exp(as.numeric((coef_output)[1])+1.96*as.numeric((coef_output)[2]))
      )

      df = rbind(df, df_out)

    }

  }

  res = rbind(res, df)

}
fwrite(res, paste0(dir,"/results/combined_results_predwcovarperf_alttraitcoef_071025.csv"))
