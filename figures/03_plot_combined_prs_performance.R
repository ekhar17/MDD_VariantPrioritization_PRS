# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpattern)
  library(colorspace)
})

alttraitcoef_results = fread(mdd_path("MDD_PROJECT", "/06.Elastic_Net/muti_trait_MDDset_alttraitcoef//results/combined_results_covar0_alttraitcoef_071125.csv")) %>% as.data.frame()
alttraitset_results = fread(mdd_path("MDD_PROJECT", "/06.Elastic_Net/muti_trait_MDDcoef_alttraitset/results/combined_results_covar0_alttraitset_071125.csv"))%>% as.data.frame()
multitrait_results = fread(mdd_path("MDD_PROJECT", "/06.Elastic_Net/multi_trait/results/combined_results_covar0_multitrait_071125.csv"))%>% as.data.frame()
singletrait_results = fread(mdd_path("MDD_PROJECT", "/06.Elastic_Net/single_trait/results/combined_results_covar0_071125.csv"))%>% as.data.frame()

singletrait_results$type  = "MDD"
multitrait_results$type  = "Multi Trait"
alttraitset_results$type = "AltTraitSet-MDDCoef"
alttraitcoef_results$type = "MDDSet-AltTraitCoef"

df = rbind(singletrait_results %>% filter(model %in% c("cr95set",
                                                       "ALLCovar_NoTWAS_caret_roc")),
           multitrait_results %>% filter(model %in% c("cr95setlogbf4_ALLCovar_caret_roc",
                                                      "cr95setonly_ALLCovar_caret_roc")),
           alttraitset_results %>% filter(model %in% c("cr95setlogbf4_ALLCovar_caret_roc",
                                                       "cr95setonly_ALLCovar_caret_roc")),
           alttraitcoef_results %>% filter(model %in% c("cr95setlogbf4_ALLCovar_caret_roc",
                                                        "cr95setonly_ALLCovar_caret_roc")))
#df $OR[df$OR < 1] = 1/df$OR[df$OR < 1]
df_long = pivot_longer(df, cols = colnames(df)[3:30])

df_long = df_long %>%
  filter(name %in% c("cor_pheno_adj", "nag_r2", "test_auc",
                     "top1_rest", "top10_rest","top5_rest", "top2_rest",
                     "OR", "AUC_fullmodel_wcovar", "top10_rest_mod",           "top5_rest_mod"        ,    "top2_rest_mod"         ,   "top1_rest_mod"  ))
df_long = df_long %>%
  filter(name %in% c("cor_pheno_adj", "nag_r2",
                     "top10_rest","top1_rest"))
df_long$name[df_long$name == "cor_pheno_adj"] = "Correlation"
df_long$name[df_long$name == "nag_r2"] = "Nag R2"
df_long$name[df_long$name == "top10_rest"] = "Top 10% OR"
df_long$name[df_long$name == "top1_rest"] = "Top 1% OR"
df_long$name = factor(df_long$name, levels = c("Correlation", "Nag R2", "Top 10% OR", "Top 1% OR"))
df_long$type = factor(df_long$type, levels = c("MDD", "Multi Trait", "AltTraitSet-MDDCoef", "MDDSet-AltTraitCoef"))

df_long$set_type = "cr95set"
df_long$set_type[grepl("logbf4", df_long$model)] = "cr95set+logbf4"
df_long$set_type[grepl("NoTWAS", df_long$model)] = "cr95set+logbf4"

df_long_summary = df_long %>% group_by(model, name, type, set_type) %>% summarise(mean = mean(value), max = max(value))
df_long_summary$model1 = paste(df_long_summary$model, df_long_summary$type)
df_long_summary %>% group_by(name) %>% summarise(model1[which.max(mean)])
cr95set_summary = df_long_summary[df_long_summary$model1 == "cr95set MDD",]

table(df_long$type)
table(df_long$set_type)
table(df_long$name)
df_long[df_long$type == "Multi Trait"&
          df_long$set_type == "cr95set+logbf4",]
df_long[df_long$type == "MDD" &
          df_long$set_type == "cr95set+logbf4",]


t.test(df_long[df_long$type == "Multi Trait"&
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Correlation",]$value -
df_long[df_long$type == "MDD" &
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Correlation",]$value)


t.test(
df_long[df_long$type == "Multi Trait"&
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Top 1% OR",]$value -
df_long[df_long$type == "MDD" &
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Top 1% OR",]$value)


t.test(
df_long[df_long$type == "Multi Trait"&
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Top 10% OR",]$value -
df_long[df_long$type == "MDD" &
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Top 10% OR",]$value)


t.test(
df_long[df_long$type == "Multi Trait"&
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Nag R2",]$value -
df_long[df_long$type == "MDD" &
          df_long$set_type == "cr95set+logbf4" &
          df_long$name == "Nag R2",]$value)
df_long
df_long_summary %>% group_by(name) %>% summarise()
library(ggpattern)
library(colorspace)
ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/EUR_test_set_comboplot_cor_nagr2_top10_top1.png"),
       df_long %>%
         ggplot(aes(fill = set_type, y = value, x = type))+facet_wrap(~name, scales = "free_y")+geom_boxplot_pattern(aes(pattern = type)) +
         scale_fill_manual(values = c("pink1", lighten("purple1", amount = 0.2), "#2A9D8F", "snow1", "white"))  +
         scale_pattern_manual(values = c("none","stripe", "circle", "crosshatch"))+ geom_text(
           data = df_long_summary,
           aes( y = max, label = round(mean, 3), x = type, color = set_type),
           position = position_dodge(width = 0.75), # aligns with boxplot dodging
           vjust = -1
           #   size = 6
         )+
         theme_bw()+
         geom_hline(data = cr95set_summary,
                    aes(yintercept = mean), col = "red", lty = 2)+
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               strip.text.x = element_text(size = 14, colour = "black"),
               text = element_text(size = 14, color = "black"),
               legend.key.size = unit(1, 'cm')) +
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.2)))+scale_color_manual(values = c("black", "black"),
                                                                                                            guide = "none")+
         #  theme(axis.title.x=element_blank(),
         #       axis.text.x=element_blank(),
         #      axis.ticks.x=element_blank())+
         labs(pattern = "Class", fill = "Set Type", x = "PRS")+
         guides(
           fill = guide_legend(override.aes = list(pattern = "none")),
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black")))
       ,

       height = 8, width = 14.5)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/EUR_test_set_comboplot_cor_nagr2_top10_top1_small.png"),
       df_long %>%
         ggplot(aes(fill = set_type, y = value, x = type))+facet_wrap(~name, scales = "free_y")+geom_boxplot_pattern(aes(pattern = type)) +
         scale_fill_manual(values = c("pink1", lighten("purple1", amount = 0.2), "#2A9D8F", "snow1", "white"))  +
         scale_pattern_manual(values = c("none","stripe", "circle", "crosshatch"))+ geom_text(
           data = df_long_summary,
           aes( y = max, label = round(mean, 3), x = type, color = set_type),
           position = position_dodge(width = 0.75), # aligns with boxplot dodging
           vjust = -1
           #   size = 6
         )+
         theme_bw()+
         geom_hline(data = cr95set_summary,
                    aes(yintercept = mean), col = "red", lty = 2)+
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               strip.text.x = element_text(size = 14, colour = "black"),
               text = element_text(size = 14, color = "black"),
               legend.key.size = unit(1, 'cm')) +
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.2)))+scale_color_manual(values = c("black", "black"),
                                                                                                            guide = "none")+
         #  theme(axis.title.x=element_blank(),
         #       axis.text.x=element_blank(),
         #      axis.ticks.x=element_blank())+
         labs(pattern = "Class", fill = "Set Type", x = "PRS")+
         guides(
           fill = guide_legend(override.aes = list(pattern = "none")),
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black")))
       ,

       height = 8, width = 8)


df_long$model_full = paste(df_long$model, df_long$type, df_long$set_type)

table(df_long$model_full)
table(df_long$name)
df_res = data.frame()
for(cur_val in unique(df_long$name)){
  for(model_cur in unique(df_long$model_full)){
    df1 = df_long[df_long$model_full == "cr95set MDD cr95set" & df_long$name == cur_val,]
    df2 = df_long[df_long$model_full == model_cur & df_long$name == cur_val,]
    t_test = t.test(df1$value - df2$value)
    df = data.frame(measure = cur_val ,
                    model = model_cur,
                    model_t = df2$model[1],
                    type = df2$type[1],
                    set = df2$set_type[1],
                    p = t_test$p.value)
    df_res = rbind(df_res,df)
  }
}

df_res = df_res[df_res$model!= "cr95set MDD cr95set", ]
df_res1 = df_res[df_res$measure == "Nag R2",]
df_res2 = df_res[df_res$measure == "Correlation",]
df_res3 = df_res[df_res$measure == "Top 10% OR",]
df_res4 = df_res[df_res$measure == "Top 1% OR",]

df_res1[7,]
df_res2[7,]
df_res3[7,]
df_res4[7,]


table(df_long$name)
df_res11 = data.frame()
for(cur_val in unique(df_long$name)){
  for(model_cur in unique(df_long$model_full)){
    df1 = df_long[df_long$model_full == "ALLCovar_NoTWAS_caret_roc MDD cr95set+logbf4" & df_long$name == cur_val,]
    df2 = df_long[df_long$model_full == model_cur & df_long$name == cur_val,]
    t_test = t.test(df1$value - df2$value)
    df = data.frame(measure = cur_val ,
                    model = model_cur,
                    model_t = df2$model[1],
                    type = df2$type[1],
                    set = df2$set_type[1],
                    p = t_test$p.value)
    df_res11 = rbind(df_res11,df)
  }
}

df_res11 = df_res11[df_res11$model!= "ALLCovar_NoTWAS_caret_roc MDD cr95set+logbf4", ]
df_res11 = df_res11[df_res11$model!= "cr95set MDD cr95set", ]
df_res11 = df_res11[df_res11$set == "cr95set+logbf4", ]

df_res111 = df_res11[df_res11$measure == "Nag R2",]
df_res112 = df_res11[df_res11$measure == "Correlation",]
df_res113 = df_res11[df_res11$measure == "Top 10% OR",]
df_res114 = df_res11[df_res11$measure == "Top 1% OR",]


df_res1x = rbind(df_res1, df_res111)
df_res2x = rbind(df_res2, df_res112)
df_res3x = rbind(df_res3, df_res113)
df_res4x = rbind(df_res4, df_res114)

df_res1$p_adj = p.adjust(df_res1$p, method = "fdr")
df_res2$p_adj = p.adjust(df_res2$p, method = "fdr")
df_res3$p_adj = p.adjust(df_res3$p, method = "fdr")
df_res4$p_adj = p.adjust(df_res4$p, method = "fdr")

df_res1x$p_adj = p.adjust(df_res1x$p, method = "fdr")
df_res2x$p_adj = p.adjust(df_res2x$p, method = "fdr")
df_res3x$p_adj = p.adjust(df_res3x$p, method = "fdr")
df_res4x$p_adj = p.adjust(df_res4x$p, method = "fdr")

## cor
## Top 10
df_res
