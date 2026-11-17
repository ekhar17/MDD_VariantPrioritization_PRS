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

library(ggpattern)
library(colorspace)

library(data.table)
library(tidyverse)
df_SAS = fread( mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/performance_SAS_071025.csv"))
df_AFR = fread( mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/performance_AFR_071025.csv"))

ggplot(df_SAS, aes(y = test_auc, x = model))+geom_point()

df_SAS$split = NA
df_AFR$split = NA
df_SAS$type = "NoCombo"
df_AFR$type = "NoCombo"


dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/single_trait/")
df_SAS_singlcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_SAS_071025.csv"))
df_AFR_singlcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_AFR_071025.csv"))

df_SAS_singlcombo$ancestry = "SAS"
df_AFR_singlcombo$ancestry = "AFR"
df_SAS_singlcombo$type = "Single"
df_AFR_singlcombo$type = "Single"

dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/multi_trait/")
df_SAS_multcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_SAS_multitrait_071025.csv"))
df_AFR_multcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_AFR_multitrait_071025.csv"))
df_SAS_multcombo$ancestry = "SAS"
df_AFR_multcombo$ancestry = "AFR"
df_SAS_multcombo$type = "Multi"
df_AFR_multcombo$type = "Multi"

dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/muti_trait_MDDcoef_alttraitset/")
df_SAS_altraitsetcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_SAS_alttraitset_071025.csv"))
df_AFR_altraitsetcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_AFR_alttraitset_071025.csv"))
df_SAS_altraitsetcombo$ancestry = "SAS"
df_AFR_altraitsetcombo$ancestry = "AFR"
df_SAS_altraitsetcombo$type = "AltTraitSet"
df_AFR_altraitsetcombo$type = "AltTraitSet"


dir = mdd_path("MDD_PROJECT", "/06.Elastic_Net/muti_trait_MDDset_alttraitcoef/")
df_SAS_altraitcoefcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_SAS_alttraitcoef_071025.csv"))
df_AFR_altraitcoefcombo = fread( paste0(dir,"/results_altancestry/combined_results_covar0_AFR_alttraitcoef_071025.csv"))
df_SAS_altraitcoefcombo$ancestry = "SAS"
df_AFR_altraitcoefcombo$ancestry = "AFR"
df_SAS_altraitcoefcombo$type = "AltTraitCoef"
df_AFR_altraitcoefcombo$type = "AltTraitCoef"

df_SAS_all = rbind(df_SAS, df_SAS_singlcombo, df_SAS_multcombo,
                   #df_SAS_altraitsetcombo,
                   df_SAS_altraitcoefcombo)
df_AFR_all = rbind(df_AFR, df_AFR_singlcombo, df_AFR_multcombo,
                   #df_AFR_altraitsetcombo,
                   df_AFR_altraitcoefcombo)
df_SAS_all_long = df_SAS_all %>% pivot_longer(colnames(df_SAS_all)[2:25])
df_AFR_all_long = df_AFR_all %>% pivot_longer(colnames(df_AFR_all)[2:25])

df_SAS_all_long %>% group_by(name) %>% summarise(best_mod = model[which.max(value)],
                                                 best_type = type[which.max(value)])
df_AFR_all_long %>% group_by(name) %>% summarise(best_mod = model[which.max(value)],
                                                 best_type = type[which.max(value)])

df_SAS = df_SAS[grepl("MDD",df_SAS$model),]
df_AFR = df_AFR[grepl("MDD",df_AFR$model),]
df_SAS = df_SAS[!grepl("twas",df_SAS$model),]
df_AFR = df_AFR[!grepl("twas",df_AFR$model),]


df_SAS_singlcombo = df_SAS_singlcombo[df_SAS_singlcombo$model %in% c("ALLCovar_NoTWAS_caret_roc"),]
df_AFR_singlcombo = df_AFR_singlcombo[df_AFR_singlcombo$model %in% c("ALLCovar_NoTWAS_caret_roc"),]

df_SAS_multcombo = df_SAS_multcombo[df_SAS_multcombo$model %in% c("cr95setlogbf4_ALLCovar_caret_roc", "cr95setonly_ALLCovar_caret_roc"),]
df_AFR_multcombo = df_AFR_multcombo[df_AFR_multcombo$model %in% c("cr95setlogbf4_ALLCovar_caret_roc", "cr95setonly_ALLCovar_caret_roc"),]

df_SAS_altraitsetcombo = df_SAS_altraitsetcombo[df_SAS_altraitsetcombo$model %in% c("cr95setlogbf4_ALLCovar_caret_roc", "cr95setonly_ALLCovar_caret_roc"),]
df_AFR_altraitsetcombo = df_AFR_altraitsetcombo[df_AFR_altraitsetcombo$model %in% c("cr95setlogbf4_ALLCovar_caret_roc", "cr95setonly_ALLCovar_caret_roc"),]

df_SAS_altraitcoefcombo = df_SAS_altraitcoefcombo[df_SAS_altraitcoefcombo$model %in% c("cr95setlogbf4_ALLCovar_caret_roc", "cr95setonly_ALLCovar_caret_roc"
),]

df_AFR_altraitcoefcombo = df_AFR_altraitcoefcombo[df_AFR_altraitcoefcombo$model %in% c("cr95setlogbf4_ALLCovar_caret_roc", "cr95setonly_ALLCovar_caret_roc"),]
df_SAS_all = rbind(df_SAS, df_SAS_singlcombo, df_SAS_multcombo,
                   df_SAS_altraitsetcombo,
                   df_SAS_altraitcoefcombo)
df_AFR_all = rbind(df_AFR, df_AFR_singlcombo, df_AFR_multcombo,
                   df_AFR_altraitsetcombo,
                   df_AFR_altraitcoefcombo)
df_SAS_all$OR[df_SAS_all$OR<1] = 1/df_SAS_all$OR[df_SAS_all$OR<1]
df_AFR_all$OR[df_AFR_all$OR<1] = 1/df_AFR_all$OR[df_AFR_all$OR<1]
df_SAS_all_long = df_SAS_all %>% pivot_longer(colnames(df_SAS_all)[2:25])
df_AFR_all_long = df_AFR_all %>% pivot_longer(colnames(df_AFR_all)[2:25])

df_SAS_all_long %>% group_by(name) %>% summarise(best_mod = model[which.max(value)],
                                                 best_type = type[which.max(value)])
df_AFR_all_long %>% group_by(name) %>% summarise(best_mod = model[which.max(value)],
                                                 best_type = type[which.max(value)])

df_SAS_all_long = df_SAS_all_long %>%
  filter(name %in% c("cor_pheno_adj", "nag_r2", "test_auc",
                     "top1_rest", "top10_rest","top5_rest", "top2_rest",
                     "OR",  "AUC_fullmodel_wcovar"    ))

df_AFR_all_long = df_AFR_all_long %>%
  filter(name %in% c("cor_pheno_adj", "nag_r2", "test_auc",
                     "top1_rest", "top10_rest","top5_rest", "top2_rest", "OR",  "AUC_fullmodel_wcovar"))
df_SAS_all_long = df_SAS_all_long %>%
  filter(name %in% c("cor_pheno_adj", "nag_r2", "top10_rest","top1_rest" ))

df_AFR_all_long = df_AFR_all_long %>%
  filter(name %in% c("cor_pheno_adj", "nag_r2",
                     "top1_rest", "top10_rest"))


df_SAS_all_long$model[df_SAS_all_long$model == "ALLCovar_NoTWAS_caret_roc"] = "cr95setlogbf4_ALLCovar_caret_roc"

df_SAS_all_long$type[df_SAS_all_long$type == "NoCombo"] = "MDD"
df_SAS_all_long$type[df_SAS_all_long$type == "Single"] = "MDD Combo"
df_SAS_all_long$type[df_SAS_all_long$type == "AltTraitCoef"] = "MDDSet-AltTraitCoef"
df_SAS_all_long$type[df_SAS_all_long$type == "AltTraitSet"] = "AltTraitSet-MDDCoef"
df_SAS_all_long$type[df_SAS_all_long$type == "Multi"] = "Multi Trait"
df_SAS_all_long$type[df_SAS_all_long$type == "MDD Combo"] = "MDD"
df_AFR_all_long$type[df_AFR_all_long$type == "NoCombo"] = "MDD"
df_AFR_all_long$type[df_AFR_all_long$type == "Single"] = "MDD Combo"
df_AFR_all_long$type[df_AFR_all_long$type == "AltTraitCoef"] = "MDDSet-AltTraitCoef"
df_AFR_all_long$type[df_AFR_all_long$type == "Multi"] = "Multi Trait"
df_AFR_all_long$type[df_AFR_all_long$type == "AltTraitSet"] = "AltTraitSet-MDDCoef"
df_AFR_all_long$type[df_AFR_all_long$type == "MDD Combo"] = "MDD"

df_AFR_all_long$type = factor(df_AFR_all_long$type, levels = c("MDD","Multi Trait",
                                                               "MDD Combo", "AltTraitSet-MDDCoef" ,
                                                               "MDDSet-AltTraitCoef"
                                                              ))
df_SAS_all_long$type = factor(df_SAS_all_long$type, levels = c("MDD","Multi Trait",
                                                               "MDD Combo",  "AltTraitSet-MDDCoef" ,
                                                               "MDDSet-AltTraitCoef" ))


df_SAS_all_long$class = "cr95set"
df_SAS_all_long$class[grepl("hapmap", df_SAS_all_long$model)] = "hapmap"
df_SAS_all_long$class[grepl("logbf4", df_SAS_all_long$model)] = "cr95set+logbf4"
df_SAS_all_long$class[grepl("logbf_over4", df_SAS_all_long$model)] = "logbf4"


df_AFR_all_long$class = "cr95set"
df_AFR_all_long$class[grepl("hapmap", df_AFR_all_long$model)] = "hapmap"
df_AFR_all_long$class[grepl("logbf4", df_AFR_all_long$model)] = "cr95set+logbf4"
df_AFR_all_long$class[grepl("logbf_over4", df_AFR_all_long$model)] = "logbf4"


df_SAS_all_long$name[df_SAS_all_long$name == "cor_pheno_adj"] = "Correlation"
df_SAS_all_long$name[df_SAS_all_long$name == "nag_r2"] = "Nag R2"
df_SAS_all_long$name[df_SAS_all_long$name == "top10_rest"] = "Top 10% OR"
df_SAS_all_long$name[df_SAS_all_long$name == "top1_rest"] = "Top 1% OR"
df_SAS_all_long$name = factor(df_SAS_all_long$name, levels = c("Correlation", "Nag R2", "Top 10% OR", "Top 1% OR"))


df_SAS_all_long_summary = df_SAS_all_long %>% group_by(model, name, type, class) %>% summarise(mean = mean(value), max = max(value))

df_SAS_all_long_summary %>% group_by(name) %>% summarise(model[which.max(mean)], type[which.max(mean)])


df_SAS_all_long_summary_MDD = df_SAS_all_long_summary[df_SAS_all_long_summary$type == "MDD" & df_SAS_all_long_summary$class != "cr95set+logbf4",]

df_SAS_all_long_summary_MDD$mean2 = df_SAS_all_long_summary_MDD$mean
df_SAS_all_long_summary_MDD$mean2[df_SAS_all_long_summary_MDD$name == "Nag R2" & df_SAS_all_long_summary_MDD$class == "cr95set"] =

  df_SAS_all_long_summary_MDD$mean2[df_SAS_all_long_summary_MDD$name == "Nag R2" & df_SAS_all_long_summary_MDD$class == "cr95set"] + 0.0003


df_SAS_all_long_summary_MDD$mean2[df_SAS_all_long_summary_MDD$name == "Top 1% OR" & df_SAS_all_long_summary_MDD$class == "cr95set"] =

  df_SAS_all_long_summary_MDD$mean2[df_SAS_all_long_summary_MDD$name == "Top 1% OR" & df_SAS_all_long_summary_MDD$class == "cr95set"] + 0.2

df_SAS_all_long_summary_MDD$mean2[df_SAS_all_long_summary_MDD$name == "Top 1% OR" & df_SAS_all_long_summary_MDD$class == "hapmap"] =

  df_SAS_all_long_summary_MDD$mean2[df_SAS_all_long_summary_MDD$name == "Top 1% OR" & df_SAS_all_long_summary_MDD$class == "hapmap"] + 0.4

df_SAS_all_long_summary$mean[df_SAS_all_long_summary$mean>1] = round(df_SAS_all_long_summary$mean[df_SAS_all_long_summary$mean>1],2)
df_SAS_all_long_summary$mean[df_SAS_all_long_summary$mean<1] = round(df_SAS_all_long_summary$mean[df_SAS_all_long_summary$mean<1],3)
df_SAS_all_long_summary$max[df_SAS_all_long_summary$mean==1.33] = df_SAS_all_long_summary$max[df_SAS_all_long_summary$mean==1.33] + 0.01
df_SAS_all_long_summary$max[df_SAS_all_long_summary$mean==1.32] = df_SAS_all_long_summary$max[df_SAS_all_long_summary$mean==1.32] - 0.01

ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/SAS_test_set_comboplot_cor_nagr2_top10_top1.png"),


       df_SAS_all_long %>%
         ggplot(aes(x = type, fill =class, y = value)) +
         facet_wrap(~name, scales = "free_y") +
         geom_boxplot_pattern( data = df_SAS_all_long[
           !(df_SAS_all_long$type == "MDD" & df_SAS_all_long$class == "hapmap" )&
             !(df_SAS_all_long$type == "MDD" & df_SAS_all_long$class == "logbf4" ),],
           aes(pattern = type))+
         scale_fill_manual(values = c("pink1", lighten("purple1", amount = 0.2), "red", "snow1", "white", "#2A9D8F","#AEC7D3")) +
         scale_pattern_manual(values = c("none","stripe", "circle", "crosshatch"))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +


         geom_hline(data = df_SAS_all_long_summary_MDD,

                    aes(color = class, yintercept = mean), lty = 2)+
         geom_text(data = df_SAS_all_long_summary_MDD,
                   aes(0, mean2, label = sub("logbf_over4", "logbf4", sub("_"," ",sub("_prs","",model))), vjust = -1))+theme_bw() +
         labs(title = "SAS Ancestry", fill = "Set Type", pattern = "Class", x = "PRS")+
         scale_color_manual(values = c("maroon3", "gold4", "darkblue"), guide = "none")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.2)))+ scale_x_discrete(labels = c("MDD", "Multi Trait",
                                                                                                                      "AltTraitSet-MDDCoef", "MDDSet-AltTraitCoef"),
                                                                                                           expand = expansion(mult = c(0.55, 0.18)))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               strip.text.x = element_text(size = 14, colour = "black"),
               text = element_text(size = 14, color = "black"),
               legend.key.size = unit(1, 'cm'))       + guides(
                 fill = guide_legend(override.aes = list(pattern = "none")),
                 pattern = guide_legend(override.aes = list(fill = "white",
                                                            color = "black")))+
         expand_limits(x=-0.5)+
         geom_text( data = df_SAS_all_long_summary[!df_SAS_all_long_summary$class %in% c("logbf4", "hapmap") &
                                                     !df_SAS_all_long_summary$name %in% c("Correlation","Nag R2"),],
                    aes( y = max, label = mean, x = type),
                    position = position_dodge(width = 0.75), # aligns with boxplot dodging
                    vjust = -1
                    #   size = 6
         )+
         geom_text( data = df_SAS_all_long_summary[!df_SAS_all_long_summary$class %in% c("logbf4", "hapmap") &
                                                     df_SAS_all_long_summary$name %in% c("Correlation","Nag R2"),],
                    aes( y = max, label = mean, x = type),
                    position = position_dodge(width = .9), # aligns with boxplot dodging
                    vjust = -1
                    #   size = 6
         )


       ,

       height = 8, width = 14.5)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/SAS_test_set_comboplot_cor_nagr2_top10_top1_small.png"),


       df_SAS_all_long %>%
         ggplot(aes(x = type, fill =class, y = value)) +
         facet_wrap(~name, scales = "free_y") +
         geom_boxplot_pattern( data = df_SAS_all_long[
           !(df_SAS_all_long$type == "MDD" & df_SAS_all_long$class == "hapmap" )&
             !(df_SAS_all_long$type == "MDD" & df_SAS_all_long$class == "logbf4" ),],
           aes(pattern = type))+
         scale_fill_manual(values = c("pink1", lighten("purple1", amount = 0.2), "red", "snow1", "white", "#2A9D8F","#AEC7D3")) +
         scale_pattern_manual(values = c("none","stripe", "circle", "crosshatch"))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +


         geom_hline(data = df_SAS_all_long_summary_MDD,

                    aes(color = class, yintercept = mean), lty = 2)+
         geom_text(data = df_SAS_all_long_summary_MDD,
                   aes(0, mean2, label = sub("logbf_over4", "logbf4", sub("_"," ",sub("_prs","",model))), vjust = -1))+theme_bw() +
         labs(title = "SAS", fill = "Set Type", pattern = "Class", x = "PRS")+
         scale_color_manual(values = c("maroon3", "gold4", "darkblue"), guide = "none")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.2)))+scale_x_discrete(labels = c("MDD", "Multi Trait",
                                                                                                                     "AltTraitSet-MDDCoef", "MDDSet-AltTraitCoef"),
                                                                                                          expand = expansion(mult = c(0.55, 0.18)))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               strip.text.x = element_text(size = 14, colour = "black"),
               text = element_text(size = 14, color = "black"),
               legend.key.size = unit(1, 'cm'))       + guides(
                 fill = guide_legend(override.aes = list(pattern = "none")),
                 pattern = guide_legend(override.aes = list(fill = "white",
                                                            color = "black")))+
         expand_limits(x=-0.5)
       ,

       height = 8, width = 8)


df_AFR_all_long$name[df_AFR_all_long$name == "cor_pheno_adj"] = "Correlation"
df_AFR_all_long$name[df_AFR_all_long$name == "nag_r2"] = "Nag R2"
df_AFR_all_long$name[df_AFR_all_long$name == "top10_rest"] = "Top 10% OR"
df_AFR_all_long$name[df_AFR_all_long$name == "top1_rest"] = "Top 1% OR"
df_AFR_all_long$name = factor(df_AFR_all_long$name, levels = c("Correlation", "Nag R2", "Top 10% OR", "Top 1% OR"))


df_AFR_all_long$class[df_AFR_all_long$model == "ALLCovar_NoTWAS_caret_roc"] = "cr95set+logbf4"

df_AFR_all_long[df_AFR_all_long$model == "MDD",]


df_AFR_all_long_summary = df_AFR_all_long %>% group_by(model, name, type, class) %>% summarise(mean = mean(value), max = max(value))

df_AFR_all_long_summary %>% group_by(name) %>% summarise(model[which.max(mean)], type[which.max(mean)])


df_AFR_all_long_summary_MDD = df_AFR_all_long_summary[df_AFR_all_long_summary$type == "MDD" & df_AFR_all_long_summary$model != "ALLCovar_NoTWAS_caret_roc",]


df_AFR_all_long_summary_MDD$mean2 = df_AFR_all_long_summary_MDD$mean
df_AFR_all_long_summary_MDD$mean2[df_AFR_all_long_summary_MDD$name == "Top 10% OR" & df_AFR_all_long_summary_MDD$class == "hapmap"] =

  df_AFR_all_long_summary_MDD$mean2[df_AFR_all_long_summary_MDD$name == "Top 10% OR" & df_AFR_all_long_summary_MDD$class == "hapmap"] + 0.03


df_AFR_all_long_summary_MDD$mean2[df_AFR_all_long_summary_MDD$name == "Top 1% OR" & df_AFR_all_long_summary_MDD$class == "cr95set"] =

  df_AFR_all_long_summary_MDD$mean2[df_AFR_all_long_summary_MDD$name == "Top 1% OR" & df_AFR_all_long_summary_MDD$class == "cr95set"] + 0.1

df_AFR_all_long_summary$mean[df_AFR_all_long_summary$mean>1] = round(df_AFR_all_long_summary$mean[df_AFR_all_long_summary$mean>1],2)
df_AFR_all_long_summary$mean[df_AFR_all_long_summary$mean<1] = round(df_AFR_all_long_summary$mean[df_AFR_all_long_summary$mean<1],3)
df_AFR_all_long_summary$max[df_AFR_all_long_summary$mean==1.43] = 1.4
df_AFR_all_long_summary$max[df_AFR_all_long_summary$mean==0.874] = 0.868
df_AFR_all_long_summary$max[df_AFR_all_long_summary$mean==1.1 &
                              df_AFR_all_long_summary$name == "Top 10% OR"] = df_AFR_all_long_summary$max[df_AFR_all_long_summary$mean==1.1 &
                                                                                                            df_AFR_all_long_summary$name == "Top 10% OR"] - 0.005
df_AFR_all_long_summary$max[df_AFR_all_long_summary$model == "MDD_cr95set_prs" &
                              df_AFR_all_long_summary$name == "Nag R2"] = df_AFR_all_long_summary$max[df_AFR_all_long_summary$model == "MDD_cr95set_prs" &
                                                                                                        df_AFR_all_long_summary$name == "Nag R2"] + 0.00009

df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set+logbf4" &
                              df_AFR_all_long_summary$type == "Multi Trait" &
                              df_AFR_all_long_summary$name == "Nag R2"] = df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set+logbf4" &
                                                                                                        df_AFR_all_long_summary$type == "Multi Trait" &
                                                                                                        df_AFR_all_long_summary$name == "Nag R2"] + 0.00009

df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set" &
                              df_AFR_all_long_summary$type == "MDD SNP Set: AltTrait Coef" &
                              df_AFR_all_long_summary$name == "Nag R2"] = df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set" &
                                                                                                        df_AFR_all_long_summary$type == "MDD SNP Set: AltTrait Coef" &
                                                                                                        df_AFR_all_long_summary$name == "Nag R2"] + 0.00011


df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set" &
                              df_AFR_all_long_summary$type == "AltTrait SNP Set: MDD Coef" &
                              df_AFR_all_long_summary$name == "Nag R2"] = df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set" &
                                                                                                        df_AFR_all_long_summary$type == "AltTrait SNP Set: MDD Coef" &
                                                                                                        df_AFR_all_long_summary$name == "Nag R2"] - 0.00004
df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set+logbf4" &
                              df_AFR_all_long_summary$type == "AltTrait SNP Set: MDD Coef" &
                              df_AFR_all_long_summary$name == "Nag R2"] = df_AFR_all_long_summary$max[df_AFR_all_long_summary$class == "cr95set+logbf4" &
                                                                                                        df_AFR_all_long_summary$type == "AltTrait SNP Set: MDD Coef" &
                                                                                                        df_AFR_all_long_summary$name == "Nag R2"] - 0.00010

ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/AFR_test_set_comboplot_cor_nagr2_top10_top1.png"),


       df_AFR_all_long %>%
         ggplot(aes(x = type, fill =class, y = value)) +
         facet_wrap(~name, scales = "free_y") +
         geom_boxplot_pattern( data = df_AFR_all_long[
           !(df_AFR_all_long$type == "MDD" & df_AFR_all_long$class == "hapmap" )&
             !(df_AFR_all_long$type == "MDD" & df_AFR_all_long$class == "logbf4" ),],
           aes(pattern = type))+
         scale_fill_manual(values = c("pink1", lighten("purple1", amount = 0.2), "red", "snow1", "white", "#2A9D8F","#AEC7D3")) +
         scale_pattern_manual(values = c("none","stripe", "circle", "crosshatch"))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +


         geom_hline(data = df_AFR_all_long_summary_MDD,

                    aes(color = class, yintercept = mean), lty = 2)+
         geom_text(data = df_AFR_all_long_summary_MDD,
                   aes(0, mean2, label = sub("logbf_over4", "logbf4", sub("_"," ",sub("_prs","",model))), vjust = -1))+theme_bw() +
         labs(title = "AFR Ancestry", fill = "Set Type", pattern = "Class", x = "PRS")+
         scale_color_manual(values = c("maroon3", "gold4", "darkblue"), guide = "none")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.2)))+ scale_x_discrete(labels = c("MDD", "Multi Trait",
                                                                                                                      "AltTraitSet-MDDCoef", "MDDSet-AltTraitCoef"),
                                                                                                           expand = expansion(mult = c(0.55, 0.18)))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               strip.text.x = element_text(size = 14, colour = "black"),
               text = element_text(size = 14, color = "black"),
               legend.key.size = unit(1, 'cm'))       + guides(
                 fill = guide_legend(override.aes = list(pattern = "none")),
                 pattern = guide_legend(override.aes = list(fill = "white",
                                                            color = "black")))+
         expand_limits(x=-0.5)+
         geom_text( data = df_AFR_all_long_summary[!df_AFR_all_long_summary$class %in% c("logbf4", "hapmap") &
                                                     !df_AFR_all_long_summary$name %in% c("Correlation","Nag R2"),],
                    aes( y = max, label = mean, x = type),
                    position = position_dodge(width = 0.75), # aligns with boxplot dodging
                    vjust = -1
                    #   size = 6
         )+
         geom_text( data = df_AFR_all_long_summary[!df_AFR_all_long_summary$class %in% c("logbf4", "hapmap") &
                                                     df_AFR_all_long_summary$name %in% c("Correlation","Nag R2"),],
                    aes( y = max, label = mean, x = type),
                    position = position_dodge(width = .9), # aligns with boxplot dodging
                    vjust = -1
                    #   size = 6
         )

       ,

       height = 8, width = 14.5)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/AFR_test_set_comboplot_cor_nagr2_top10_top1_small.png"),


       df_AFR_all_long %>%
         ggplot(aes(x = type, fill =class, y = value)) +
         facet_wrap(~name, scales = "free_y") +
         geom_boxplot_pattern( data = df_AFR_all_long[
           !(df_AFR_all_long$type == "MDD" & df_AFR_all_long$class == "hapmap" )&
             !(df_AFR_all_long$type == "MDD" & df_AFR_all_long$class == "logbf4" ),],
           aes(pattern = type))+
         scale_fill_manual(values = c("pink1", lighten("purple1", amount = 0.2), "red", "snow1", "white", "#2A9D8F","#AEC7D3")) +
         scale_pattern_manual(values = c("none","stripe", "circle", "crosshatch"))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +


         geom_hline(data = df_AFR_all_long_summary_MDD,

                    aes(color = class, yintercept = mean), lty = 2)+
         geom_text(data = df_AFR_all_long_summary_MDD,
                   aes(0, mean2, label = sub("logbf_over4", "logbf4", sub("_"," ",sub("_prs","",model))), vjust = -1))+theme_bw() +
         labs(title = "AFR Ancestry", fill = "Set Type", pattern = "Class", x= "PRS")+
         scale_color_manual(values = c("maroon3", "gold4", "darkblue"), guide = "none")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.2)))+scale_x_discrete(labels = c("MDD", "Multi Trait",
                                                                                                                     "AltTraitSet-MDDCoef", "MDDSet-AltTraitCoef"),
                                                                                                          expand = expansion(mult = c(0.55, 0.18)))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               strip.text.x = element_text(size = 14, colour = "black"),
               text = element_text(size = 14, color = "black"),
               legend.key.size = unit(1, 'cm'))+ guides(
                 fill = guide_legend(override.aes = list(pattern = "none")),
                 pattern = guide_legend(override.aes = list(fill = "white",
                                                            color = "black"))),

       height = 8, width = 8)
