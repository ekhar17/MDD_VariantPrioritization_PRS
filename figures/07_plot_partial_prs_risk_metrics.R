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


results = fread(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/nagr2_diff_eur_res.csv"))


results$class = gsub("_.*", "", results$type)
results$prs_type = gsub(".*Coef_", "", results$type)
clean_trait = function(df, col = "column",column_input){
  df$model = df[,column_input]

  df[grepl("MDD", df$model), col] = "MDD"
  df[grepl("BIP", df$model), col] = "BIP"
  df[grepl("SCZ", df$model), col] = "SCZ"
  df[grepl("ADHD", df$model), col] = "ADHD"
  df[grepl("ANX", df$model), col] = "ANX"
  return(df)
}
clean_type = function(df, col = "column", column_input){
  df$model = df[,column_input]
  df[grepl("hapmap", df$model), col] = "hapmap"
  df[grepl("logbf4", df$model), col] = "logbf4"
  df[grepl("logbf_over4", df$model), col] = "logbf4"
  df[grepl("cr95set", df$model), col] = "cr95set"
  df[grepl("twassig", df$model), col] = "twassig"
  df[grepl("TWAS", df$model), col] = "twassig"
  return(df)
}
results = clean_type(results %>% as.data.frame(), col = "group_type", column_input = "type")
results = clean_trait(results %>% as.data.frame(), col = "coef", column_input = "type")
results = clean_trait(results %>% as.data.frame(), col = "set", column_input = "type")
results = clean_trait(results %>% as.data.frame(), col = "trait", column_input = "type")
results[results $class == "AltTraitSet-MDDCoef","coef"] = "MDD"
results[results $class == "MDDSet-AltTraitCoef","set"] = "MDD"

performance = results
performance$type_type = "Trait SNP Set: Trait Coef"
performance$type_type[performance$class == "AltTraitSet-MDDCoef"] = "AltTrait SNP Set : MDD Coef"
performance$type_type[performance$class == "MDDSet-AltTraitCoef"] = "MDD SNP Set : AltTrait Coef"

performance$coef = factor(performance$coef, levels = c("MDD", "ADHD", "ANX", "BIP", "SCZ"))
performance_stats = performance %>% group_by(type,class,trait,prs_type, set, coef, group_type,type_type) %>% summarise(
  mean_auc = mean(nag_r2), max_auc = max(nag_r2)
)

performance$class = factor(performance$class,
                           levels = c("MDDSet-MDDCoef","AltTraitSet-AltTraitCoef",
                                      "AltTraitSet-MDDCoef","MDDSet-AltTraitCoef"   ))
library(ggpattern)

ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/TraitSNPSet_TraitCoef_AllEUR_barplot_wpattern_r2.png"),
       performance %>% filter( class == "MDDSet-MDDCoef") %>%
         ggplot( aes(y = nag_r2*100, fill = group_type, x = group_type))+geom_bar(stat = "identity",position = "dodge") +
         scale_pattern_manual(values = c( "none","stripe"), labels = c( "MDD","Alt Trait"))+
         scale_fill_manual(values = c("pink1", "gold2", "dodgerblue3", "darkgreen"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats[
             (performance_stats$class %in% c("MDDSet-MDDCoef") )
             , ],
           aes(y = max_auc*100, label = round(mean_auc*100, 3)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y =bquote('Nagelkerke' ~~ R^'2'), x = "Set", fill = "Set", pattern = "Class")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+
         guides(
           fill = "none",
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black"))
         ),
       height = 6, width = 8)


df_keep_wcovar = fread( paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_1kg/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_QCset_prediction_wMDDstatus_wcovar.txt")))
df_keep_wcovar$type = "MDDSet-MDDCoef"
LD_ref = "1kg"
colnames(df_keep_wcovar)[3:6] = paste("MDDSet-MDDCoef",colnames(df_keep_wcovar), sep = "_")[3:6]

pred_wcov = df_keep_wcovar %>% as.data.frame()

df_out_full = data.frame()
for(k in 3:6){

  pred_p_quant_10 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.10))),
                            labels=c("0-10", "10-20","20-30",
                                     "30-40","40-50", "50-60",
                                     "60-70","70-80","80-90",
                                     "90-100"), include.lowest=TRUE))
  pred_p_quant_100 = try(cut(pred_wcov[,k], breaks=c(quantile(pred_wcov[,k], probs = seq(0, 1, by = 0.01))),
                             labels=c(1:100), include.lowest=TRUE))

  pred_p_quant_10_table = table(pred_p_quant_10, pred_wcov$MDD)

  pred_p_quant_10_table_top_odds = pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][2]/pred_p_quant_10_table[rownames(pred_p_quant_10_table) == "90-100",][1]
  pred_p_quant_10_table_nottop_odds = colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "90-100",])[2]/colSums(pred_p_quant_10_table[rownames(pred_p_quant_10_table) != "90-100",])[1]
  top10_rest = pred_p_quant_10_table_top_odds/pred_p_quant_10_table_nottop_odds


  pred_p_quant_100_table = table(pred_p_quant_100, pred_wcov$MDD)
  pred_p_quant_100_table_top_odds = pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "100",][2]/ pred_p_quant_100_table[rownames(pred_p_quant_100_table) == "100",][1]
  pred_p_quant_100_table_nottop_odds = colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "100",])[2]/colSums(pred_p_quant_100_table[rownames(pred_p_quant_100_table) != "100",])[1]
  top1_rest = pred_p_quant_100_table_top_odds/pred_p_quant_100_table_nottop_odds
  df_out = data.frame(model = colnames(pred_wcov)[k],  top10_rest, top1_rest)
  df_out_full = rbind(df_out_full, df_out)
}

df_out_full$group_type = sub("_prs","", df_out_full$model)
df_out_full$group_type[df_out_full$group_type == "logbf_over4"] = "logbf4"


performance_stats = df_out_full %>% group_by( group_type) %>% summarise(
  mean_top10 = mean(top10_rest), max_top10 = max(top10_rest),
  mean_top1 = mean(top1_rest), max_top1 = max(top1_rest)
)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/TraitSNPSet_TraitCoef_AllEUR_barplot_top10.png"),
       df_out_full %>%
         ggplot( aes(y = top10_rest, fill = group_type, x = group_type))+geom_bar(stat = "identity",position = "dodge") +
         scale_pattern_manual(values = c( "none","stripe"), labels = c( "MDD","Alt Trait"))+
         scale_fill_manual(values = c("pink1", "gold2", "dodgerblue3", "darkgreen"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats,
           aes(y =mean_top10, label = round(mean_top10, 3)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y ="Top 10% OR", x = "Set", fill = "Set", pattern = "Class")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.15)))+
         guides(
           fill = "none",
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black"))
         ),
       height = 6, width = 8)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/TraitSNPSet_TraitCoef_AllEUR_barplot_top1.png"),
       df_out_full %>%
         ggplot( aes(y = top1_rest, fill = group_type, x = group_type))+geom_bar(stat = "identity",position = "dodge") +
         scale_pattern_manual(values = c( "none","stripe"), labels = c( "MDD","Alt Trait"))+
         scale_fill_manual(values = c("pink1", "gold2", "dodgerblue3", "darkgreen"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats,
           aes(y =mean_top1, label = round(mean_top1, 3)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y ="Top 1% OR", x = "Set", fill = "Set", pattern = "Class")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.15)))+
         guides(
           fill = "none",
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black"))
         ),
       height = 6, width = 8)
