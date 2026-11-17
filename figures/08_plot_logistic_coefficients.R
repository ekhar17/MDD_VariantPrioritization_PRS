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

library(fmsb)
LD_ref = "1kg"
df_keep_wcovar_pversion =  fread(paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/P_version_PRSCS/data/clean_output_pversion/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_QCset_prediction_wMDDstatus_wcovar.txt"))

df_keep_wcovar_update =  fread(paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_output_redo/ldref_"),LD_ref,"/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_QCset_prediction_wMDDstatus_wcovar.txt"))

df_keep_wcovar_OLD =  fread(paste0(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/OLD/clean_output/ldref_1kg/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_QCset_prediction_wMDDstatus_wcovar.txt")))

all.equal(df_keep_wcovar_OLD[,7:19],df_keep_wcovar_update[,7:19])

all.equal(df_keep_wcovar_OLD[,7:19],df_keep_wcovar_pversion[,7:19])

df_keep_wcovar_update$assessment_center = as.factor(df_keep_wcovar_update$assessment_center)

null_mod = glm(MDD ~ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+
                 sex+age+as.factor(assessment_center)+genotype_axiom,family = "binomial",data = df_keep_wcovar_update)
pheno_adj = null_mod$residuals
pheno_adj_norm = qnorm((rank(pheno_adj ,na.last="keep")-0.5)/sum(!is.na(pheno_adj )))


df = data.frame()
pred_wcov = df_keep_wcovar_update %>% as.data.frame()
k = 3
mod_hapmap = try(glm(MDD ~ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+
                       sex+age+as.factor(assessment_center)+genotype_axiom+scale(pred_wcov[,k]),family = "binomial",data = pred_wcov))

AIC(mod_hapmap)

k = 4
mod_cr95set = try(glm(MDD ~ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+
                        sex+age+as.factor(assessment_center)+genotype_axiom+scale(pred_wcov[,k]),family = "binomial",data = pred_wcov))
AIC(mod_cr95set)

k = 5
mod_lbf4 = try(glm(MDD ~ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+
                     sex+age+as.factor(assessment_center)+genotype_axiom+scale(pred_wcov[,k]),family = "binomial",data = pred_wcov))
k = 6
mod_twassig = try(glm(MDD ~ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+
                        sex+age+as.factor(assessment_center)+genotype_axiom+scale(pred_wcov[,k]),family = "binomial",data = pred_wcov))

saveRDS(mod_hapmap, mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/hapmap.rds"))
saveRDS(mod_cr95set, mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/cr95set.rds"))
saveRDS(mod_lbf4, mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/lbf4.rds"))
saveRDS(mod_twassig , mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/twassig.rds"))


mod_hapmap= readRDS( mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/hapmap.rds"))
mod_cr95set= readRDS( mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/cr95set.rds"))
mod_lbf4= readRDS( mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/lbf4.rds"))
mod_twassig = readRDS( mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/EUR_model/twassig.rds"))


library(pscl)
vuong(m1 = mod_hapmap, m2 = mod_cr95set)

h = summary(mod_hapmap)$coefficients["scale(pred_wcov[, k])",] %>% t() %>% as.data.frame()
c = summary(mod_cr95set)$coefficients["scale(pred_wcov[, k])",]%>% t() %>% as.data.frame()
l = summary(mod_lbf4)$coefficients["scale(pred_wcov[, k])",] %>% t()%>% as.data.frame()
t = summary(mod_twassig)$coefficients["scale(pred_wcov[, k])",]%>% t() %>% as.data.frame()

df = data.frame(type = c("hapmap", "cr95set", "logbf4", "twassig"),
                est = c(h$Estimate, c$Estimate, l$Estimate, t$Estimate),
                max = c(h$Estimate+h$`Std. Error`*1.96, c$Estimate+1.96*c$`Std. Error`, l$Estimate+1.96*l$`Std. Error`,
                        t$Estimate+1.96*l$`Std. Error`),
                min = c(h$Estimate-h$`Std. Error`*1.96, c$Estimate-1.96*c$`Std. Error`, l$Estimate-1.96*l$`Std. Error`,
                        t$Estimate-1.96*l$`Std. Error`))
df$type = fadctor()


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/AllEur_CoefCompare.png"),
       ggplot(df, aes(x = type,y = est, color = type))+geom_point()+
         geom_errorbar(aes(ymax = max, ymin = min))+scale_color_manual(values = c("maroon1","gold4", "dodgerblue3","darkgreen" ))+
         geom_text(aes(y = max+0.015, label = paste0("(", round(min,3)," , ", round(max,3), ")")),
                   color = "black", size = 5)+
         theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

                          axis.ticks.x = element_blank(),

                          text = element_text(size = 20, color = "black"),

                          axis.text = element_text(colour="black"),

                          strip.text.x = element_text(size = 22, colour = "black"),

                          legend.key.size = unit(1.5, "cm"))+
         labs(y = "Coefficient Estimates", x = "Set")+
         guides(color = "none")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.1, 0.1)))
       ,
       height = 8, width = 8)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/AllEur_CoefCompare_continuous_OR.png"),
       ggplot(df, aes(x = type,y = exp(est), color = type))+geom_point()+
         geom_errorbar(aes(ymax = exp(max), ymin = exp(min)))+scale_color_manual(values = c("maroon1","gold4", "dodgerblue3","darkgreen" ))+
         geom_text(aes(y = exp(max)+0.015, label = paste0("(", round(exp(min),3)," , ", round(exp(max),3), ")")),
                   color = "black", size = 5)+
         theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

                          axis.ticks.x = element_blank(),

                          text = element_text(size = 20, color = "black"),

                          axis.text = element_text(colour="black"),

                          strip.text.x = element_text(size = 22, colour = "black"),

                          legend.key.size = unit(1.5, "cm"))+
         labs(y = "Continuous Odds Ratio", x = "Set")+
         guides(color = "none")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.1, 0.1)))
       ,
       height = 8, width = 8)
