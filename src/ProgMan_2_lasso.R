## Developed by WoLin DM Team ###
### Used for prognosis dataming ###

# 该代码主要用于文章的第一部分描述性分析

library(renv) #加载renv
renv::activate("./") #激活当前项目，用renv管理环境
library(ProgMan)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(readxl)
library(writexl)
library(data.table)
library(ggpubr)

# 5. lasso-cox建模 -----
genelist <- genelist1 #初始基因集-较少数目未筛选
genelist <- Prog_DEGs #经差异表达分析和单因素cox分析后筛选得到的基因集

clin_all <- read_rds(file.path(data_dir, 'TCGA.clin.info.rds')) #生存结局
log2tpm <- read_rds(file.path(data_dir, 'TCGA.expr.dats.rds'))[[2]] #基因表达数据
#只使用肿瘤样本开展预后分析
tumor_samples <- colnames(log2tpm)[str_detect(colnames(log2tpm), '01A|03A|05A')]
# 调用函数获取cox分析输入,同上面单因素cox分析
surv_data <- get_gene_surv(log2tpm, clin_all, tumor_samples, genelist)
# 调用函数进行lasso-cox建模
?lasso_cox
lasso_res <- lasso_cox(surv_data, time.col = "OS.time",
                       status.col = "OS", nfold = 20, seed = 10,
                       plot.dir = file.path(out_dir, '3.Lasso_plots'))
lasso_res
# lasso_cox输出一个列表，包含筛选后的基因及其系数，以及模型的lambda值。并且在画图目录下输出两个筛选过程图
# 调用函数计算模型C-index指标
?C_index
C_index(lasso_res$coef$gene, surv.dats = surv_data,
        time.col = "OS.time",status.col = "OS")
# 不同seed下结果可能不同，可多次筛选，C_index更高者更优，不设置seed则随机出结果
lasso_res <- lasso_cox(surv_data, time.col = "OS.time",status.col = "OS",
                       nfold = 20, seed = 200,
                       plot.dir = file.path(out_dir, '3.Lasso_plots'))
C_index(lasso_res$coef$gene, surv.dats = surv_data, time.col = "OS.time",status.col = "OS")

# 将最终lasso_cox模型结果写入结果文件夹
write_xlsx(lasso_res,file.path(out_dir, '3.lasso_res.xlsx'))
# 该模型为risk score模型，用于后续区分预后及下游分析

# 森林图展示多基因模型：基因和系数
lasso_res <- read_xlsx(file.path(out_dir, '3.lasso_res.xlsx'))
# 由于森林图一般同时展示每个基因的风险值，读入之前unicox结果
genelist_unicox <- read_xlsx(file.path(out_dir, '2.genelist_unicox.xlsx'))
# 森林图输入
lasso_res_unicox <- lasso_res %>%
  left_join(genelist_unicox, by = c('gene' = 'Var')) %>%
  arrange(HR)
# 调用函数生成森林图
lasso_forest <- plotutil.cox.forest(lasso_res_unicox,gene.col = 'gene',
                                    HR.col = 'HR', HR.CI95.col = 'CI95',
                                    pval.col = 'PValue', coef.col = 'coef')
pdf(file.path(out_dir,'3.lasso_forest.pdf'), width = 8, height = 8)
lasso_forest
dev.off()

