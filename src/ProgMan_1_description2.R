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


# 3. 功能基因集在肿瘤 vs 正常表达量差异 -----
# 情况二：基因较多如genelist2，可先开展差异表达分析做一次筛选
count_dats <- read_rds(file.path(data_dir, 'TCGA.expr.dats.rds'))[[1]]
sample.group <- tibble(
  sample = colnames(count_dats)[-1],
  group = lapply(colnames(count_dats)[-1], function(x){ifelse(str_starts(str_sub(x,start = -3, end = -1), '0'), 'tumor', 'normal')}) %>% unlist()
)
# 调用函数开展差异表达分析
# ?run_deseq2
deseq2_res <-
  run_deseq2(count_dats, sample.group, control.group = 'normal')

# 调用函数标记显著上调和下调基因，并使用火山图展示全部基因差异分析结果
?plotutil.expr.volcano
sig_deseq2 <-
  plotutil.expr.volcano(
    deseq2.res = deseq2_res,
    log2FC.cutoff = c(-1, 1),
    padj.cutoff = 0.05,
    show.num = 5,
    color.pal = c("#fa6e57", "#4695d6")
  )
# 可保存deseq2结果!
View(sig_deseq2$deseq2.sig)
write_xlsx(sig_deseq2$deseq2.sig, file.path(out_dir, '2.TN_deseq2.xlsx'))
table(sig_deseq2$deseq2.sig$group)
#只查看功能基因集在肿瘤差异表达情况
genelist_deseq2 <- deseq2_res %>%
  filter(gene %in% genelist2) %>%
  mutate(group = case_when(
    padj <= 0.05 & log2FoldChange >= 1 ~ "up",
    padj <= 0.05 & log2FoldChange <= -1 ~ "down",
    T ~ "none"
  ))
table(genelist_deseq2$group)
# 火山图展示功能基因集差异分析结果
genelist_sig_deseq2 <- plotutil.expr.volcano(deseq2.res = genelist_deseq2,
                                             log2FC.cutoff = c(-1, 1),
                                             padj.cutoff = 0.05,
                                             show.num = 5,
                                             color.pal = c("#fa6e57", "#4695d6"))
genelist_sig_deseq2$volcano.plot
ggsave(file.path(out_dir, "3.genelist_sig_deseq2.pdf"), genelist_sig_deseq2$volcano.plot, width = 16, height = 8, units = "cm")

# 4. 单因素cox分析  -----
# 在功能基因中筛选与预后相关的基因
clin_all <- read_rds(file.path(data_dir, 'TCGA.clin.info.rds')) #含生存结局
log2tpm <- read_rds(file.path(data_dir, 'TCGA.expr.dats.rds'))[[2]] #表达矩阵
#只使用肿瘤样本开展预后分析
tumor_samples <- colnames(log2tpm)[str_detect(colnames(log2tpm), '01A|03A|05A')]
# 调用函数获取特定基因姐的cox分析输入
?get_gene_surv
surv_data <- get_gene_surv(log2tpm, clin_all, tumor_samples, genelist2)
# 将基因表达量转换为高低表达二分类
surv_data[,seq(2,ncol(surv_data)-2)]=apply(surv_data[,seq(2,ncol(surv_data)-2)],2,con2cat)
# 调用函数批量运行单因素cox回归
?batch_unicox
genelist_unicox <- batch_unicox(surv_data,
                                cols = colnames(surv_data)[2:(ncol(surv_data)-2)],
                                time.col = 'OS.time',
                                status.col = 'OS')
# 筛选与预后显著相关基因
genelist_unicox <- genelist_unicox %>%
  # 常用p值0.05或0.1为阈值筛选，可设置
  mutate(group = case_when(
    (PValue < 0.1) & (HR > 1) ~ 'poor', #基因高表达患者预后差
    (PValue < 0.1) & (HR < 1) ~ 'good', #基因高表达患者预后好
    TRUE ~ 'none'
  )) %>%
  arrange(factor(group, level = c('good', 'poor', 'none')),
          PValue,
          desc(HR)
  )
write_xlsx(genelist_unicox, file.path(out_dir, '2.genelist_unicox.xlsx'))

# 若预后相关基因较少，可用森林图展示结果，此处为示例
test_unicox <- genelist_unicox %>%
  filter(group %in% c('poor', 'good'))

# test_unicox <- genelist_unicox %>% filter(Var %in% genelist1)

test_unicox <- test_unicox[c(1:10, (nrow(test_unicox)-9):nrow(test_unicox)),] #示例20个基因
# 调用函数画forestplot
?plotutil.cox.forest
unicox_forest <- plotutil.cox.forest(test_unicox,
                                     gene.col = 'Var',
                                     HR.col = 'HR',
                                     HR.CI95.col = 'CI95',
                                     pval.col = 'PValue')
pdf(file.path(out_dir,'2.unicox_forest.pdf'), width = 8, height = 8)
unicox_forest
dev.off()

# 具有预后作用的差异表达基因（取上述两步交集基因）
DEGs <- genelist_deseq2 %>%
  filter(group %in% c('up', 'down'))
Prog_genes <- genelist_unicox %>%
  filter(group %in% c('poor', 'good')) %>%
  pull(Var)
Prog_DEGs <- intersect(DEGs$gene, Prog_genes)
# 韦恩图展示结果
venn_input <- list(DEGs = DEGs$gene,
                   UniCox = Prog_genes)
venn.plot <-
  plotutil.venn(
    venn_input,
    outfile =  file.path(out_dir, '2.Prog_DEGs_venn.pdf'),
    fontsize = 2,
    pic.size = 8
  )
