## Developed by WoLin DM Team ###
### Used for prognosis dataming ###

# 该代码主要用于文章的第三部分模型验证

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

# 6.建模后验证 ----
lasso_res <- read_xlsx(file.path(out_dir, '3.lasso_res.xlsx')) #模型基因和系数
log2tpm <- read_rds(file.path(data_dir, 'TCGA.expr.dats.rds'))[[2]] #表达数据
clin_all <- read_rds(file.path(data_dir, 'TCGA.clin.info.rds')) #生存情况
## 6.0 根据模型计算每个样本的risk score ----
# 只针对肿瘤样本
tumor_samples <- colnames(log2tpm)[str_detect(colnames(log2tpm), '01A|03A|05A')]
#调用函数获取模型基因表达和样本生存数据
surv_data <- get_gene_surv(log2tpm, clin_all,
                           tumor_samples,
                           genelist = lasso_res$gene)
# 每个样本risk score和risk group(需指定高低分组的分位值，常用中位值分组-分位于50，也可指定高风险为top四分之一，即分位于25)
lasso_risk <- risk.score.calulator(lasso.res = lasso_res,
                                   expr.surv = surv_data,
                                   group.cutoff = 50)
write_xlsx(lasso_risk, file.path(out_dir, '4.lasso_risk.xlsx'))
## 6.1 模型基因+riskscore+生存情况展示图 ----
# 调用函数生成模型展示图
risk.plot <- plotutil.lasso.model(lasso.risk = lasso_risk,
                                  time.col = "OS.time",
                                  status.col = "OS",
                                  color_pal1 = c("#DF8F44FF", "#374E55FF"),
                                  color_pal2 = c("#79BEDB", "#FD7013")
)
risk.plot
pdf(file.path(out_dir, '4.lasso_model_plot.pdf'), height=7, width=7, onefile = FALSE)
risk.plot
dev.off()

## 6.2 TCGA训练集KM曲线 -----
lasso_risk <- read_xlsx(file.path(out_dir, '4.lasso_risk.xlsx'))
# 调用函数针对n年生存率画KM曲线
p1 <- plotutil.nyear.KM.cruve(df.surv = lasso_risk,
                              nyear = 10, #一般3年，5年或10年
                              legend.ref = 'high_risk',
                              hr.ref = 'low_risk',
                              time.col = "OS.time",
                              status.col = "OS",
                              survival.var = 'OS',
                              surfit.var = c('risk_group'),
                              leg.pos = c(0.8, 0.85),
                              color.pal = c("#ED0000FF", "#00468BFF"))
p1
pdf(file.path(out_dir, '4.TCGA_KM_10year.pdf'),height = 6, width = 6, onefile = FALSE)
p1
dev.off()

## 6.3 独立验证集读入、计算risk score和画KM曲线 ----
# 每个癌种验证集可提供，形式为rds数据，大数据使用rds形式保存便于读入
# 读入数据，查看数据格式，sample列+各基因列+OS信息两列
testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE12945_expr_surv.rds'))
testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE17536_expr_surv.rds'))
testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE17537_expr_surv.rds'))
testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE29621_expr_surv.rds'))
testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE39582_expr_surv.rds'))
testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE41258_expr_surv.rds'))
# testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE71187_expr_surv.rds'))
# testset <- read_rds(file.path(wd, '0.rawdata/CEO_Validation_Set/GSE87211_expr_surv.rds'))

# 多个验证集逐一计算risk score，注意查看是否缺基因
# 调用函数根据lasso模型计算每个样本的risk score
testset.risk <- risk.score.calulator(lasso.res = lasso_res,
                                     expr.surv = testset, #输入验证集
                                     group.cutoff = 50
)

setdiff(lasso_res$gene,colnames(testset))

# "SUCLG2P2"为假基因，在多个GEO数据集中表达数据不存在，因此优先推荐使用"protein_coding"基因建模，
# 如前所述，若在多个验证队列中均存在模型基因表达不存在，则应重新建模
# 参考前面从genelist2到genelist3，应使用genelist3建模，此处暂将模型中SUCLG2P2剔除做验证集的演示
lasso_res <- lasso_res %>%
  filter(gene != 'SUCLG2P2')
# 剔除后可选的验证集有GSE17536，GSE17537，GSE29621，GSE39582，其他验证集损失基因更多，不可用
# 对每个数据集计算risk score并画KM曲线
testsets <-
  c(
    file.path(wd, '0.rawdata/CEO_Validation_Set/GSE17536_expr_surv.rds'),
    file.path(wd, '0.rawdata/CEO_Validation_Set/GSE17537_expr_surv.rds'),
    file.path(wd, '0.rawdata/CEO_Validation_Set/GSE29621_expr_surv.rds'),
    file.path(wd, '0.rawdata/CEO_Validation_Set/GSE39582_expr_surv.rds')
  )
# 批量每一个验证集画KM曲线，并输出到结果文件夹
for (i in seq_along(testsets)) {
  testset <- read_rds(testsets[i])
  testset.risk <- risk.score.calulator(lasso.res = lasso_res, #剔除一个基因后，用于演示的模型基因
                                       expr.surv = testset, #输入验证集
                                       group.cutoff = 50
  ) %>%
    dplyr::select(sample, lasso_res$gene, OS.time, OS.status, risk_score, risk_group)
  write_xlsx(testset.risk, file.path(out_dir, paste0('4.testset', as.character(i), '_lasso_risk.xlsx')))
  p1 <- plotutil.nyear.KM.cruve(df.surv = testset.risk,
                                nyear = 10,
                                legend.ref = 'high_risk',
                                hr.ref = 'low_risk',
                                time.col = "OS.time",
                                status.col = "OS.status",
                                survival.var = 'OS',
                                surfit.var = c('risk_group'),
                                leg.pos = c(0.75, 0.85),
                                color.pal = c("#ED0000FF", "#00468BFF"))
  pdf(file.path(out_dir, paste0('4.testset', as.character(i), '_KM_10year.pdf')),height = 6, width = 6, onefile = FALSE)
  print(p1)
  dev.off()
}
# 预期结果：在独立验证集中，高低风险KM曲线分得很开，p<0.05
# 示例：在独立验证集中模型表现不好，应重新训练模型（可尝试调整基因集，分组分位值等）
# 此处不影响后续代码演示所以暂不重新训练

## 6.4 训练集和验证集时间相关AUC ----
# 训练集数据
lasso_risk <- read_xlsx(file.path(out_dir, '4.lasso_risk.xlsx'))
# 调用函数生成时间相关ROC曲线，并显示AUC
p1 <- plotutil.time.AUC(df.surv = lasso_risk,
                        surv.year = c(1,3,5), #预测n年生存率的AUC
                        time.col = 'OS.time',
                        status.col = 'OS',
                        marker.var = 'risk_score', #预测因子
                        color.pal = c("#BC3C29FF", "#0072B5FF", "#E18727FF")
)
p1
ggsave(filename = file.path(out_dir, '4.TCGA_AUC.pdf'), p1, height = 5, width = 5)
# 验证集数据获取AUC
for (i in seq_along(testsets)) {
  testset.risk <- read_xlsx(file.path(out_dir, paste0('4.testset', as.character(i), '_lasso_risk.xlsx')))
  p2 <- plotutil.time.AUC(df.surv = testset.risk,
                          surv.year = c(1,3,5),
                          time.col = "OS.time",
                          status.col = "OS.status",
                          marker.var = 'risk_score',
                          color.pal = c("#9b2226", "#2c6e49", "#7400b8")
  )
  ggsave(filename = file.path(out_dir, paste0('4.testset', as.character(i), '_AUC.pdf')), p2, height = 5, width = 5)
}

## 6.5 多因素cox风险回归 ----
# 一般针对TCGA训练集，需要临床信息较完整，用于判断分子预后模型是否独立于临床因素
lasso_risk <- read_xlsx(file.path(out_dir, '4.lasso_risk.xlsx')) # 样本对应的risk score和risk group
clin_all <- read_rds(file.path(data_dir, 'TCGA.clin.info.rds')) #临床信息
# 多因素cox输入数据框，列包含样本、risk score和risk group、样本对应的临床信息列、生存时间和状态信息
multicox.input <- lasso_risk %>%
  dplyr::select(sample, risk_score, risk_group) %>%
  left_join(clin_all) %>%
  mutate(tumor.stage = ifelse(tumor.stage %in% c('I', 'II', 'III', 'IV'),
                              tumor.stage, NA_character_))
# 可打开multicox.input表，选取纳入多因素回归分析的临床因素
clin_factors <- c("age","LVI","MSI.status" ,"gender" ,"tumor.stage","race")
multic_res <- Multicox(
  df.surv = multicox.input,
  clin.factors = clin_factors,
  marker.var = 'risk_score', #关注的预后因子（分子预后模型结果）
  time.col = 'OS.time',
  status.col = 'OS'
)
# 运行时出现warning，表示多因素风险回归无法拟合，纳入的临床变量不合适，打开结果查看
multic_res$tbl
multic_res$plot
# 结果表中Characteristics列代表变量，连续变量直接显示变量名，分类变量显示变量名+类别名如tumor.stageII，表示tumor.stage这一临床因素中，tumor.stage II相比tumor.stage I(基准类别)的风险

# 可以看到MSI.status,tumor.stage和race这几个临床因素的CI95列出现了Inf或者NA值，因此这些临床因素纳入多因素cox回归模型会导致模型无法拟合，一般情况去掉这几个临床因素再做多因素cox回归
clin_factors <- c("age","LVI","gender")
multic_res <- Multicox(
  df.surv = multicox.input,
  clin.factors = clin_factors,
  marker.var = 'risk_score',
  time.col = 'OS.time',
  status.col = 'OS'
)
multic_res$tbl #查看结果表
multic_res$plot #查看图，可输出到结果文件夹

# 注意：一般tumor.stage是需要纳入多因素cox回归的，若由于分多类别导致模型无法拟合时
# 将tumor.stage这种有程度递增关系的分类变量转换为等级变量是合理的
multicox.input[['tumor.stage']] <- as.numeric(factor(multicox.input[['tumor.stage']],ordered = T))
multicox.input[['tumor.stage']] %>% table() #查看转换后tumor.stage表示
# 重新拟合回归模型
clin_factors <- c("age","gender","tumor.stage")
multic_res <- Multicox(
  df.surv = multicox.input,
  clin.factors = clin_factors,
  marker.var = 'risk_score',
  time.col = 'OS.time',
  status.col = 'OS'
)
multic_res$tbl #查看结果表
multic_res$plot #查看图，可输出到结果文件夹

#另外，分子预后模型得到的预后因子除了risk score外，还有risk group，
# 因此也可将risk group作为目标预后因子，查看是否相比其他临床因素独立预后
clin_factors <- c("age","LVI","gender","tumor.stage")
# 先指定是high_risk vs low_risk
multicox.input$risk_group <- factor(multicox.input$risk_group, levels = c('low_risk', 'high_risk'))
multic_res <- Multicox(
  df.surv = multicox.input,
  clin.factors = clin_factors,
  marker.var = 'risk_group',
  #修改关注的预后因子
  time.col = 'OS.time',
  status.col = 'OS'
)
multic_res$tbl #查看结果表
write_xlsx(multic_res$tbl, file.path(out_dir, '4.multicox_res.xlsx'))
multic_res$plot #查看图，可输出到结果文件夹
# 同时保存multicox.input用于后续复用
write_xlsx(multicox.input,file.path(out_dir, '4.multicox_input.xlsx'))

## 6.6 nomogram构建+校准曲线 ----
multic_df <- read_xlsx(file.path(out_dir, '4.multicox_input.xlsx'))
# multic_df <- read_xlsx('~/Downloads/4.multicox_input.xlsx')
# 同样可以基于risk_score或risk_group构建nomogram，其他纳入的临床因素应采用上面多因素cox模型中临床因素，纳入不合理临床因素（如多因素cox分析中导致模型无法拟合的因素）会报错
clin_factors <- c("age","LVI","gender","tumor.stage")
nomo.risk <-
  plotutil.nomogram(
    df.surv = multic_df,
    surv.year = c(3,5),
    clin.factors = clin_factors,
    marker.var = 'risk_score',
    time.col = 'OS.time',
    status.col = 'OS',
    out_dir = out_dir, #设置图的输出目录
    nomo.h = 7,
    nomo.w = 9.5 #设置nomogram图输出的宽和高，保证合理显示
  )
# 函数输出1：输出nomogram图和校准曲线到文件夹
# 函数输出2：打印Nomogram模型预测值的C-index，用于评价nomogram预测性能
# "Nomogram C-index; 0.783594084869267"
# 函数输出3：Nomogram模型预测值，可保存文件
write_xlsx(nomo.risk,file.path(out_dir, '5.nomogram_risk.xlsx'))

## 6.7 nomogram模型与临床因素/分子预后riskscore之间的AUC比较 ----
multic_df <- read_xlsx(file.path(out_dir, '4.multicox_input.xlsx'))
auc.factors <- c("age","LVI","gender","tumor.stage", "risk_score")
# 用于比较的因素为纳入nomogram的临床因素和risk score
auc_plotlist <- plotutil.nomo.compare.AUC(df.surv = multic_df,
                                          auc.factors,
                                          surv.year = c(3,5,7),  #想要比较n年生存率的预测性能
                                          time.col = 'OS.time',
                                          status.col = 'OS')
auc_plotlist #查看输出，可保存到输出
ggsave(filename = file.path(out_dir, '5.nomogram_risk_roc.pdf'),auc_plotlist,width = 12,height = 6)
