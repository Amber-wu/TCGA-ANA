### Developed by WoLin DM Team ###
### Used for prognosis dataming ###

# 该代码主要用于原始数据的整理和清洗

# library(renv) #加载renv
# renv::activate("./") #激活当前项目，用renv管理环境
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

# 设置工作目录 ----
# 一般即为项目目录，之后所有数据和结果都位于该目录下
wd <- './ProgTest/'
dir.create(wd, recursive = T)
# 推荐前期下载的TCGA表达文件等原始数据文件存入wd下的0.rawdata文件夹
# 设置清洗后的数据存储目录
data_dir <- file.path(wd, '1.cleandata')
dir.create(data_dir, recursive = T)
# 设置分析结果输出目录
out_dir <- file.path(wd, '2.results')
dir.create(out_dir, recursive = T)

# 1. TCGA数据准备 ----
## 1.1 表达数据处理 -----
# 输入三个下载文件，输出count和log2(tpm+1)矩阵
?Xena.expr.dats()
# 输入文件
probemap.file <- file.path(wd, '0.rawdata/TCGA/gencode.v22.annotation.gene.probeMap')
count.file <- file.path(wd, '0.rawdata/TCGA/TCGA-COAD.htseq_counts.tsv.gz')
fpkm.file <- file.path(wd, '0.rawdata/TCGA/TCGA-COAD.htseq_fpkm.tsv.gz')

# 调用函数清洗表达矩阵
expr_dats <- Xena.expr.dats(count.file, fpkm.file, probemap.file)
# count矩阵
View(expr_dats$count_df)
# TPM矩阵
View(expr_dats$log2tpm)
# 存储
write_rds(expr_dats, file.path(data_dir, 'TCGA.expr.dats.rds'))

## 1.2 临床数据处理 -----
# 导入临床数据
clin.file <- file.path(wd, '0.rawdata/TCGA/TCGA-COAD.GDC_phenotype.tsv.gz')
#clin.file <- file.path(wd, '0.rawdata/TCGA/TCGA-LAML.GDC_phenotype.tsv.gz')
# 选择部分重点临床信息列——依癌种而定
clin_dats <- read_tsv(clin.file) %>%
  dplyr::select(submitter_id.samples, submitter_id,
                age_at_initial_pathologic_diagnosis,
                lymphatic_invasion, microsatellite_instability,
                pathologic_M, pathologic_N, pathologic_T,
                perineural_invasion_present,
                preoperative_pretreatment_cea_level,
                gender.demographic, race.demographic,
                tumor_stage.diagnoses) %>%
  rename(sample = submitter_id.samples,
         id = submitter_id,
         age = age_at_initial_pathologic_diagnosis,
         LVI = lymphatic_invasion,
         MSI.status = microsatellite_instability,
         pathologic.M = pathologic_M,
         pathologic.N = pathologic_N,
         pathologic.T = pathologic_T,
         PNI = perineural_invasion_present,
         CEA = preoperative_pretreatment_cea_level,
         gender = gender.demographic,
         race = race.demographic,
         stage = tumor_stage.diagnoses
  )
# 读入生存数据
surv.file <- file.path(wd, '0.rawdata/TCGA/TCGA-COAD.survival.tsv')
surv_dats <- read_tsv(surv.file)
# 整合临床和生存数据
expr_dats <- read_rds(file.path(data_dir, 'TCGA.expr.dats.rds'))[[1]]
clin_all <- clin_dats %>%
  left_join(surv_dats, by = 'sample') %>%
  filter(sample %in% colnames(expr_dats)) %>%  #具有表达数据的样本
  distinct(id, .keep_all = TRUE) %>%  #去重一个患者的多个样本
  mutate(
    OS.status = case_when(
      is.na(OS) ~ 'missing',
      OS == 1 ~ 'dead',
      OS == 0 ~ 'alive'
    ), #将生存状态变为分类变量
    tumor.stage = case_when(
      str_detect(stage, 'iv') ~ 'IV',
      str_detect(stage, 'iii') ~ 'III',
      str_detect(stage, 'ii') ~ 'II',
      str_detect(stage, 'i') ~ 'I',
      TRUE ~ 'not reported'
    ) # 整合肿瘤分期信息
  )
# 存储临床信息用于后续分析
write_rds(clin_all, file.path(data_dir, 'TCGA.clin.info.rds'))


## 答疑补充代码
table(clin_all$stage)
clin_all <- clin_all %>%
  # filter(stage == "stage iii") %>%
  filter(stage %in% c("stage iii","stage iiia","stage iiib","stage iiic"),
         age > 50)

# 统计用于分析的患者，输出table1
summary_get_table_one(
  tb = clin_all,
  var_all = c("age", "CEA","gender","race", "LVI", "PNI", "tumor.stage", "OS.status"), # 全部变量
  var_cat = c("gender","race", "LVI", "PNI", "tumor.stage", "OS.status"), # 分类变量
  var_con = c("age", "CEA" ), #连续变量
  output_xlsx = file.path(out_dir, '1.table1.xlsx')
)

# 2.1 功能基因集读入01 ----
# 一般来源于文献收集，第一列为基因名的表格文件，有表头
genelist1 <- read_xlsx(file.path(wd, '0.rawdata/genelist/Pyroptosis_33_genes.xlsx'))[[1]] %>% unique()

# 2.2 功能基因集读入02 ----
# 用关键词从数据库筛选基因集
keywords <- c('necrosis','apoptosis', 'pyroptosis')
# keywords <- c("epithelial")
# MsigDB功能基因集合，可下载多个文件输入下面函数
gmt1 <- file.path(wd, '0.rawdata/MsigDB/c2.cp.v2023.2.Hs.symbols.gmt')
gmt2 <- file.path(wd, '0.rawdata/MsigDB/c5.go.v2023.2.Hs.symbols.gmt')
gmt3 <- file.path(wd, '0.rawdata/MsigDB/h.all.v2023.2.Hs.symbols.gmt')

# 调用函数获取基因集
?geneset.selector
# genesets <- geneset.selector(keywords, pattern = '&', gmt1, gmt2, gmt3)
genesets <- geneset.selector(keywords, pattern = '|', gmt1, gmt2, gmt3)
genesets <- geneset.selector(keywords, pattern = '|', gmt3)
genelist2 <- genesets %>%
  distinct(gene) %>%
  pull(gene) #去重得到目标基因集

# 非必需，可筛选蛋白编码基因用于建模，后续GEO芯片数据集更方便验证
gene_type <- read.csv(system.file("extdata", "gencode.v22.hg38.gene.type.csv",package = "ProgMan")) #文件已包含在ProgMan包中
protein_coding <- gene_type %>%
  filter(gene_type == 'protein_coding')
genelist3 <- intersect(genelist2, protein_coding$gene_name) #3738 to 3525
