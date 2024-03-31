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
# # 情况一：基因较少如genelist1，调用函数展示两组间表达差异
expr_dats <- read_rds(file.path(data_dir, 'TCGA.expr.dats.rds'))
log2tpm <- expr_dats$log2tpm
# 制作样本分组表
sample.group <- tibble(
  sample = colnames(log2tpm)[-1],
  group = lapply(colnames(log2tpm)[-1], function(x){ifelse(str_starts(str_sub(x,start = -3, end = -1), '0'), 'tumor', 'normal')}) %>% unlist()
)
sample.group$group %>% table()  #如果正常只有几个，需要借助GETx数据，非必需
head(sample.group) #查看分组表格式
tail(sample.group)
#调用函数画表达差异boxplot
?plotutil.group.expr.boxplot
expr_bxp <- plotutil.group.expr.boxplot(log2tpm,
                                        genelist1,
                                        # genelist1[1:10],
                                        #c("CASP4","NLRP2"), 
                                        sample.group,
                                         color_pal = c("#B22222","#A7D8DE"))
                                        # color_pal = get_palette('jco',2))
expr_bxp 
 # 存储画图文件
ggsave(file.path(out_dir, '1.genelist1.TN.expr.pdf'),
       expr_bxp, width = 14, height = 6)

