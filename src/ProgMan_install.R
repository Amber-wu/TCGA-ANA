#### 配置电脑环境以及必要的package ####
#使用renv来进行环境管理，先安装renv包
library(renv) #加载renv
renv::activate("./") #激活当前项目，用renv管理环境

# 安装包
renv::install("readxl")
# 测试是否安装成功
library(readxl)

# 安装包
renv::install("writexl")
# 测试是否安装成功
library(writexl)

# 安装课程核心包（安装的同时会把其中依赖的其他包一起安装）
renv::install('./ProgMan_0.0.0.9000.tar.gz')
renv::install('MatrixModels')

#加载多基因预后模型R包
library(ProgMan)
#加载常用数据处理R包
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(data.table)
library(ggpubr)

# 如果上述package在加载的过程中均为出现error提示，后续课程的代码则可以正常运行