rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input.Rdata')
a[1:4,1:4]
head(df) 

group_list = df$g
plate = df$plate
table(plate)

a[1:4,1:4]
library(scran)
# https://mp.weixin.qq.com/s/nFSa5hXuKHrGu_othopbWQ
sce <- SingleCellExperiment(list(counts=dat)) 
#list() 创建列表

