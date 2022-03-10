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

library(org.Mm.eg.db)
mm.pairs <- readRDS(system.file("exdata","mouse_cycle_markers.rds",
                                package = "scran"))
ensembl <- mapIds(org.Mm.eg.db, keys = rownames(sce),
                  keytype = "SYMBOL", column = "ENSEMBL")
#取探针名创建一个向量
#rownames(sce) 取行名（即实验检测到的基因）

##if(F)就是将{}代码注释掉，需要单独运行
if(F){
  assigned <- cyclone(sce,pairs = mm.pairs, gene.names = ensembl)
  save(assigned, file = 'cell_cycle_assigned.Rdata')
}

load(file = 'cell_cycle_assigned.Rdata')
head(assigned$scores)
table(assigned$phases)

draw = cbind()
