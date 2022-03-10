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

draw = cbind(assigned$score,assigned$phases)  #合并assigned$score列和assigned$phase列
colnames(draw)
attach(draw)

library(scatterplot3d)
scatterplot3d(G1, S, G2M, angle = 20,
              color = rainbow(3)[as.numeric(as.factor(assigned$phases))],
              grid = TRUE, box=FALSE)
detach(draw)

library(pheatmap)
cg=names(tail(sort(apply(dat,1,sd)),100))
n = t(scale(t(dat[cg,])))

#pheatmap(n, show_colnames =F, show_rownames=F)
library(pheatmap)
df$cellcycle = assigned$phases
ac=df
rownames(ac) = colnames(n)

pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,
         filename = 'all_cells_top_100_sd_all_infor.png')
dev.off()
head(ac)
table(ac[,c(1,5)])
