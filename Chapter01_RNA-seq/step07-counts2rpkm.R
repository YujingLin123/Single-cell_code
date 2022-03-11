setwd("J:/C盘/知识点汇总/单细胞转录组/scRNA_smart_seq2-master/section01-RNA-seq")

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input.Rdata')
a[1:4,1:4]
head(df)

## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
# 注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
?TxDb

## 下面是定义基因长度为 非冗余exon长度之和
if(F){
  exon_txdb = exons(txdb)
  genes_txdb = genes(txdb)
  genes_txdb
  ?GRanges
  o = findOverlaps(exon_txdb,genes_txdb)
  o
  ## exon - 1 : chr1 4807893-4807982
  ## 1        6523
  #  genes_txdb[6523]  # chr1 4807893-4846735 , 18777
  
  t1 = exon_txdb[queryHits(o)]   ##queryHits表示交集的A的位点序列
  t2 = genes_txdb[subjectHits(o)]    ##subjectHits表示交集的B的位点序列
  t1 = as.data.frame(t1)
  t1$geneid = mcols(t2)[,1]  ##mcols(x, use.names=FALSE), mcols(x) <- value: Get or set the metadata columns. value can be NULL, or a data.frame-like object (i.e. DataFrame or data.frame) holding element-wise metadata.
  
  # 如果觉得速度不够，就参考R语言实现并行计算
  # http://www.bio-info-trainee.com/956.html
  #lapply : 遍历列表向量内的每个元素，并且使用指定函数来对其元素进行处理。返回列表向量。
  #函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组;
  #它的返回值是一个列表，代表分组变量每个水平的观测。
  
  #一个基因的外显子=所有exon的加和-该基因所有exon交集的长度
  g_1 = lapply(split(t1,t1$geneid),function(x){
    # x = split(t1,t1$geneid)[[1]] ,split函数使t1的data.frame根据geneid对t1进行拆分，根据拆分的每个基因进行第二个函数的操作，如此循环下去
    head(x)
    tmp = apply(x,1,function(y){
      y[2]:y[3]    ##y[2]:start,y[3]:end，取出每个基因的外显子的坐标
     })
    length(unique(unlist(tmp)))   ##unlist就是返回该data.frame的非list元素，统计每个基因的非冗余外显子的长度
    #sum(x[,4])
  })
  head(g_1)
  g_1 = data.frame(gene_id = names(g_1),length=as.numeric(g_1))
  
  save(g_1, file = 'step7-g_1.Rdata')
  }
## 下面是定义基因长度为 最长转录本长度
if(F){
  t_1 = transcriptLengths(txdb)
  head(t_1)
  t_1 = na.omit(t_1)
  head(t_1)
  t_1 = t_1[order(t_1$gene_id,t_1$tx_len,decreasing=T),]
  head(t_1)
  str(t_1)
  t_1 = t_1[!duplicated(t_1$gene_id),]
  head(t_1)
  g_l = t_1[,c(3,5)]
}
head(g_l)
library(org.Mm.eg.db)
s2g = toTable(org.Mm.egSYMBOL)
head(s2g)
g_l = merge(g_l,s2g,by='gene_id')   #把g_l,s2g两个数据框以'gene_id'为连接进行拼接

# 参考counts2rpkm，定义基因长度为非冗余CDS之和
# http://www.bio-info-trainee.com/3298.html  
a[1:4,1:4]
ng = intersect(rownames(a),g_l$symbol)   #取a数据框的行名与g_l数据框的symbol列的交集

# 有了counts矩阵和对应的基因长度信息，就很容易进行各种计算了：
exprSet = a[ng,]
lengths=g_l[match(ng,g_l$symbol),2]
head(lengths)
head(rownames(exprSet))
# http://www.biotrainee.com/thread-1791-1-1.html
exprSet[1:4,1:4]    ##获得每个细胞每个基因的count数
total_count<- colSums(exprSet)    ##计算获得每个细胞的所有的count数
head(total_count)
head(lengths)
total_count[4]
lengths[1]
#rpkm: 1*10^9/(1122*121297)
rpkm <- t(do.call(rbind,
                  lapply(1:length(total_count),
                         function(i){
  10^9*exprSet[,i]/lengths/total_count[i]
})))
rpkm[1:4,1:4]
# 下面可以比较一下 自己根据counts值算出来的RPKM和作者提供的RPKM区别。
a=read.table('../GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt.gz',
             header = T ,sep = '\t')
# 每次都要检测数据
a[1:4,1:4]
rpkm_paper=a[ng,] 
rpkm_paper[1:4,1:4]

rpkm[1:4,1:4]

  
  
  
