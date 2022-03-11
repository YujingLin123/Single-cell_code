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
    # x = split(t1,t1$geneid)[[1]]
    head(x)
  
  
  
