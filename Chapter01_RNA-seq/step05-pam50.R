rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input.Rdata')
a[1:4,1:4]    ##行：样本；列：基因
head(df) 

## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的细胞基因）
# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。

group_list=df$g
plate=df$plate
table(plate)

rownames(dat) = toupper(rownames(dat))
dat[1:4,1:4]

library(genefu)

if(T) {
  ddata = t(dat)
  ddata[1:4,1:4]
  s = colnames(ddata); head(s) ; tail(s)  ##把实验检测到的基因赋值给S
  library(org.Hs.eg.db)  ##人类基因信息的包
  s2g = toTable(org.Hs.egSYMBOL)
  g = s2g[match(s,s2g$symbol),1] ; head(g) ##取出实验检测到的基因所对应的基因名
  dannot = data.frame(probe = s,
                      "Gene.Symbol" = s,
                      "EntrezGene.ID" = g)
  View(dannot)
  ddata = ddata[, !is.na(dannot$EntrezGene.ID)] ##ID转换
  dim(ddata)
  dannot = dannot[!is.na(dannot$EntrezGene.ID),]     ##去除有NA的行，既剔除无对应的基因
  View(dannot)
  ddata[1:4,1:4]
  library(genefu)
  data("pam50.robust")
  data(pam50)
  data(scmod2.robust)
  # c("scmgene", "scmod1", "scmod2","pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow")
  s <- molecular.subtyping(sbt.model = "scmod2", data = ddata, annot = dannot, do.mapping=TRUE)
  
   table(s$subtype)
   tmp = as.data.frame(s$subtype)
   subtypes = as.character(s$subtype)
}
head(df)
df$subtypes = subtypes
table(df[,c(1,5)])

library(genefu)
 str(pam50)    ##其中centroids富含基因名
 """
 List of 7
 $ method.cor      : chr "spearman"
 $ method.centroids: chr "mean"
 $ std             : chr "none"
 $ rescale.q       : num 0.05
 $ mins            : num 5
 $ centroids       : num [1:50, 1:5] 0.718 0.537 -0.575 -0.119 0.3 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:50] "ACTR3B" "ANLN" "BAG1" "BCL2" ...
  .. ..$ : chr [1:5] "Basal" "Her2" "LumA" "LumB" ...
 $ centroids.map   :'data.frame':	50 obs. of  3 variables:
  ..$ probe          : chr [1:50] "ACTR3B" "ANLN" "BAG1" "BCL2" ...
  ..$ probe.centroids: chr [1:50] "ACTR3B" "ANLN" "BAG1" "BCL2" ...
  ..$ EntrezGene.ID  : int [1:50] 57180 54443 573 596 332 644 891 898 991 990 ...
 """

pam50genes = pam50$centroids.map[c(1,3)]
pam50genes[pam50genes$probe == 'CDCA1', 1] = 'NUF2'
pam50genes[pam50genes$probe == 'KNTC2', 1] = 'NDC80'
pam50genes[pam50genes$probe == 'ORC6L', 1] = 'ORC6'
x = dat
dim(x)   ###genecards存有人类所有基因ID

x = x[pam50genes$probe[pam50genes$probe %in% rownames(x)],]
table(group_list)

tmp = data.frame(group = group_list,
                 subtypes = subtypes)
rownames(tmp) = colnames(x)

rownames(x)

library(pheatmap)
pheatmap(x,show_rownames = T,show_colnames=F,
         annotation_col = tmp,
         filenames = 'ht_by_pam50.png')




   

                      
