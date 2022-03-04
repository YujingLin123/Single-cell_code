rm(list = ls())
options(stringsAsFactors = F)

if (F) {
  load(file = '../input.Rdata')
  a[1:4, 1:4]
  head(metadata)
  ## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
  # 注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。
  group_list = metadata$g
