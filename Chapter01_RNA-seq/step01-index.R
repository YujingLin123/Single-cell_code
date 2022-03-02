rm(list = ls())  ##魔幻操作，一键清空
Sys.setenv(R_MAX_NUM_DLLS=999)  ##Sys.setenv 修改环境设置，R的namespace 是有上限的，如果导入包时超过这个上次就会报错，R_MAX_NUM_DLLS可以修改这个上限
options(stringAsFactors = F)  ##options:允许用户对工作空间进行全局设置，stringsAsFactors放置R自动把字符串string的列辨认成factor

options()$repos  ##查看使用install.packages安装时的默认镜像
options()$BioC_mirror  ##查看使用bioconductor的默认镜像
options(BioC_mirror = "https://mirrors.ustc.edu.cn/bios")  ##指定镜像，这个时中国科技大学镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))  ##指定install.packages安装镜像，这个是清华镜像
options()$repos
options()$BioC_mirror

if ( ! requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

if(F) {
  BiocManager::install(c( 'scran'),ask = F,update = F)
  
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene",ask = F,update = F)
  BiocManager::install("org.Mm.eg.db",ask = F,update = F)
  BiocManager::install("genefu",ask = F,update = F)
  BiocManager::install("org.Hs.eg.db",ask = F,update = F)
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene",ask = F,update = F)
  install.packages("ggfortify")
  install.packages("FactoMineR")
  install.packages("factoextra")
  install.packages("ggfortify")
  }

if (T) {
  a = read.table('../GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz', header = T, sep="\t")
  a[1:6,1:4] 
  dat = a[apply(a, 1, function(x) sum(x>1) > floor(ncol(a)/50)), ]  ##读取RNA-seq 的counts定量结果， 表达矩阵需要进行简单的过滤
  #筛选表达量合格的行（基因）， 列（细胞）数不变
    ##ncol()返回矩阵的列数值；floor()四舍五入取整数；
  ##function() 定义一个函数；sum()求和
  #上面的apply()指令代表对矩阵a进行行计算，判断每行表达量>1的样本总个数，并筛选出细胞表达量合格的基因（行）
  #第一个参数是指要参与计算的矩阵——a
  #第二个参数是指按行计算还是按列计算，1——表示按行计算，2——按列计算；
  #第三个参数是指具体的运算参数,定义一个函数x（即表达量为x）
  
  #对每行中x>1的列（即样本数）求和，即得出的是每行中表达量大于1的样本数，
  #然后再筛选出大于floor(ncol(a)/50)的行，这样的行（基因）的细胞表达量才算合格
  #因为2%的细胞有表达量，所以对于768个细胞样本，每个行(基因)在细胞中的表达至少要有15.36（约等于15）个样本表达才算合格
  ## 2 % 的细胞有表达量

dat[1:4, 1:4]
sum(dat[,3])
#  0610007P14Rik in  SS2_15_0048_A5 
  log2(18*1000000/sum(dat[,3])+1)
  ## 18 -- > 6.459884  ## SS2_15_0048_A5
  dat=log2(edgeR::cpm(dat)+1)
  ##归一化的一种选择，这里是CPM(count-per-million，每百万碱基中每个转录本的count值)
  ###CPM只对read count相对总reads数做了数量的均一化，去除文库大小差异。
  dat[1:4,1:4]            

  # 总结：
  # 
  # - dist函数计算行与行（样本）之间的距离
  # - cor函数计算列与列（样本）之间的相关性
  # - scale函数默认对每一列（样本）内部归一化
  # - 计算dist之前，最好是对每一个样本（列）进行scale一下
  # 
 #层次聚类，因为近 800细胞，非常耗时。
  hc=hclust(dist(t(dat))) ##样本间层次聚类
                
  # 原始表达矩阵转置后，细胞在行，所以计算的是细胞与细胞之间的距离。
  clus = cutree(hc, 4) #对hclust()函数的聚类结果进行剪枝，即选择输出指定类别数的系谱聚类结果。
  group_list= as.factor(clus) ##转换为因子属性
  table(group_list) ##统计频数

  #提取批次信息                
  colnames(dat) #取列名
  library(stringr)
  plate=str_split(colnames(dat),'_',simplify = T)[,3] #取列名，以'_'号分割，提取第三列。
  #str_split()函数可以分割字符串
  table(plate)
                
n_g = apply(a, 2, function(x) sum(x>1))  ##统计每个样本有表达的有多少行（基因）
#reads数量大于1的那些基因为有表达，一般来说单细胞转录组过半数的基因是不会表达的。
            
df = data.frame(g = group_list, plate = plate, n_g = n_g)      ##新建数据框（细胞的属性信息）
#样本为行名，列分别为：样本分类信息，样本分组，样本表达的基因数   
df$all = 'all'            
metadata=df
save(dat, metadata, file = '../input_rpkm.Rdata') ##保持a, dat, df这变量到上级目录
}
            
load(file = '../input.Rdata') ##从上级目录载入input.Rdata

## 每次载入以前的变量，都是可以简单检查一下。
a[1:4,1:4] 
dat[1:4,1:4] 
head(df)

load(file = '../input_rpkm.Rdata') ##从上级目录载入 input_rpkm.Rdata

## 每次载入以前的变量，都是可以简单检查一下。
#a[1:4,1:4] 
dat[1:4,1:4] 
head(metadata)
            
