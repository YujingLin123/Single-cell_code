rm(list = ls())  ##魔幻操作，一键清空
Sys.setenv(R_MAX_NUM_DLLS=999)  ##Sys.setenv 修改环境设置，R的namespace 是有上限的，如果导入包时超过这个上次就会报错，R_MAX_NUM_DLLS可以修改这个上限
options(stringAsFactors = F)  ##options:允许用户对工作空间进行全局设置，stringsAsFactors放置R自动把字符串string的列辨认成factor

options()$repos  ##查看使用install.packages安装时的默认镜像
options()$BioC_mirror  ##查看使用bioconductor的默认镜像
options(BioC_mirror = "https://mirrors.ustc.edu.cn/bios")  ##指定镜像，这个时中国科技大学镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))  ##指定install.packages安装镜像，这个是清华镜像
options()$repos
options()$BioC_mirror


