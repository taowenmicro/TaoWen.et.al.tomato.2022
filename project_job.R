


# 备忘： 数据是一个有机的组合，代码也是一个有机的组合，proj是一个总览
pap_path <- "a0_paper_proj"
dir.create(pap_path)

cod_path <- "a1_coding_proj"
dir.create(cod_path)

dat_path <- "a2_dat_proj"
dir.create(dat_path)

result_path <- "a2_result_proj"
dir.create(result_path)


#---library R  package-------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyStat)
library(EasyMicrobiome)
library(ggthemes)
library(RColorBrewer)#调色板调用包

library(igraph)
library(sna)
# 设置主题
mytheme1 = theme_bw() + theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  legend.position="right",
  
  legend.title = element_blank(),
  legend.background=element_blank(),
  legend.key=element_blank(),
  # legend.text= element_text(size=7),
  # text=element_text(),
  # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
  plot.title = element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
  axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
  axis.text = element_text(size = 10,face = "bold"),
  axis.text.x = element_text(colour = "black",size = 8),
  axis.text.y = element_text(colour = "black",size = 8),
  legend.text = element_text(size = 15,face = "bold")
)


#----set color or fill#------------

#调用所有这个包中的调色板
display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
display.brewer.pal(9,"Dark2")
# colset1 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))
colset1 <- c(brewer.pal(12,"Dark2"))
colset2 <- brewer.pal(12,"Dark2")
colset3 <- c(brewer.pal(12,"Dark2"),brewer.pal(9,"Pastel1"))


data(iris)
head(iris)
p <- ggplot(iris) + geom_point(aes(x = Sepal.Length,y = Sepal.Width ,color = Species ))
p + mytheme1 + scale_color_manual(values = colset1)

