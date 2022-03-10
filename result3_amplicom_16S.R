


#---result1 base_diversity analyses-------------
#--creat path save result1 base_diversity
res1path <- paste(result_path,"/result3_amplicon16S",sep = "")
dir.create(res1path)

#--数据导入
ps = readRDS("./a2_dat_proj/缺-多物质-RS-番茄生长-定量16s-RS-发病率-测序16S/16S/result/ps构建/ps.rds")
map = sample_data(ps)
map$Group
# write.csv(map ,"./a2_dat_proj/缺-多物质-RS-番茄生长-定量16s-RS-发病率-测序16S/16S/mapNew.csv")
# map  = read.csv("./a2_dat_proj/缺-多物质-RS-番茄生长-定量16s-RS-发病率-测序16S/16S/mapNew.csv",row.names =1 )
# unique(map$Group)
# sample_data(ps) = map
# saveRDS(ps,"./a2_dat_proj/缺-多物质-RS-番茄生长-定量16s-RS-发病率-测序16S/16S/result/ps构建/ps.rds")
#--有些分析需要重复数量相同
ps <- subset_samples(ps,!Group%in% c("GP"))

#  order of x axis #------------
axis_order = c("PP","WP","PW","WW")
unique(map$Group)

#--1--基础多样性分析#——------

#--alpha多样性#---------
alppath = paste(res1path,"/alpha/",sep = "")
dir.create(alppath)

#---多种指标alpha多样性分析加出图-标记显著性
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")
index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
for (name_i in index) {
  # 计算alpha多样性和出图
  alp = alpha(ps = ps,group = "Group",inde=name_i,Plot = TRUE )
  p = alp[[1]]
  ## 作图数据是否提取
  plotdata = alp[[2]]
  FileName <- paste(alppath,name_i,"alpha_diversity.csv", sep = "")
  write.csv(plotdata,FileName,sep = "")
  
  FileName <- paste(alppath,name_i,"aov_bar", ".pdf", sep = "")
  # library("Cairo")
  ggsave(FileName, p, width = 12, height = 8)
  FileName1 <- paste(alppath,name_i,"_aov_bar", ".jpg", sep = "")
  ggsave(FileName1, p, width = 12, height = 8)
  if (name_i == "Shannon") {
    result = alp[[3]]
    FileName <- paste(alppath,"DATA_Alpha_diversity.csv", sep = "")
    write.csv(result,FileName)
  }
  
}

#--多种组合alpha分析和差异分析出图

# 使用分面出图，得益于我开发的R包EasyStat
alp = alpha(ps = ps,inde="Shannon",group = "Group",Plot = TRUE )
index= alp[[3]]
head(index)
#--从这里发现bof第六个样本存在问题：BOF_insect_6
sel = c(match("Shannon",colnames(index)),match("Richness",colnames(index)),match("Pielou_evenness",colnames(index)))


data = cbind(data.frame(ID = rep(1:length(index$Group)),group = index$Group),index[sel])
head(data)
#
result = MuiKwWlx(data = data,num = c(3:5))
result

FileName <- paste(alppath,"/alpha_diversity_different_label.csv", sep = "")
write.csv(result,FileName,sep = "")
FileName <- paste(alppath,"/alpha_diversity_7_index.csv", sep = "")
write.csv(index,FileName,sep = "")

result1 = FacetMuiPlotresultBox(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]] + 
  scale_x_discrete(limits = axis_order) + 
  theme_bw() + mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1


res = FacetMuiPlotresultBar(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1)
p1_2 = res[[1]]+ 
  scale_x_discrete(limits = axis_order) + 
  guides(color = FALSE) + theme_bw()  + 
  mytheme1+ guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = FacetMuiPlotReBoxBar(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1)
p1_3 = res[[1]]+ 
  scale_x_discrete(limits = axis_order) +
  theme_bw()  + 
  mytheme1 + guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(alppath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 5, height =9)

FileName <- paste(alppath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = 5, height = 9)

FileName <- paste(alppath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = 5, height = 9)


# 
# 
# #--alpha稀释曲线绘制#------
# 
# source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\alpha_rare_all.R",encoding = "utf-8")
# result = alpha_rare_all(ps = ps, group = "Group", method = "chao1", start = 200, step = 200)
# 
# #--提供单个样本溪稀释曲线的绘制
# p2_1 <- result[[1]] + mytheme1 +
#   scale_color_manual(values = colset1)
# ## 提供数据表格，方便输出
# raretab <- result[[2]]# data table
# head(raretab)
# #--按照分组展示稀释曲线
# p2_2 <- result[[3]]# output group curve plot
# #--按照分组绘制标准差稀释曲线
# p2_3 <- result[[4]]# output group curve with CI plot
# 
# FileName <- paste(alppath,"Alpha_rare_sample", ".pdf", sep = "")
# ggsave(FileName, p2_1, width = 8, height =6)
# FileName <- paste(alppath,"Alpha_rare_group", ".pdf", sep = "")
# ggsave(FileName, p2_2, width = 8, height =6)
# FileName <- paste(alppath,"Alpha_rare_groupwithSD", ".pdf", sep = "")
# ggsave(FileName, p2_3, width = 8, height =6)
# 
# FileName <- paste(alppath,"/Alpha_rare_data.csv", sep = "")
# write.csv(raretab,FileName,sep = "")


#---beta-diversity#-------------
betapath = paste(res1path,"/beta/",sep = "")
dir.create(betapath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")


# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA

methodlist = c("NMDS","PCoA", "PCA")

for (method in methodlist) {
  result = BetaDiv(ps = ps, group = "Group", dist = "bray", method = method, Micromet = "adonis", pvalue.cutoff = 0.05)
  p3_1 = result[[1]] + scale_fill_manual(values = colset1,guide = F)+
    scale_color_manual(values = colset1,guide = F) + mytheme1 + theme(legend.position = c(0.2,0.2))
  p3_1
  #带标签图形出图
  p3_2 = result[[3]] + scale_fill_manual(values = colset1,guide = F)+
    scale_color_manual(values = colset1,guide = F) + mytheme1 + theme(legend.position = c(0.2,0.2))
  p3_2
  
  FileName <- paste(betapath,"/a2_",method,"bray.pdf", sep = "")
  ggsave(FileName, p3_1, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"",method,"bray.jpg", sep = "")
  ggsave(FileName1 , p3_1, width = 12, height = 12)
  
  FileName <- paste(betapath,"/a2_",method,"bray_label.pdf", sep = "")
  ggsave(FileName, p3_2, width = 12, height = 12)
  FileName1 <- paste(betapath,"/a2_",method,"bray_label.jpg", sep = "")
  ggsave(FileName1 , p3_2, width = 12, height = 12)
  
  # 提取出图数据
  plotdata = result[[2]]
  FileName <-  paste(betapath,"/a2_",method,"bray.csv", sep = "")
  write.csv(plotdata,FileName)
  #---------排序-精修图
  plotdata =result[[2]]
  head(plotdata)
  # 求均值
  cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
  cent
  # 合并到样本坐标数据中
  segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
                by = 'Group', sort = FALSE)
  
  # p2$layers[[2]] = NULL
  # library(ggcor)
  library(ggsci)
  p3_3 = p3_1 +geom_segment(data = segs,
                            mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
    geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = FALSE) + mytheme1 + theme(legend.position = c(0.2,0.2))
  p3_3
  
  FileName <- paste(betapath,"/a2_",method,"bray_star.pdf", sep = "")
  ggsave(FileName, p3_3, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"bray_star.jpg", sep = "")
  ggsave(FileName1 , p3_3, width = 8, height = 8)
  
}

#提取总体比较
TResult =result[[5]]
head(TResult)

# 提取两两检测结果
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_adonis.csv", sep = "")
write.csv(pair,FileName)
FileName <- paste(betapath,"Total_adonis.csv", sep = "")
write.csv(TResult,FileName)


#---微生物组成分析#——-------


#-------门类水平展示
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")
barpath = paste(res1path,"/Microbial_composition/",sep = "")
dir.create(barpath)
j = "Phylum"
result = barMainplot(ps = ps,j = "Phylum",axis_ord = NULL,label = FALSE ,sd = FALSE,Top = 10)
p4_1 <- result[[1]] + scale_fill_brewer(palette = "Paired") + 
  scale_x_discrete(limits = axis_order) + 
  mytheme1
p4_1
p4_2  <- result[[3]] +
  scale_fill_brewer(palette = "Paired") + 
  scale_x_discrete(limits = axis_order) + mytheme1
p4_2

databar <- result[[2]] 

FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
ggsave(FileName1, p4_2, width = 8, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
ggsave(FileName2, p4_2, width = 8, height =8 )

FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
ggsave(FileName1, p4_1, width = 8, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
ggsave(FileName2, p4_1, width = 8, height =8 )

FileName <- paste(barpath,"/a2_",j,"_bar_data",".csv", sep = "")
write.csv(databar,FileName)


#--距离和丰度合并

source("G:\\Shared_Folder\\Function_local\\R_function\\micro/cluMicro.bar.R")
library("ggdendro")
library(phyloseq)
library(tidyverse)
library(ggtree)
library( ggstance)


result <-  cluMicro.bar (dist = "bray",
                         ps = ps,
                         j = "Phylum",
                         Top = 10, # 提取丰度前十的物种注释
                         tran = TRUE, # 转化为相对丰度值
                         hcluter_method = "complete",
                         Group = "Group",
                         cuttree = 2
)


p5_1 <- result[[1]]
p5_1

p5_2 <- result[[2]]
p5_2
clubardata <- result[[3]]

FileName1 <- paste(barpath,"/a2_",j,"_cluster",".pdf", sep = "")
ggsave(FileName1, p5_1, width = 12, height =8 )

FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar",".pdf", sep = "")
ggsave(FileName1, p5_2, width = 12, height =8 )

FileName <- paste(barpath,"/a2_",j,"_cluster_bar_data",".csv", sep = "")
write.csv(clubardata,FileName)


#---共有微生物特有微生物
#---flower plot#-------
flowpath = paste(res1path,"/flowplot/",sep = "")
dir.create(flowpath)





# source("G:\\Shared_Folder\\Function_local\\R_function\\micro/flowerplot.R")
# flowerplot(ps = ps,rep = 6,path =flowpath )
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/ggflowerplot.R")
p0_1 <- ggflower(ps = ps,
                 # rep = 1,
                 group = "Group",
                 start = 1, # 风车效果
                 m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                 a = 0.2, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "yellow"
)

# p + scale_fill_brewer(palette = "Paired")
FileName1 <- paste(flowpath,"ggflowerID.pdf", sep = "")
ggsave(FileName1, p0_1, width = 8, height = 8)
FileName2 <- paste(flowpath,"ggflowerID.jpg", sep = "")
ggsave(FileName2, p0_1, width = 8, height = 8 )

#---Ven-Upset#----------
Venpath = paste(res1path,"/Ven_Upset/",sep = "")
dir.create(Venpath)
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/Ven-Upset.R")
# otutab = as.data.frame(otu_table(ps_sub))
# map = as.data.frame(sample_data(ps_sub))


result = VenUpset(ps = ps,
                  group = "Group",
                  # rep = 14,
                  path = Venpath
                  
)
grid.draw(T)

filename3 <- paste(Venpath,"Upset.pdf", sep = "")
pdf(file=filename3,width = 12, height = 12)
## upset 我不会直接输出一个图形对象，我现在先在这里把数据提出来在出图
upset(result[[2]], sets = colnames(result[[2]]),
      number.angles = 30, point.size = 2, line.size = 1,
      mainbar.y.label = "OTU", sets.x.label = "OTU Per Treatment",
      text.scale = c(2, 2, 2,2, 2, 2),mb.ratio = c(0.7, 0.3),order.by = "freq",keep.order = TRUE,
      queries = list(list(query = intersects, params =
                            list(colnames(result[[2]])), color = "red", active = T),
                     list(query = intersects, params =
                            list(colnames(result[[2]])), color = "red", active = T),
                     list(query = intersects, params =
                            list(colnames(result[[2]])), color = "red", active = T)))

dev.off()

#--二分网络#-------
biospath = paste(res1path,"/biospr_network_Ven/",sep = "")
dir.create(biospath)

ps_sub1 = filter_taxa(ps, function(x) sum(x ) > 800 , TRUE)
N = 0.5
result = div_network(ps_sub1)
edge = result[[1]]
head(edge)
data = result[[3]]


result <- div_culculate(table = result[[3]],distance = 1.1,distance2 = 1.5,distance3 = 1.3,order = FALSE)
# result <- div_culculate(table = result[[3]],distance = 1,distance2 = 1.2,distance3 = 1.1,order = FALSE)
edge = result[[1]]
head(edge)
plotdata = result[[2]]
head(plotdata)
#--这部分数据是样本点数据
groupdata <- result[[3]]
# table(plotdata$elements)
node =  plotdata[plotdata$elements == unique(plotdata$elements), ]

otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
res = merge(node,tax_table,by = "row.names",all = F)
dim(res)
head(res)
row.names(res) = res$Row.names
res$Row.names = NULL
plotcord = res

xx = data.frame(mean  =rowMeans(otu_table))
head(xx)
plotcord = merge(plotcord,xx,by = "row.names",all = FALSE)
head(plotcord)
# plotcord$Phylum
row.names(plotcord) = plotcord$Row.names
plotcord$Row.names = NULL
head(plotcord)
library(ggrepel)
p = ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
                            data = edge, size = 0.3,color = "yellow") +
  geom_point(aes(X1, X2,fill = Phylum,size =mean ),pch = 21, data = plotcord) +
  geom_point(aes(X1, X2),pch = 21, data = groupdata,size = 5,fill = "blue",color = "black") +
  geom_text_repel(aes(X1, X2,label = elements ), data = groupdata) +
  theme_void()

p

filename = paste(biospath,"/","biostr_Ven_network.pdf",sep = "")
ggsave(filename,p,width = 8,height = 6)
filename = paste(biospath,"/","biostr_Ven_network.png",sep = "")
ggsave(filename,p,width = 8,height = 6)


#--差异分析edger#----
diffpath = paste(res1path,"/diff_tax/",sep = "")
dir.create(diffpath)
# 准备脚本
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\Plot.CompareWithCK.R",encoding = "utf-8")

unique(map$Group)
# group1 = c("GP","WW")
# group2 = c("WP","WW")
# group3 = c("PW","WW")
# group4 = c("PP","WW")
# b= data.frame(group1,group2,group3,group4)
# ,artGroup = b
res = EdgerSuper(ps = ps,group  = "Group")
head(res)

filename = paste(diffpath,"/","edger.csv",sep = "")
write.csv(res,filename)


# p <- Plot.CompareWithCK(ps = ps,CK = "WW",j = "Genus",abun = 0.001,result = res )
# p <- p +
#   scale_fill_brewer(palette = "Paired") + 
#   scale_x_discrete(limits = axis_order) + mytheme1
# p
# 
# filename = paste(diffpath,"/","edger_001_diff_bio_plot.pdf",sep = "")
# ggsave(filename,p,width = 14,height = 8)

#----差异热图----------
#----差异气泡图--------

heatpath = paste(res1path,"/diff_heapmap_boplot/",sep = "")
dir.create(heatpath)

source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")

#提取丰度最高的前20个OTU做展示
ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps_rela
otu = as.data.frame(otu_table(ps_rela))
otu$mean <- rowMeans(as.data.frame(otu_table(ps_rela)))

idtab <- otu %>% arrange(desc(mean)) %>%
  head(n = 20)
id = row.names(idtab)


result <- Microheatmap(ps_rela = ps_rela,id = id)

p1 <- result[[1]]
p2 <- result[[2]]

filename = paste(heatpath,"/","Top20ggheatmap.pdf",sep = "")
ggsave(filename,p1,width = 14,height = 8)

filename = paste(heatpath,"/","Top20ggbubble.pdf",sep = "")
ggsave(filename,p2,width = 14,height = 8)

#----lefse-------------
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\lefse_py_pre.R",encoding = "utf-8")
lefpath = paste(res1path,"/lefse_py/",sep = "")
dir.create(lefpath)

library(phyloseq)
library(EasyMicrobiome)
library("tidyverse")

tablefse <- lefse_py_pre(ps = ps,taxGlomRank = "Genus",filter = 400)
head(tablefse)

filename = paste(lefpath,"/LEFSE_to_run_G_level.txt",sep = "")
write.table(tablefse,filename,append = F, quote = F,col.names= F,sep = "\t")

# #文件预处理
# format_input.py LEFSE_to_run_G_level.txt pri_lefse.in -c 1 -u 2 -o 1000000
# # 注意这里 –c用来指定分组信息-u 1指定样品信息
# #文件分析,这里-l设置LDA阈值，默认为2，我们使用4 会更加严格
# ~/src/nsegata-lefse/run_lefse.py pri_lefse.in pri_lefse_2.res  -l 2
# #柱状图绘制
# plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
# #树状图绘制
# plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf
# #做每个差异的柱状图
# mkdir biomarkers_raw_images
# plot_features.py pri_lefse.in pri_lefse_2.res biomarkers_raw_images/


#--机器学习#-------
matpath = paste(res1path,"/Machine_learing/",sep = "")
dir.create(matpath )
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\MicroMachine_learning.R")

library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
# #--三种机器学习方法评测
# result = MicroRoc( ps = ps,group  = "Group")
# #--提取roc曲线
# p <- result[[1]]
# p
# #提取AUC值
# data <- result[[2]]
# 
# filename = paste(matpath,"/three_method_AUCvalue.csv",sep = "")
# write.csv(data,filename,quote = F)
# 
# data <- result[[3]]
# filename = paste(matpath,"/three_method_AUCdata.csv",sep = "")
# write.csv(data,filename,quote = F)
# 
# filename = paste(matpath,"/three_method_AUC_plot.pdf",sep = "")
# ggsave(filename,p,width = 8,height = 8)

mapping = as.data.frame(sample_data(ps))

#--随机森林全套
result = MicroRF(ps = ps,group  = "Group",optimal = 40,rfcv = TRUE,nrfcvnum = 5,min = -1,max = 5)
#火柴图展示前二十个重要的OTU
p <- result[[1]]

filename = paste(matpath,"/randonforest_loading.pdf",sep = "")
ggsave(filename,p,width = 8,height = 12)
# 圈图展示
p <- result[[2]]
filename = paste(matpath,"/randonforest_loading_circle.pdf",sep = "")
ggsave(filename,p,width = 8,height = 12)
# 展示交叉验证结果
p <- result[[3]]
filename = paste(matpath,"/randonforest_cross_check.pdf",sep = "")
ggsave(filename,p,width = 8,height = 12)



data <- result[[5]]
filename = paste(matpath,"/randomforest_data.csv",sep = "")
write.csv(data,filename,quote = F)

data <- result[[4]]
filename = paste(matpath,"/randomforest_cross_data.csv",sep = "")
write.csv(data,filename,quote = F)



#-- network #-----------
netpath = paste(res1path,"/network/",sep = "")
dir.create(netpath)

library(igraph)
library(sna)

result = ggClusterNet::network(ps = ps,
                               N = 120,
                               r.threshold=0.8,
                               p.threshold=0.05,
                               label = FALSE,
                               path = netpath,
                               zipi = TRUE,
                               ncol = 5)

# 全部样本的网络比对
p4_1 = result[[1]] + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p4_1 + mytheme1
# 全部样本网络参数比对
data = result[[2]]
# plotname1 = paste(netpath,"/network_all.jpg",sep = "")
# ggsave(plotname1, p4_1,width = 64,height = 16,limitsize = FALSE)
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p4_1,width = 78,height = 16,limitsize = FALSE)

tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)

p4_1 = result[[3]] + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p4_1 + mytheme1
# 全部样本网络参数比对
data = result[[2]]
# plotname1 = paste(netpath,"/network_all.jpg",sep = "")
# ggsave(plotname1, p4_1,width = 64,height = 16,limitsize = FALSE)
plotname1 = paste(netpath,"/network_all2.pdf",sep = "")
ggsave(plotname1, p4_1,width = 78,height = 16,limitsize = FALSE)




#---对多物质验证发病率的研究进行深入挖掘
difpath = paste(res1path,"/RS_abundancee/",sep = "")
dir.create(difpath)

#--计算病原菌丰度
ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 

# tax = as.data.frame(vegan_tax(ps_rela))

rank_names(ps_rela)
pssub <- ps_rela %>%
  subset_taxa(
    Species %in% "Ralstonia_solanacearum"
    # row.names(tax_table(ps_rela ))%in%c("")
  )
pssub

otu = as.data.frame(vegan_otu(pssub))
head(otu)
otu$ID = row.names(otu)
dat <- as.tibble(sample_data(pssub)) %>% inner_join(otu) %>% as.data.frame()
colnames(dat)[2] <- "group"

result = MuiKwWlx(data = dat,num = c(3))
result

FileName <- paste(difpath,"/RS_diff.csv", sep = "")
write.csv(result,FileName,sep = "")


result1 = FacetMuiPlotresultBox(data = dat,num = c(3),result = result,sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]] + 
  # scale_x_discrete(limits = axis_order) + 
  theme_bw() + mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1


res = FacetMuiPlotresultBar(data = dat,num = c(3 ),result = result,sig_show ="abc",ncol = 1)
p1_2 = res[[1]]+ 
  # scale_x_discrete(limits = axis_order) + 
  guides(color = FALSE) + theme_bw()  + 
  mytheme1+ guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = FacetMuiPlotReBoxBar(data = dat,num = c(3 ),result = result,sig_show ="abc",ncol = 1)
p1_3 = res[[1]]+ 
  # scale_x_discrete(limits = axis_order) +
  theme_bw()  + 
  mytheme1 + guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(difpath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 8, height =5)

FileName <- paste(difpath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = 8, height = 5)

FileName <- paste(difpath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = 8, height = 5)

#四个大门和RS的相关性质
ps_rela =  subset_samples(ps_rela,Group %in% c("PW"))

ps_rela2 = filter_taxa(ps_rela, function(x) sum(x ) > 0.006, TRUE);ps_rela2 #筛选序列数量大于1的

pssub <- ps_rela %>%
  subset_taxa(
    # !Phylum %in% c("Proteobacteria") |
    Species %in% "Ralstonia_solanacearum"|
    row.names(tax_table(ps_rela ))%in% row.names(tax_table( pssub ))
  )




pssub <- pssub %>%
  subset_taxa(
    # Phylum %in% c("Firmicutes","Proteobacteria","Actinobacteria","Bacteroidetes") |
    !Phylum %in% c("Proteobacteria") 
    #   Species %in% "Ralstonia_solanacearum" |
    # row.names(tax_table(ps_rela ))%in% row.names(tax_table(ps_rela2 ))
  )

pssub

tax = as.data.frame(vegan_tax(pssub))
head(tax)
A = tax$Phylum
A[tax$Species %in% "Ralstonia_solanacearum"] = "Ralstonia_solanacearum"
unique(A)
tax$filed = A
tax_table(pssub) = as.matrix(tax)

group2 <- data.frame(SampleID = row.names(tax),Group = A)
group2$taxg = gsub("g__","",tax$Family)
colnames(group2) <- c("ID","group","taxg")
group2$group  =as.factor(group2$group)

result <- corBiostripe(ps = pssub,
                       # group = group2,
                       r.threshold = 0.9,
                       p.threshold = 0.05
                       )

#--提取相关矩阵
cor = result[[1]]

result2 = PolygonRrClusterG (cor = cor,nodeGroup = group2,zoom = 1,zoom2 = 1 )
node = result2[[1]]

### nodeadd 节点注释的简单封装，便捷实用otu表格和分组文件进行注释
library(dplyr)
nodes <- node %>%
  inner_join(group2,by =c("elements" = "ID") )
#-----计算边#--------
edge = edgeBuild(cor = cor,plotcord = node)
head(edge)
### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group),pch = 21, size = 5,data = nodes) +
  scale_colour_brewer(palette = "Set1") + theme_void()
pnet





a = as.data.frame(tax_table(pssub))
unique(a$Phylum)
#--网络

data = as.data.frame(vegan_otu(pssub))
head(data)



occor = psych::corr.test(data,use="pairwise",method= method,adjust="fdr",alpha=.05)
print("over")
occor.r = occor$r #
occor.p = occor$p #

r.threshold  = 0
p.threshold  = 0

# occor.r[occor.p > p.threshold&abs(occor.r)< r.threshold] = 0


A <- levels(as.factor(group2$Group))
A
Gru = group2
i = 1
for (i in 1:length(A)) {
  fil <- intersect(row.names(occor.r),as.character(Gru$SampleID[Gru$Group == A[i]]))
  a <- row.names(occor.r) %in% fil
  occor.r[a,a] = 0
  occor.p[a,a] = 1
}

length(as.vector(occor.r["ASV_198",]))

as.vector(occor.r["ASV_198",])[as.vector(occor.r["ASV_198",]) < 0] %>% length()
