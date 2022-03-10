



res1path <- paste(result_path,"/result1_GCMS",sep = "")
dir.create(res1path)


psG  =readRDS("./a2_dat_proj/GC_MS/ps_normIS.rds")
map = sample_data(psG)
map$Group = c(rep("Health",4),rep("Disease",4))
sample_data(psG) = map
rank_names(psG)


tax = as.data.frame(vegan_tax(psG))
table(tax$classcify)

#--------result1  soil GC analyses plot#------------
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")


barpath = paste(res1path,"/GC_soil_tax_class_barplot/",sep = "")
dir.create(barpath)

tax = as.data.frame(vegan_tax(psG))
table(tax$classcify)


ps_rela  = transform_sample_counts(psG, function(x) x / sum(x) )


library(EasyMicrobiome)

psdata <- tax_glom_wt(ps = ps_rela ,ranks = "classcify")

otu = as.data.frame((vegan_otu(psdata)))
map = as.data.frame(sample_data(psdata))
map
mapp = map[,c(2,4)]
colnames(mapp) = c("ID","group")
data = cbind(mapp,otu)
result = MuiKwWlx(data = data,num = c(3:ncol(data)))

data = data
i =c(3:ncol(data)) 
result = result

res <- MuiPlotStackBar(data = data,i =c(3:ncol(data)) ,result = result,errbar = F)
#--提取图片
g1_1 = res[[1]] + scale_fill_manual(values = colset3) +  
  # scale_x_discrete(limits = axis_order) +
  mytheme1 + labs (y = "Relative abundance (%)")

g1_1
#--提取数据
plotdata = res[[2]]
j = "GCsoil"

FileName1 <- paste(barpath,"/a2_",j,"_wlx_sbar",".pdf", sep = "")
ggsave(FileName1, g1_1, width = 8, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_wlx_sbar",".jpg", sep = "")
ggsave(FileName2, g1_1, width = 8, height =8 )

FileName2 <- paste(barpath,"/a2_",j,"_wlx_sbar_data",".csv", sep = "")
write.csv( plotdata,FileName2 )

result1 = FacetMuiPlotresultBox(data = data,num = c(3:ncol(data)),result = result,sig_show ="abc",ncol = 4 )
g1_2 <- result1[[1]] + scale_fill_manual(values = colset1) +  
  # scale_x_discrete(limits = axis_order) + 
  mytheme1
g1_2


FileName1 <- paste(barpath,"/a2_",j,"_wlx_sbar_facet",".pdf", sep = "")
ggsave(FileName1, g1_2, width = 15, height =10 )
FileName2 <- paste(barpath,"/a2_",j,"_wlx_sbar_facet",".jpg", sep = "")
ggsave(FileName2, g1_2, width = 15, height =10)


#--------beta ordanation plot #--------
betapath = paste(res1path,"/GCbeta/",sep = "")
dir.create(betapath)

source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

result = BetaDiv(ps = psG, group = "Group", dist = "bray", method = "PCA", Micromet = "adonis", pvalue.cutoff = 0.05)

#不带标签默认出图
# ,guide = NULL,guide = NULL
g2_1 = result[[1]]  + scale_fill_manual(values = colset1)  + mytheme1 +
  theme(legend.position = c(0.6,0.2))
g2_1
#带标签图形出图
g2_2 = result[[3]] + scale_fill_manual(values = colset1) + 
  # scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  theme(legend.position = c(0.6,0.2))
g2_2


#--提取出图数据
plotdata =result[[2]]
head(plotdata)
#--提取两两比较结果
pairResult =result[[4]]
head(pairResult)

#提取总体比较
TResult =result[[5]]
head(TResult)



FileName <- paste(betapath,"/PCA_plotdata.csv", sep = "")
write.csv(plotdata,FileName)

FileName <- paste(betapath,"/pairResult.csv", sep = "")
write.csv(pairResult,FileName)

FileName <- paste(betapath,"/TResult.csv", sep = "")
write.csv(TResult,FileName)

FileName <- paste(betapath,"/a2_PCA.pdf", sep = "")
ggsave(FileName, g2_1, width = 6, height = 6)


FileName1 <- paste(betapath,"/a2_PCA_label.pdf", sep = "")
ggsave(FileName1 , g2_2, width = 6, height = 6)


plotdata =result[[2]]
head(plotdata)


# 求均值
cent <- aggregate(cbind(x,y) ~ Group, data = plotdata, FUN = mean)
cent
# 合并到样本坐标数据中
segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
              by = 'Group', sort = FALSE)


g2_1$layers[[2]] = NULL
# library(ggcor)
library(ggsci)

g2_3 = g2_1 +geom_segment(data = segs,
                          mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group)) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 20,color = "blue",fill = "blue") +
  scale_fill_manual(values = colset1)+
  scale_color_manual(values = colset1,guide = F) +
  mytheme1 + theme(legend.position = c(0.4,0.2))
g2_3



# 提取两两检测结果
pair = result[[4]]
FileName <- paste(betapath,"Pair_adonis.csv", sep = "")
write.csv(pair,FileName)


FileName <- paste(betapath,"/a2_bray_PCA_adjust.pdf", sep = "")
ggsave(FileName, g2_3, width = 6, height = 6)
FileName1 <- paste(betapath,"/a2_bray_PCA_adjust.jpg", sep = "")
ggsave(FileName1 , g2_3, width = 6, height = 6)

#----------GC different soil heatmap #------------


#-----人工指定分组信息
group1 = c("Disease","Health")
# group2 = c("CK","CF")
b= data.frame(group1)

source("G:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")

# path = getwd()
wlxpath = paste(res1path,"/GCwlxtest/",sep = "")
dir.create(wlxpath)
ps_rela  = transform_sample_counts(psG, function(x) x / sum(x) )
result = statSuper(ps = ps_rela,group  = "Group",artGroup = b,method = "wilcox")
head(result)



fileName =  paste(wlxpath,"/wlxtest_arti.csv",sep = "")
write.csv(result,fileName)

#----筛选差异大的分泌物可视化
#----筛选差异大的分泌物可视化
head(result)
difres <- result %>%  filter(Disease_Health_fdr < 0.05 )
dim(difres)

fileName=  paste(wlxpath,"/wlxtest_arti_fdr05.csv",sep = "")
write.csv(difres,fileName)

#---差异分泌物可视化--这里的是标准化到100的丰度，注意不是到1.
#---
data <- difres[,c(sample_names(psG))]
head(data)
data$ID = paste(difres$id,difres$Metabolite.name,sep = "_")

row.names(data) = data$ID
data$ID = NULL

# data <-sqrt (data)
data[data > 0.4]<-0.4
# wt2<-sqrt(wt2)

#,fontsize_col = 10 修改行大小
#fontsize_row =  10 修改列自己大小
library(pheatmap)
map = as.data.frame(sample_data(psG))
map
# annotation_row = data.frame(map$Group)
# rownames(annotation_row) = row.names(map)



annotation_col = data.frame(map$Group)
rownames(annotation_col) =  row.names(map)



g3_1 = pheatmap(data,fontsize=6,cellwidth = 10, cellheight = 6,cluster_rows = TRUE,
                color = colorRampPalette(brewer.pal(11,"Spectral"))(60),
                display_numbers = F,fontsize_col = 10,fontsize_row =  10,
                annotation_col = annotation_col)
g3_1

fileName =  paste(wlxpath,"/pheartmap_diff.pdf",sep = "")
ggsave(fileName, g3_1, width = 12, height =20)

#--气泡图

library(reshape2)
data <- difres[,c(sample_names(psG))]
head(data)
data$ID = paste(difres$id,difres$Metabolite.name,sep = "_")

row.names(data) = data$ID
data$ID = NULL
head(data)
data$id = row.names(data)
pcm = melt(data, id = c("id"))
head(pcm)

# pcm$ID <- factor(pcm$ID,levels=unique(pcm$ID))

colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
#----样本在y轴上
g3_2 = ggplot(pcm, aes(y = id, x = variable)) + 
  geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
  labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
  # scale_fill_manual(values = colours, guide = FALSE) + 
  scale_x_discrete(limits = rev(levels(pcm$variable))) + mytheme1 

g3_2


fileName =  paste(wlxpath,"/bubbleplot_diff.pdf",sep = "")
ggsave(fileName,g3_2, width = 12, height =20 )

#------- network of soil GC MS result#-------------
library(tidyverse)
library(ggrepel)


netpath = paste(res1path,"/GCnetwork/",sep = "")
dir.create(netpath)
library(ggClusterNet)
library(igraph)
library(ggpubr)
result =ggClusterNet:: network(ps = ps_rela,
                               N = 0,
                               r.threshold=0.6,
                               p.threshold=0.05,
                               label = T,
                               path = netpath ,
                               zipi = F,
                               fill = "classcify",
                               ncol = 2
)



# 全部样本的网络比对
g4_1 = result[[1]] + mytheme1
g4_1
# 全部样本网络参数比对
data = result[[2]]
# plotname1 = paste(netpath,"/network_all.jpg",sep = "")
# ggsave(plotname1, g4_1,width = 15,height = 5)
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, g4_1,width = 32,height = 15)

g4_2 = result[[3]] + mytheme1
g4_2
# 全部样本网络参数比对
data = result[[2]]
# plotname1 = paste(netpath,"/network_all2.jpg",sep = "")
# ggsave(plotname1, g4_2,width = 15,height = 5)
plotname1 = paste(netpath,"/network_all2.pdf",sep = "")
ggsave(plotname1, g4_2,width = 32,height = 15)


tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)


#-------目标代谢物筛选
id <- c("melibiose 1",
"sucrose",
"ribose",
"lactic acid",
"xylose 1",
"myo-inositol",
"mannose 1",
"fructose 1",
"maltose",
"Gluconic lactone 1",
"ribitol")

head(difres)

difres$tax1[match(id,difres$tax1)]

tax = as.data.frame(vegan_tax(psG))

setdiff(id,tax$tax1)

intersect(id,tax$tax1)


# 差异代谢物的分析
difpath = paste(res1path,"/GC_diff_10/",sep = "")
dir.create(difpath)

data = read.csv("./a2_dat_proj/GC_MS/Diff_11_table/Diff_11_table.csv")
head(data)

dat = data[,c(2:9)]
head(dat)
row.names(dat) = data$tax
unique(data$tax)

dat = as.data.frame(t(dat))

dat$ID = row.names(dat)
dat$group = c(rep("Health",4),rep("Disease",4))

dat = dat %>% select(ID,group,everything())

#
result = MuiKwWlx(data = dat,num = c(3:13))
result

FileName <- paste(difpath,"/diff_10.csv", sep = "")
write.csv(result,FileName,sep = "")
FileName <- paste(difpath,"/diff_10.csv", sep = "")
write.csv(index,FileName,sep = "")

result1 = FacetMuiPlotresultBox(data = dat,num = c(3:13),result = result,sig_show ="abc",ncol = 3 )
p1_1 = result1[[1]] + 
  # scale_x_discrete(limits = axis_order) + 
  theme_bw() + mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1


res = FacetMuiPlotresultBar(data = dat,num = c(3:13),result = result,sig_show ="abc",ncol = 3)
p1_2 = res[[1]]+ 
  # scale_x_discrete(limits = axis_order) + 
  guides(color = FALSE) + theme_bw()  + 
  mytheme1+ guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = FacetMuiPlotReBoxBar(data = dat,num = c(3:13),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]+ 
  # scale_x_discrete(limits = axis_order) +
  theme_bw()  + 
  mytheme1 + guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(difpath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 9, height = 7)

FileName <- paste(difpath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = 9, height = 7)

FileName <- paste(difpath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = 9, height = 7)
