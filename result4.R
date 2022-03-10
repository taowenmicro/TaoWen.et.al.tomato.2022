res1path <- paste(result_path,"/result4",sep = "")
dir.create(res1path)

difpath = paste(res1path,"/8_RS_number_day135/",sep = "")
dir.create(difpath)

dat = read.csv("./a2_dat_proj/八株菌-多物质-无菌番茄-培养/RS和八株菌涂布结果/RS和八株菌涂布结果.csv")
head(dat)


result = MuiKwWlx(data = dat,num = c(3:5))
result

FileName <- paste(difpath,"/RS_diff.csv", sep = "")
write.csv(result,FileName,sep = "")


result1 = FacetMuiPlotresultBox(data = dat,num = c(3:5),result = result,sig_show ="abc",ncol = 3 )
p1_1 = result1[[1]] + 
  # scale_x_discrete(limits = axis_order) + 
  theme_bw() + mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1


FileName <- paste(difpath,"RS", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 10, height =3)

result = MuiKwWlx(data = dat,num = c(6:8))
result

FileName <- paste(difpath,"/RS_diff.csv", sep = "")
write.csv(result,FileName,sep = "")


result1 = FacetMuiPlotresultBox(data = dat,num = c(6:8),result = result,sig_show ="abc",ncol = 3 )
p1_2 = result1[[1]] + 
  # scale_x_discrete(limits = axis_order) + 
  theme_bw() + mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_2
FileName <- paste(difpath,"bac_8", ".pdf", sep = "")
ggsave(FileName, p1_2, width = 10, height =3)

library(patchwork)

p <- p1_1/p1_2

FileName <- paste(difpath,"bac_RS", ".pdf", sep = "")
ggsave(FileName, p, width = 10, height =6)

#--Rs 和八株菌的相关

head(dat)
dat2  =dat[,c(2,3,6)]


head(dat2)
colnames(dat2) = c("group","RS","SynCom")
dat3 <- dat2 %>% filter(group %in%c("PP100","PP10","PP1"))
library(ggpubr)

p = ggplot2::ggplot(dat3,aes(x = RS,y = SynCom)) + 
  ggplot2::geom_point(aes(color = group)) + 
  ggpubr::stat_cor(label.y=lab*0.01)+
  ggpubr::stat_regline_equation(label.y=lab*5) +
  geom_smooth( method=lm, se=T) + 
  scale_color_manual(values = colset1) + 
  mytheme1 +  theme(legend.position = c(0.2,0.2))
p 
library("ggExtra")
p1 = ggMarginal(p, groupColour = TRUE, groupFill = TRUE)  
p1


FileName <- paste(difpath,"cor_PP", ".pdf", sep = "")
ggsave(FileName, p1, width = 6, height =6)

#--葡萄糖加菌
dat2$group
dat3 <- dat2 %>% filter(group %in%c("GP10","GP1","GP100"))
library(ggpubr)

p = ggplot2::ggplot(dat3,aes(x = RS,y = SynCom)) + 
  ggplot2::geom_point(aes(color = group)) + 
  ggpubr::stat_cor(label.y=lab*0.01)+
  ggpubr::stat_regline_equation(label.y=lab*5) +
  geom_smooth( method=lm, se=T) + 
  scale_color_manual(values = colset1) + 
  mytheme1 +  theme(legend.position = c(0.2,0.2))
p 
library("ggExtra")
p2 = ggMarginal(p, groupColour = TRUE, groupFill = TRUE)  
p2


FileName <- paste(difpath,"cor_GP", ".pdf", sep = "")
ggsave(FileName, p2, width = 6, height =6)



#--水
dat2$group
dat3 <- dat2 %>% filter(group %in%c("WW"))
library(ggpubr)
  
p3 = ggplot2::ggplot(dat3,aes(x = RS,y = SynCom)) + 
  ggplot2::geom_point(aes(color = group)) + 
  ggpubr::stat_cor(label.y=lab*0.01)+
  ggpubr::stat_regline_equation(label.y=lab*5) +
  geom_smooth( method=lm, se=T) + 
  scale_color_manual(values = colset1) + 
  mytheme1 +  theme(legend.position = c(0.2,0.2))
 
library("ggExtra")
p3 = ggMarginal(p3, groupColour = TRUE, groupFill = TRUE)  
p3


FileName <- paste(difpath,"cor_WW", ".pdf", sep = "")
ggsave(FileName, p3, width = 6, height =6)


