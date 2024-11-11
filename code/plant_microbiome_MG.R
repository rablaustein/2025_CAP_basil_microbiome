######################################
### Plant Microbiome Analysis ###
######################################

# Mairui Gao
# mrgao@umd.edu
# last edited: 11/08/2024

######################################

# Load library
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(tidyverse)
library(vegan)
library(SCRuB)
library(stats)

### set working directory
#setwd("<path>/CAP_basil_microbiome")

###
### (1) Process taxa tables
###

### load files
### SCRuB read counts from negative control DNA extractions and library preps
###

## metadata file for scrub (contaminant removal)
meta_scrub = read.delim("metadata/R_meta_scrub.txt", h=T, row.names=1)
meta_scrub = meta_scrub[order(rownames(meta_scrub)),]
rownames(meta_scrub) = str_replace_all(rownames(meta_scrub), "-", ".")
head(meta_scrub)

## phylum
# set table
p<- read.delim("taxa_tables/p.txt", h=T, row.names=1, check.names = FALSE)
# remove singletons
p <- p[which(apply(p,1,sum) > 0),]
colnames(p) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(p), meta_scrub,
                 control_order = c("control_DNA_extraction"))
taxa_phy_s = as.data.frame(t(scr_out$decontaminated_samples))

## family
# set table
f <- read.delim("taxa_tables/f.txt", h=T, row.names=1, check.names = FALSE)
# remove singletons
f <- f[which(apply(f,1,sum) > 0),]
colnames(f) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(f), meta_scrub,
                 control_order = c("control_DNA_extraction"))
taxa_fam_s <- as.data.frame(t(scr_out$decontaminated_samples))

## genus
# set table
g <- read.delim("taxa_tables/g.txt", h=T, row.names=1, check.names = FALSE)
# remove singletons
g <- g[which(apply(g,1,sum) > 0),]
colnames(g) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(g), meta_scrub,
                 control_order = c("control_DNA_extraction"))
taxa_gen_s <- as.data.frame(t(scr_out$decontaminated_samples))

## ASV
# set table
asv <- read.delim("taxa_tables/feature-table.txt", h=T, row.names=1, check.names=FALSE)
asv <- asv[which(apply(asv,1,sum) > 0),]
colnames(asv) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(asv), meta_scrub,
                 control_order = c("control_DNA_extraction"))
taxa_asv_s <- as.data.frame(t(scr_out$decontaminated_samples))

### remove singletons
###
taxa_phy_s1 = taxa_phy_s[which(apply(taxa_phy_s, 1, sum) > 1),]
taxa_fam_s1 = taxa_fam_s[which(apply(taxa_fam_s, 1, sum) > 1),]
taxa_gen_s1 = taxa_gen_s[which(apply(taxa_gen_s, 1, sum) > 1),]
taxa_asv_s1 = taxa_gen_s[which(apply(taxa_asv_s, 1, sum) > 1),]

#write.csv(taxa_phy_s,'taxa_tables/p.scrub.csv')
#write.csv(taxa_fam_s,'taxa_tables/f.scrub.csv')
#write.csv(taxa_gen_s,'taxa_tables/g.scrub.csv')
#write.csv(taxa_asv_s,'taxa_tables/feature-table.scrub.csv')


###
### (2) alpha diversity
###

rm(list=ls())

df<-read.csv("taxa_tables/feature-table.scrub.csv", row.names = 1)
df_per<-sweep(df,2,colSums(df), '/')
#write.csv(df_per,"taxa_tables/feature-table.scrub.RA.csv")

df<-read.csv('taxa_tables/feature-table.scrub.RA.csv', header = T, row.names = 1)
Shannon<-diversity(df, index = "shannon", MARGIN = 2)
Shannon<-as.data.frame(Shannon)
Shannon<-data.frame(sample = rownames(Shannon), Shannon)

mygroup<-read.delim('metadata/R_meta.txt', header = T)
colnames(mygroup)[1]<-c('sample')

df1<-Shannon
df2<-merge(df1, mygroup, by='sample')
df2$Treatment<-factor(df2$Treatment,levels = unique(df2$Treatment))
mycol<-c("Control"='#CCCCCC',"50w"='#74C476',"90w"='#FDAE6B', "120w"='#FB9A99')

ggboxplot(df2, x = "Treatment", y = "Shannon",fill="Treatment",
          add = c('jitter'),  
          shape = 21,  
          bxp.errorbar = TRUE,  
          bxp.errorbar.width = 0.25,  
          add.params = list(size = 3, alpha = 0.5, width = 0.5),  
          outlier.shape = NA  # Omits plotting outliers
) + 
  scale_fill_manual(values = mycol) +  
  labs(y = "Shannon", x = "Treatment") +  
  theme_bw() +  
  theme(
    legend.title = element_blank(),  
    legend.position = 'none',  
    axis.title = element_text(face = 'bold', size = 16),  
    axis.text = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

#ggsave('alpha.diversity.pdf', width=6, height=4)

kruskal.test(Shannon~Treatment, data=df2)


###
### (3) beta diversity
###

rm(list=ls())
family<- read.csv('taxa_tables/f.scrub.RA.csv',header = T, row.names = 1,check.names = F)
family<-as.data.frame(t(family))
dist <- vegdist(family, method="bray")
dist <- as.matrix(dist)
nmds_result<-  metaMDS(dist,k=2)
stress <- paste0("Stress=", format(nmds_result$stress, scientific = TRUE, digits = 3))
nmds12 <- as.data.frame(nmds_result$points)
nmds12$samples<-rownames(nmds12)
groups <-read.delim('metadata/R_meta.txt', header = T)
colnames(groups)[1]<-'samples'
nmds12<-merge(nmds12,groups,by='samples')
nmds12$Distance <- factor(nmds12$Distance)
nmds12$Treatment<-factor(nmds12$Treatment,levels = c("Control", "50w", "90w","120w"))

mycol<-c('#CCCCCC','#74C476','#FDAE6B', '#FB9A99')

ggscatter(nmds12, x="MDS1", y="MDS2", 
          color = "Treatment", shape = "Distance", size = "ROS", palette = mycol,
          ellipse = FALSE, conf.int.level = 0.95) +
  labs(x = 'NMDS1', y = 'NMDS2',
       color = 'Treatment', 
       shape = 'Distance',  
       size = 'ROS') +       
  geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "solid") +
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "solid") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_size_continuous(
    range = c(2, 6),
    breaks = c(15, 30, 45, 60),
    labels = c("15%", "30%", "45%", "60%")
  ) +  
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.position = "right", 
        text = element_text(size=14),
        legend.title = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  guides(color = guide_legend(order = 1, title = "Treatment"), 
         shape = guide_legend(order = 2, title = "Distance"), 
         size = guide_legend(order = 3, title = "Normalized Fluorescence(%)"))

#ggsave('beta2.pdf', width=8.5, height=5)


###
### (4) Permanova
###

rm(list=ls())
family<- read.csv('taxa_tables/f.scrub.RA.csv',header = T, row.names = 1,check.names = F)
family<-as.data.frame(t(family))
family<-family[-c(grep("Ct_qPCR_pos_zymo", rownames(family))),]
dist <- vegdist(family, method="bray")
dist <- as.matrix(dist)
groups <-read.delim('metadata/R_meta.txt',header = T)
colnames(groups)[1]<-'samples'
adonis2(dist ~ groups$ROS, permutations = 999)

# abundance: filter at 10% prevalence to remove rare taxa
myinput<-"taxa_tables/f.scrub.csv"
present<-0.1 #10% prevalence
df<-read.csv(file=myinput,header=T,row.names=1,check.names=FALSE)
df<-na.omit(df)
sample<-colnames(df)[sapply(df, is.numeric)]
df1<-df[,sample]

#df1$not0<-apply(df1, MARGIN = 1, function(x) sum(x))
df1$not0<-apply(df1, MARGIN = 1, function(x) length(which(x>0)))
df2<-df[which(df1$not0>round(length(sample)*present)),]
#df3<-df[which(df1$not0<round(length(sample)*present)),]
#write.csv(df2, "taxa_tables/f.scrub.0.1.csv", row.names=TRUE)

# calculate relative abundance
myinput<-"taxa_tables/f.scrub.0.1.csv"
df<-read.csv(file=myinput, header=T, check.names=FALSE, row.names=1)
df_per<-sweep(df,2,colSums(df), '/')
#write.csv(df_per, "taxa_tables/f.scrub.0.1.RA.csv", row.names=TRUE)

# compute mean value 
rm(list=ls())
myinput<-"taxa_tables/f.scrub.0.1.RA.csv"
groupfile<-"metadata/R_meta.csv"
options(stringsAsFactors = F)
df<-read.csv(myinput,header = T,
             row.names = 1,check.names=F)
group<-read.csv(groupfile,header = T)
outord<-group$Trt_Dis[!duplicated(group$Trt_Dis)]
df<-df[,group$sample.id]
colnames(df)<-group$Trt_Dis

df1<-as.data.frame(t(apply(df, 1, function(x) tapply(x, colnames(df), mean))))
df1<-df1[,outord,drop=F]
df2<-cbind(c(rownames(df1)),df1)
colnames(df2)[1]<-'tax'
#write.table(df2,paste(myinput,".mean.txt",sep=''),row.names = F,sep = '\t',quote = F)


###
### (5) taxonomic bar charts of relative abundance using mean values
###

rm(list=ls())
group<-read.csv("metadata/meta_mean_plot.csv")
myinput<-"taxa_tables/f.scrub.0.1.RA.mean.txt"
df1<-read.delim(myinput,header = T,row.names = 1, check.names=FALSE)
#df1<-df1[,group[,1]]
df1$tax<-rownames(df1)
df3<-melt(df1)
df3$value<-df3$value*100
head(df3)

colnames(df3)[2]<-'sample'
colset<-read.csv("metadata/color_family.csv", header=T)
colnames(colset)[1:2]<-c("tax", "Color")
colset$tax<-factor(colset$tax, levels=unique(colset$tax))
colset<-colset%%arrange(tax)
mycol<-colset$Color
names(mycol)<- unique(colset$tax)
df4<-merge(df3, colset, by="tax")
colnames(group)[1]<-'sample'
df5<-merge(df4, group, by="sample")
df5$tax<-factor(df5$tax,colset$tax)
df5$sample<-factor(df5$sample,levels = group[,1])
df5$Treatment <- factor(df5$Treatment,levels = c("50w", "90w","120w"),ordered = TRUE)
ggplot(data=df3, aes(x=sample, y=value, fill=tax)) +
  geom_bar(stat="identity",width=0.5, color="black")+
  guides(fill = guide_legend(ncol=2,title = "Family"))+
  scale_fill_manual(values = mycol) +  
  facet_wrap(~Treatment, scales= "free_x", nrow=1) +
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 65), axis.text=element_text(size=16), 
    legend.title= element_text(size=16), legend.text = element_text(size=10),
    legend.key.size = unit(6, 'mm'), panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1),axis.title = element_text(size=20,face='bold'),strip.text=element_text(size=20)
  )+
  scale_y_continuous(expand = c(0,0))+
  labs(y = "Relative abundance(%)",x="Group")

#ggsave('myinput.pdf',width = 10,height = 6)


###
### (6) statistical analysis of relative abundance using Kruskal wallis with FDR correction
###

rm(list=ls())
input<- 'taxa_tables/f.scrub.0.1.RA.csv'
groupfile<-'metadata/R_meta.csv'
df1<-read.csv(input,header = T,check.names = F)
colnames(df1)[1]<-'level3'
df2<-melt(df1,id.vars = 'level3')
df2$value<-df2$value*100
groups1<-read.csv(groupfile,header = T)
colnames(groups1)<-c('sample','mygroup','distance','trt_2')
df3<-merge(groups1,df2,by.x='sample',by.y='variable')
df3$mygroup<-factor(df3$mygroup,levels=unique(groups1$mygroup))
df3$level3<-factor(df3$level3)
allkw<-as.data.frame(compare_means(value~mygroup,data=df3,
                                   group.by = 'level3',method = 'kruskal', p.adjust.method = "fdr"))
allkw2<-allkw%% filter(p<0.05)
stat.test1<-as.data.frame(compare_means(as.formula(paste('value','~mygroup')),data = df3,
                                        group.by = 'level3', p.adjust.method = "fdr"))
stat.test2<-stat.test1 %% filter(p<0.05)

#write.csv(stat.test2,"Output_MG/stat.test2.csv")
#write.csv(stat.test1,"Output_MG/stat.test1.csv")
#write.csv(allkw,"Output_MG/allkw.csv")
#write.csv(allkw2,"Output_MG/allkw2.csv")

