dt = read.table("../00.data/reads_num.profile.clu",sep="\t",header=T, row.names=1, check.names=F)
dt = dt[-1,]
dt = t(t(dt)/colSums(dt))
map  = read.table("../00.data/sample.group.txt", sep="\t", header=T, check.names=F)

dtt = t(dt)
dtt = dtt[map[,1],]
nsample = ncol(dtt)
names = colnames(dtt)
result = rbind()

for (i in 1:nsample){
  tt = t.test(dtt[,i] ~ map[,2])
  mean_c = tt$estimate[1]
  mean_ibd = tt$estimate[2]
  temp = data.frame(name=names[i], t.test=tt$p.value, mean_c=mean_c, mean_ibd=mean_ibd)
  result = rbind(result, temp)
}

#write.table(result, "clipboard", sep="\t", quote = F)

library(vegan)
library(ggpubr)
#dtt = dtt[,-1]

ado = adonis(dtt ~ map$group2, method='bray')
r2 = round(ado$aov.tab$R2[1], digits = 6)
prf = ado$aov.tab$`Pr(>F)`[1]
n = ncol(ado$coef.sites)
m = ado$aov.tab$Df[1]
radj = RsquareAdj(r2, n=n, m=m)




#------------PCoA 计算------------------------
otu.dist <- vegdist(dtt,method="bray", binary=F)
otu_pcoa<- cmdscale(otu.dist,eig=TRUE, k=4)
pc12 <- otu_pcoa$points[,1:4]
pc_importance<-round(otu_pcoa$eig/sum(otu_pcoa$eig)*100,digits = 2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
data <- merge(pc12, map, by="samples", by.y='sample')
data$group = as.factor(data$group2)

paletteColor = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","#666666")
#------------画图------------------------
p = ggscatter(data, x='V1', y='V2',color='group2',
              ellipse = T, ellipse.level=0.8,
              star.plot=T,
              palette = paletteColor)+
  labs(x=paste("PCoA 1 (", pc_importance[1],digits=4,"%)", sep=""), 
       y=paste("PCoA 2 (", pc_importance[2],digits=4, "%)", sep=""),
       title=paste("PCoA bray_curtis Species\nR2       = ",r2, "\nPr(>F) = ", prf,sep=""))  + #绘制点图并设定大小
  theme_bw()


p

dd = read.table("clipboard", sep="\t", header=F)

library(ggsignif)
kk = list(c("HC","IBD_stage1"), c("HC", "IBD_stage2"), c("IBD_stage1", "IBD_stage2"))
p3 = ggplot(data=dd, aes(x=V1, y=V3, fill=V1))+
  geom_boxplot()+
  theme_bw()+
  stat_compare_means(comparisons = kk)+
  theme(axis.text.x=element_text(angle=45, hjust=1))

kk = list(c("HC", "IBD"))
p4 = ggplot(data=dd, aes(x=V2, y=V3, fill=V2))+
  geom_boxplot()+
  stat_compare_means(comparisons = kk)+
  theme_bw()

library(patchwork)
p1+p2+p4+p3
