dt = read.table("./xx.txt",sep="\t", header=T,check.names=F)
library(ggplot2)
dt1 = dt
dt1$lfdr = -log10(dt1$lfdr)
dt1$fold_change = ifelse(dt1$enriched=="Control", 0-(dt1$fold_change), (dt1$fold_change))


ggplot(dt1, aes(x=fold_change,
                y=lfdr,
                colour=color
                #alpha=color,
                #shape=color
                ))+
  geom_point(size=3)+
  geom_text(aes(x=fold_change+10, y=lfdr, label=label))+
  theme_bw()+
  #geom_vline(xintercept=log(5), linetype=5)+
  #geom_vline(xintercept=-log(5), linetype=5)+
  geom_hline(yintercept = -log10(0.05), linetype=5)+ # qval
  geom_hline(yintercept = -log10(0.001), linetype=5)+ # qval
  scale_color_manual(values=c("#d53e4f", "gray", "#3288bd"))+
  scale_alpha_manual(values=c(1,1,0.5))+
  scale_shape_manual(values=c(16,16,1))+
  xlab("log2(fold_change)")+
  ylab("-log10(lfdr)")
