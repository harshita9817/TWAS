library(data.table)
library(ggplot2)
library(dplyr)

gwas.datTWAS <- fread("TWAS.CMLM_FT_2021.csv", data.table = F)

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))


colnames(gwas.datTWAS)

nCHR <- length(unique(gwas.datTWAS$Chromosome))
gwas.datTWAS$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.datTWAS$Chromosome))){
  nbp[i] <- max(gwas.datTWAS[gwas.datTWAS$Chromosome == i,]$`Position.`)
  gwas.datTWAS[gwas.datTWAS$Chromosome == i,"BPcum"] <- gwas.datTWAS[gwas.datTWAS$Chromosome == i,"Position."] + s
  s <- s + nbp[i]
}

axis.set <- gwas.datTWAS %>% 
  group_by(Chromosome) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

# to set the limite of y axis
if (max(-log10(gwas.datTWAS$FDR)) > -log10(0.05) ) {
  limite = ceiling(max(-log10(gwas.datTWAS$FDR)))
} else {
  limite = ceiling(-log10(0.05))
}
gTWAS.fdr <- ggplot() + 
  geom_point(data=gwas.datTWAS, aes(BPcum, -log10(FDR), colour=factor(Chromosome, levels = c(1:10))), size = 2) + 
  geom_hline(yintercept = -log10(0.05), linetype=2) + 
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  #annotate("text", label="paste(italic(MADS69))", y=12, x=713290558, parse=T, size=5) +
  scale_x_continuous(label = axis.set$Chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](FDR))) + 
  xlab("Chromosome") + 
  scale_y_continuous(limits = c(0, limite),
                     breaks = seq(1, limite, 2))
#+ labs(title="", subtitle=paste0(trait))

gTWAS.fdr

ggsave(plot=gTWAS.fdr, paste0("Manhathan/Manhathan_",".png"), width = 7, height = 4.5)

dev.off()





gTWAS.fdr <- ggplot() + 
  rasterise(geom_point(data=gwas.dat.RNA, aes(BPcum, -log10(gene), colour=factor(CHROM, levels = c(1:10))), size = 2),dpi=600) + 
  geom_hline(yintercept = -log10(0.05/snps), linetype=2) + 
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  geom_point(aes(y=gene.pos.y, x=gene.pos.x), fill="red", size = 5, shape = 24) +
  geom_point(aes(y=gene.pos.y, x=hen1pos), fill="black", size = 5, shape = 24) +
  geom_point(aes(y=gene.pos.y, x=rgd1pos), fill="white", size = 5, shape = 24) +
  #annotate("text", label=gene.name2, y=gene.pos.y, x=gene.pos.x, parse=T, size=5, hjust = -0.01) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](p-value))) + 
  xlab("Chromosome") + 
  scale_y_continuous(limits = c(2, limite),
                     breaks = seq(0, limite, 5))

gTWAS.fdr

















