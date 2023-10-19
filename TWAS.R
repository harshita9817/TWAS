(still working on this)

source("http://zzlab.net/GAPIT/gapit_functions.txt") 

pheno<- read.csv("bluesHN_2021.csv")
geneinfo<-read.csv("geneinfo2.csv")
tpm<- read.csv("tpm_ranged_0-2.csv")
#in all file you have an X an the beginning that is causing trouble
colnames(tpm)[1:10]
colnames(pheno)
colnames(geneinfo)
geneinfo$X<- NULL
tpm$X<-NULL
#re-naming  the columns names so that it is easy to merge
colnames(pheno)[1]<- "GeneID"
colnames(tpm)[1]<- "GeneID"
colnames(tpm3)[1:10]
########################
library(dplyr)
mergedtpmpheno <- merge(tpm, pheno, by = "GeneID", all = FALSE)
tpm1 <- tpm[tpm$GeneID %in% mergedtpmpheno$GeneID, ]
pheno1 <- pheno[pheno$GeneID %in% mergedtpmpheno$GeneID, ]

##########################

###############################
#we will remove it
tpm1$X <- NULL
pheno$X <- NULL
geneinfo$X <- NULL

#check for NAs
sum(is.na(tpm1))
sum(is.na(pheno1))
sum(is.na(geneinfo))

#Quantile<-apply(tpm1[,-1],2,  # 2 indicates it is for column and 1 indicates it is for row
 #              function(A){min_x=as.numeric(quantile(A,0.05));
  #             max_x=as.numeric(quantile(A,0.95));
   #            out<-(2*(A-min_x)/(max_x-min_x));
    #           out[out>2]<-2;out[out< 0]<- 0;return(out)})

Quantile.t <- as.data.frame(tpm1)
Quantile.t$GeneID<- row.names(Quantile.t)
myGD <-  Quantile.t[,c(ncol(Quantile.t),1: (ncol(Quantile.t)-1))]

#myGD2 <- myGD[ , colSums(is.na(myGD)) == 0]


mergedtpmgene <- merge(tpm2, geneinfo, by = "GeneID", all = FALSE)
tpmfinal <- tpm[tpm$GeneID %in% mergedtpmgene$GeneID, ]
genefinal <- geneinformation[geneinformation$GeneID %in% mergedtpmgene$GeneID, ]

myGD2<- myGD[,-2]


myY <- pheno1
colnames(myY)[1] <- "taxa"
colnames(myY)[2] <- "flowering"
head(myY$taxa)

myGM<- geneinfo
unique(myGM$chr.no) #only cromosomes 
NROW(myY[myY$taxa %in% myGD$GeneID])
myGAPIT <- GAPIT(Y=myY,
                 GD=myGD,
                 GM=myGM,
                 PCA.total=3,
                 model= "CMLM",
                 SNP.MAF=0,
                 file.output=F
)
#getting the important genes and Manhattan plots
values <- data.frame(myGAPIT$GWAS)
values$FDR <- p.adjust(values$P.value,method = "BH")
write.csv(values, paste0("TWASoil.CMLM_","FT","_2021.csv"), row.names = F)

