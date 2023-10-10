library(rehh)
library(data.table)
library(R.utils)
library(vcfR)
library(tidyr)
library(dplyr)

#Loading Example Files
#make.example.files()
#data(haplohh_cgu_bta12)

setwd("/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data")
# My Data

for (i in 6:6) {
	hh <- data2haplohh(hap_file = paste("/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/ALL.chr",i,".phase3_v5.phased_subset_noindels_withAA.vcf.gz",sep=''),polarize_vcf = TRUE,min_perc_geno.mrk=90,remove_multiple_markers=TRUE)
	res<-scan_hh(hh,threads=4,phased = TRUE)
	res.ihs<-ihh2ihs(res)
	res.ihs.df2<-as.data.frame(res.ihs$ihs)
	res.ihs.df2<-res.ihs.df2[complete.cases(res.ihs.df2),]
	write.table(res.ihs.df2,paste("/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/2023_Aug_27_iHS_chr",i,".txt",sep=''),row.names=F,col.names=F,sep='\t',quote=F)
	print(i)
}

png(file= "/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/iHS_chr6.png")
par(mfrow=c(1,2))
manhattanplot(res.ihs.df2)
manhattanplot(res.ihs.df2, pval = TRUE)
dev.off()

sink( paste0( myDir, "iHS_sessionInfo.txt"))
sessionInfo()
sink()


#Sort PBS values and take top 1% most significant
IHS_absval<-fread("2023_Aug_27_iHS_chr6.txt",header=T)
colnames(IHS_absval)[5] <- "absIHS"
IHS_sorted <- IHS_absval %>% arrange(-absIHS)
IHS_top1percent <- quantile(IHS_sorted$absIHS, 0.99, na.rm=TRUE)
quantile(all$NumofSigSNPs.x,seq(0,1,0.01))[100]IHS_outliers <- IHS_sorted[IHS_sorted$absIHS > IHS_top1percent,]
write.table(IHS_outliers,'2023-10-04_IHS_CHR6_outliers.txt',row.names=F,sep='\t',quote=F)

#Total number of SNPS in each 50 kb window
binannotation<-fread("50KBbins.txt",header=T)

IHSSNPsbins <- matrix(0, nrow = nrow(binannotation), ncol = 1) 

for (i in 1:nrow(binannotation)){
test <- c(which(IHS_sorted$POS> binannotation$Start[i] & IHS_sorted$POS<binannotation$End[i]))
IHSSNPsbins[i,] <- length(test)
}

IHSSNPsbins<- as.data.frame(IHSSNPsbins)
IHSSNPsperbin <- cbind(binannotation,IHSSNPsbins)
colnames(IHSSNPsperbin)[5] = "TotNumofSNPs"

#Num of Significant SNPs in each 50 kb window

IHSsigSNPsbins <- matrix(0, nrow = nrow(binannotation), ncol = 1) 

for (i in 1:nrow(binannotation)){
test <- c(which(IHS_outliers$POS> binannotation$Start[i] & IHS_outliers$POS<binannotation$End[i]))
IHSsigSNPsbins[i,] <- length(test)
}

IHSsigSNPsbins<- as.data.frame(IHSsigSNPsbins)
IHSsigSNPsperbin <- cbind(IHSSNPsperbin,IHSsigSNPsbins)
colnames(IHSsigSNPsperbin)[6] = "NumofSigSNPs"

#Exclude bins with less than 20 total SNPs
IHSsigSNPsperbinFinal <- IHSsigSNPsperbin[IHSsigSNPsperbin$TotNumofSNPs>20,]
write.table(IHSsigSNPsperbinFinal,'2023-10-03_IHS_chr6_outlier_count.txt',row.names=F,sep='\t',quote=F)
