library(data.table)
library(ggplot2)

setwd('/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data')


# concat RFmix proportions
#  the most likely assignment of subpopulations per CRF point 
data1_msp=fread('/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/RFmix_chr1.msp.tsv', skip=1)

for (f in 2:22) { 
data1_part<-fread(paste('RFmix_chr',f,'.msp.tsv',sep=''),skip=1)
data1_msp=rbind(data1_msp,data1_part) }
dim(data1_msp)
#Make New Columns which sum the number of 0 or 1 (CEU Or YRI) calls per genomic region
data1_msp$count_0<-apply(data1_msp[,c(7:336),with=F],1,function(x) length(which(x==0)))
data1_msp$count_1<-apply(data1_msp[,c(7:336),with=F],1,function(x) length(which(x==1)))
write.table(data1_msp[,c(1:6,337:338),with=F],'/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/2023_Aug_26_RFmix_mspsummary.txt',row.names=F,sep='\t')

 
# admixture proportions per chromosome
data1_Q=fread('/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/RFmix_chr1.rfmix.Q',skip=1)
data1_Q$chr <- 1

for (f in 2:22) {
data1_part<-fread(paste('RFmix_chr',f,'.rfmix.Q',sep=''),skip=1)
data1_part$chr <- f
data1_Q=rbind(data1_Q,data1_part) }

write.table(data1_Q,'2023_Aug_23_RFmix_Qsummary.txt',row.names=F,sep='\t')

# compare ADMIXTURE vs RFmix proportions
#setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")
#data1=read.delim('20Aug22_RFmix_vs_ADMIXTURE.txt')
#I think I need to make a file which takes the average ancestry of all chromosomes for RFmix and compare it to ADMIXTURE
data2_Q=fread('/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/RFmix_chr1.rfmix.Q',skip=1)
data3_Q <- data2_Q[,-3]

for (f in 2:22) {
data3_part<-fread(paste('RFmix_chr',f,'.rfmix.Q',sep=''),skip=1)
data3_part <- data3_part[,2]
data3_Q=cbind(data3_Q,data3_part) }
colnames(data3_Q) <- c('Ind ID','Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19','Chr20','Chr21','Chr22')

data3_Q$mean <- rowMeans(data3_Q[,2:23])

RFmix_wholegenome_ancestry <- data3_Q[,c(1,24)]


AD_sample_order <- fread('/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/ADMIXTURE_sampleorder.txt',header=F)
AD_2_Q <- fread('/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_allFilters_subset.2.Q',header=F)
AD_wholegenome_ancestry <- cbind(AD_sample_order,AD_2_Q[,2])
colnames(AD_wholegenome_ancestry) <- c('Ind ID','CEU')

RFmix_AD_compare <- merge(RFmix_wholegenome_ancestry, AD_wholegenome_ancestry, by='Ind ID')
colnames(RFmix_AD_compare) <- c('Ind ID','RFmix','ADMIXTURE')

# plot
png(file= "/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/RFmix_vs_ADMIXTURE.png")
par(mfrow=c(1,2))
plot(RFmix_AD_compare$RFmix,RFmix_AD_compare$ADMIXTURE,xlab='European ancestry proportion (RFMix)',ylab='European ancestry proportion (ADMIXTURE)',xlim=c(0,1),ylim=c(0,1),bty='n')


cor.test(RFmix_AD_compare$RFmix,RFmix_AD_compare$ADMIXTURE)

plot(density(RFmix_AD_compare$RFmix),main='',bty='n',lwd=2,xlab='European ancestry proportion')
lines(density(RFmix_AD_compare$ADMIXTURE,na.rm=T),lty=2,lwd=2)
dev.off()

