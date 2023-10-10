########
# RUN IN R
########
setwd("/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/bayenv_files")

# BayEnv results
i=6
data1=read.delim(paste('8Jun20_BayEnv_results_chr',i,'.txt',sep=''))
data2=read.delim(paste('8Jun20_BayEnv_region_results_chr',i,'.txt',sep=''))

for (i in 2:22){
tmp1=read.delim(paste('8Jun20_BayEnv_results_chr',i,'.txt',sep=''))
tmp2=read.delim(paste('8Jun20_BayEnv_region_results_chr',i,'.txt',sep=''))
names(tmp1)<-names(data1)
names(tmp2)<-names(data2)
data1=rbind(data1,tmp1)
data2=rbind(data2,tmp2)
print(i)
}

# PBS and iHS results summarized by region
both2_PBS=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/PBS_files/21May20_PBS_ALLchr_outlier_count.txt',sep=' ')
both2_iHS=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files/26May20_iHS_ALLchr_outlier_count.txt',sep=' ')

all=merge(both2_PBS,both2_iHS,by='region',all.x=T)
all=merge(all,data2,by='region',all.x=T)

genes=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files/26May20_iHS_test_gene_info.bed',header=F)
genes$region=paste(genes$V1,genes$V2,genes$V3,sep='_')
both2=merge(all,genes,by='region')

# iHS outlier region cutoff
quantile(all$outlier_SNPs.y,seq(0,1,0.01),na.rm=T)[100]

# PBS outlier region cutoff
quantile(all$outlier_SNPs.x,seq(0,1,0.01))[100]

# BayEnv outlier region cutoff
quantile(all$n_outliers,seq(0,1,0.01))[100]

write.table(subset(both2,outlier_SNPs.x>5 & outlier_SNPs.y>14),'outlier_regions.txt',row.names=F,sep='\t')
write.table(all,'outliers_by_regions.txt',row.names=F,sep='\t')
