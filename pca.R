library(ggplot2)

setwd("/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data")

sample_meta <- read.table("2023_08_16-1000G_P3_subset_map", header = FALSE)
colnames(sample_meta) <- c("ID","POP","SUP_POP","SEX")

eigenval <- read.table("ALLCHR.phase3_v5.shapeit2_mvncall_integrated_allFilters_int_subset.eigenval", header = FALSE)
eigenvec <- read.table("ALLCHR.phase3_v5.shapeit2_mvncall_integrated_allFilters_int_subset.eigenvec", header = FALSE)
colnames(eigenvec) <- c("V1","ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

eigenvec_metadata <- merge(eigenvec, sample_meta)


plot(data= eigenvec, PC1~PC2)

pca <- ggplot(eigenvec_metadata, aes(x=PC1, y=PC2, label = POP, color = POP)) + geom_point(size=2) + ggtitle("PC1 vs PC2") +  theme(plot.title = element_text(hjust = 0.5))

jpeg(file = "/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/1000G_pca.jpeg")
pca
dev.off()
