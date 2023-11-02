#!/usr/bin/env Rscript

#load distance matrix
distances <- read.csv('data/distance_matricies/manhattan.csv')
#rename x to sampe.id
colnames(distances)[1] <- 'sample.id'
metadata <- read.csv('data/mtstp_analysis_metadata.tsv', sep='\t')
#remove infected data
infected.ids <- metadata[metadata$infection.status == 'infected',]$sample.id
distances <- distances[,!names(distances) %in% infected.ids]
distances <- distances[!distances$sample.id %in% infected.ids,]
#remove sample id column
distances <- t(distances[2:length(distances)])

#run PCA
uninf.pca <- prcomp(distances)
#print summary
print(summary(uninf.pca))

pca.coords <- data.frame(uninf.pca$x)

#write PCA results
write.csv(pca.coords, 'data/manhattan_uninfected_pca.csv')
