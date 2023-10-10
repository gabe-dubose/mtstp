#!/usr/bin/env Rscript

#load data
data <- read.csv('data/counts_tables/dpl_tpm_counts_kallisto.csv')

#get euclidian distance matrix
euclidian.matrix <- dist(data, method='euclidian')
euclidian.matrix <- as.matrix(euclidian.matrix, labels=TRUE)
colnames(euclidian.matrix) <- rownames(euclidian.matrix) <- data[['X']]
#write to file
write.csv(euclidian.matrix, 'data/distance_matricies/euclidian.csv')

#get Manhattan distance matrix
manhattan.matrix <- dist(data, method='manhattan')
manhattan.matrix <- as.matrix(manhattan.matrix, labels=TRUE)
colnames(manhattan.matrix) <- rownames(manhattan.matrix) <- data[['X']]
#write to file
write.csv(manhattan.matrix, 'data/distance_matricies/manhattan.csv')