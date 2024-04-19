#!/usr/bin/env Rscript

#load data
data <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/dpl_tpm_counts_kallisto.csv')
data <- subset(data, substr(X, 8, 8) != 'i')

#get euclidian distance matrix
euclidian.matrix <- dist(data, method='euclidian')
euclidian.matrix <- as.matrix(euclidian.matrix, labels=TRUE)
colnames(euclidian.matrix) <- rownames(euclidian.matrix) <- data[['X']]
#write to file
write.csv(euclidian.matrix, '/home/gabe/Desktop/mtstp/data/intermediate_data/distance_matricies/euclidian.csv')

#get Manhattan distance matrix
manhattan.matrix <- dist(data, method='manhattan')
manhattan.matrix <- as.matrix(manhattan.matrix, labels=TRUE)
colnames(manhattan.matrix) <- rownames(manhattan.matrix) <- data[['X']]
#write to file
write.csv(manhattan.matrix, '/home/gabe/Desktop/mtstp/data/intermediate_data/distance_matricies/manhattan.csv')

#get distance matrix for grouped genes only 
grouped.genes.data <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/grouped_genes_only_data_distant.csv')
#get Manhattan distance matrix
manhattan.matrix <- dist(grouped.genes.data, method='manhattan')
manhattan.matrix <- as.matrix(manhattan.matrix, labels=TRUE)
colnames(manhattan.matrix) <- rownames(manhattan.matrix) <- grouped.genes.data[['sample.id']]
#write to file
write.csv(manhattan.matrix, '/home/gabe/Desktop/mtstp/data/intermediate_data/distance_matricies/grouped_only_manhattan_distant.csv')

#get distance matrix for singleton genes only
singleton.genes.data <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/singleton_genes_only_data_distant.csv')
#get Manhattan distance matrix
manhattan.matrix <- dist(singleton.genes.data, method='manhattan')
manhattan.matrix <- as.matrix(manhattan.matrix, labels=TRUE)
colnames(manhattan.matrix) <- rownames(manhattan.matrix) <- singleton.genes.data[['sample.id']]
#write to file
write.csv(manhattan.matrix, '/home/gabe/Desktop/mtstp/data/intermediate_data/distance_matricies/singleton_only_manhattan_distant.csv')

