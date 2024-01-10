#!/usr/bin/env Rscript

#load distance matrix
distances <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/distance_matricies/manhattan.csv')
#rename x to sampe.id
colnames(distances)[1] <- 'sample.id'
metadata <- read.csv('/home/gabe/Desktop/mtstp/data/experiment_metadata/mtstp_analysis_metadata.tsv', sep='\t')
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
#save coordinates
pca.coords <- data.frame(uninf.pca$x)
#get loadings
loadings <- uninf.pca$rotation
most.explained.first.pc <- loadings[, 1]
most.explained.first.pc.sorted <- sort(abs(most.explained.first.pc), decreasing = TRUE)
most.explained.first.pc.sorted
#write PCA results
write.csv(pca.coords, '/home/gabe/Desktop/mtstp/data/intermediate_data/pca/manhattan_uninfected_pca.csv')



# Sample dataset (Replace this with your dataset)
data <- iris[, -5]  # Excluding the species column for this example

# Perform PCA
pca_result <- prcomp(data, scale. = TRUE)  # Scaling for standardization

# Extract the variance explained by each principal component
variance_explained <- pca_result$sdev^2
total_variance <- sum(variance_explained)

# Calculate the proportion of variance explained by each component
variance_proportion <- variance_explained / total_variance

# Get the loadings (variable contributions to principal components)
loadings <- pca_result$rotation

# Find the variables that explain the most variance in the first principal component
most_var_explained <- loadings[, 1]  # Considering the first principal component
sorted_loadings <- sort(abs(most_var_explained), decreasing = TRUE)

# Print the variables that explain the most variance in the first principal component
top_variables <- names(sorted_loadings[1:5])  # Change 5 to the desired number
print(top_variables)








#perform PCA for grouped genes
#load distance matrix
distances <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/distance_matricies/grouped_only_manhattan_distant.csv')
#rename x to sampe.id
colnames(distances)[1] <- 'sample.id'
metadata <- read.csv('/home/gabe/Desktop/mtstp/data/experiment_metadata/mtstp_analysis_metadata.tsv', sep='\t')
#remove infected data
infected.ids <- metadata[metadata$infection.status == 'infected',]$sample.id
distances <- distances[,!names(distances) %in% infected.ids]
distances <- distances[!distances$sample.id %in% infected.ids,]
#remove sample id column
distances <- t(distances[2:length(distances)])

#run PCA
grouped.pca <- prcomp(distances)
#print summary
print(summary(grouped.pca))

grouped.pca.coords <- data.frame(grouped.pca$x)
#write PCA results
write.csv(grouped.pca.coords, '/home/gabe/Desktop/mtstp/data/intermediate_data/pca/grouped_manhattan_pca_distant.csv')


#perform PCA for singleton genes
#load distance matrix
distances <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/distance_matricies/singleton_only_manhattan_distant.csv')
#rename x to sampe.id
colnames(distances)[1] <- 'sample.id'
metadata <- read.csv('/home/gabe/Desktop/mtstp/data/experiment_metadata/mtstp_analysis_metadata.tsv', sep='\t')
#remove infected data
infected.ids <- metadata[metadata$infection.status == 'infected',]$sample.id
distances <- distances[,!names(distances) %in% infected.ids]
distances <- distances[!distances$sample.id %in% infected.ids,]
#remove sample id column
distances <- t(distances[2:length(distances)])

#run PCA
single.pca <- prcomp(distances)
#print summary
print(summary(single.pca))

single.pca.coords <- data.frame(single.pca$x)
#write PCA results
write.csv(single.pca.coords, '/home/gabe/Desktop/mtstp/data/intermediate_data/pca/single_manhattan_pca_distant.csv')
