#load data
data <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/dpl_tpm_counts_kallisto.csv')
data <- subset(data, substr(X, 8, 8) != 'i')
rownames(data) <- data[,1]
data <- data[, -1]
data <- t(data)

#calculate pearson correlation matrix
pearson.correlation.matrix <- cor(data, method="pearson")
pearson.correlation.matrix <- (1 - pearson.correlation.matrix)/2
#run PCA
uninf.pca <- prcomp(pearson.correlation.matrix)
#print summary
print(summary(uninf.pca))
#save coordinates
pca.coords <- data.frame(uninf.pca$x)
#write PCA results
write.csv(pca.coords, '/home/gabe/Desktop/mtstp/data/intermediate_data/pca/pearson_uninfected_pca.csv')
