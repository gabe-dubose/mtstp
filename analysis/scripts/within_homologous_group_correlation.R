#A function to extract protein ids from NCBI gene ids
# input[id]: An NCBI gene id that has the corresponding protein id nested
#            between the second and third underscores
# output[protein.id]: The extracted protein id
get.prot.id <- function(id) {
  protein.id <- unlist(strsplit(id, "_"))[3]
  return(protein.id)
}

#A function to mean-center each row in a data frame
# input[df]: A data frame with numeric values
# output[df] : A data frame where each row has been mean centered
mean.center <- function(df) {
  #apply normalization to each row
  normalized.df <- t(apply(df, 1, function(row) {
    #get row mean and standard deviation
    row.mean <- mean(row, na.rm = TRUE)
    #return row minus mean
    (row - row.mean)
  }))
}

#A function to mean-center and standardize (by standard deviation) each row in a data frame
# input[df]: A data frame with numeric values
# output[df] : A data frame where each row has been standardized
mean.center.standardize <- function(df) {
  #apply normalization to each row
  normalized.df <- t(apply(df, 1, function(row) {
    #get row mean and standard deviation
    row.mean <- mean(row, na.rm = TRUE)
    row.sd <- sd(row, na.rm = TRUE)
    #if standard deviation is 0, throw warning and return unchanged row
    if (row.sd == 0) {
      warning("Warning: The standard deviation is zero. Returning mean centered values (1s).")
      (row - row.mean)
    } else {
      (row - row.mean) / row.sd
    }
  }))
}

#A function to compare phylogenetic divergence to their expression pattern divergence within gene families
# input[gene.clusters] : Path to JSON file, where each key is the name of a gene cluster 
#                        and the values are a list of the ids in the corresponding 
#                        cluster. File must be readable by rjson.
# input[expression.data] : A data frame with sequence ids in a column called "genes."
#                          The subsequent columns will be used for clustering.
#                          Values must be numeric.
# output[expression.pattern.diversity.metrics] : A data frame that contains the computed
#                                                expression pattern diversity metrics.
# input[phylogenies.dir.path]: The path to a directory containing phylogenetic trees in newick format
# input[file.extension]: The file extension for phylogeny files
compute.phylogeny.expression.correlations <- function(gene.clusters, expression.data, phylogenies.dir.path, file.extension) {
  #initialize vector to store data
  standard.data <- c()
  signif.values <- c()
  group.size <- c()
  group.pd <- c()
  #load sequence clusters
  sequence.clusters <- gene.clusters
  #sequence.clusters <- rjson::fromJSON(file=gene.clusters)
  cluster.names <- names(sequence.clusters)
  #get list of phylogeny files
  phylogeny.files <- list.files(phylogenies.dir.path)
  #calculate expression and phylogenetic distances
  for (id in names(sequence.clusters)) {
    #get the corresponding values from the named list
    cluster.id <- sequence.clusters[[id]]
    #select rows from dataframe
    cluster.values <- subset(expression.data, genes %in% cluster.id)
    #reassign first column
    row.names(cluster.values) <- cluster.values$genes
    cluster.values <- cluster.values[, -1]
    #remove rows with only 0s
    cluster.values <- cluster.values[rowSums(cluster.values != 0) > 0, ]
    #perform mean centering and standardization
    cluster.values.standardized <- mean.center.standardize(cluster.values)
    #if there are enough genes to cluster, do so
    if (nrow(cluster.values) > 1) {
      #get distance for standardized values
      standardized.expression.distance <- as.matrix(dist(cluster.values.standardized, diag = TRUE))
      #get phylogenetic distance
      cluster.gene.tree.file <- paste(id, file.extension, sep="")
      #check if the cluster has a gene tree associated
      if (cluster.gene.tree.file %in% phylogeny.files) {
        #assemble path to file
        cluster.gene.tree.path <- paste(phylogenies.path, "/", cluster.gene.tree.file, sep="")
        #read tree
        gene.tree <- ape::read.tree(cluster.gene.tree.path)
        #phylogenetic distance matrix
        phylogenetic.distance.matrix <- as.matrix(ape::cophenetic.phylo(gene.tree))
        
        if (length(phylogenetic.distance.matrix) == length(standardized.expression.distance)) {
          #reorder 
          primary.headers <- rownames(phylogenetic.distance.matrix)
          #standardized
          standardized.expression.distance.headers <- rownames(standardized.expression.distance)
          standardized.expression.distance.headers.re <- standardized.expression.distance.headers[match(primary.headers, standardized.expression.distance.headers)]
          standardized.expression.distance <- standardized.expression.distance[match(standardized.expression.distance.headers.re, standardized.expression.distance.headers), 
                                                                               match(standardized.expression.distance.headers.re, standardized.expression.distance.headers)]
          
          #run mantel test
          #standardized
          mantel.stadardized <- vegan::mantel(as.matrix(phylogenetic.distance.matrix), as.matrix(standardized.expression.distance), permutations = 0)
          corr.coeff.standardized <- mantel.stadardized$statistic
          significance <- mantel.stadardized$signif
          standard.data <- c(standard.data, corr.coeff.standardized)
          signif.values <- c(signif.values, significance)
          
          #calculate branch lengths
          branch.lengths <- gene.tree$edge.length
          #calculate phylogenetic diversity
          sum.of.lengths <- sum(branch.lengths)
          group.pd <- c(group.pd, sum.of.lengths)
          group.size <- c(group.size, length(gene.tree$tip.label))
        } 
      }
    }
  }
  #assemble output
  correlation.results <- list("r" = standard.data,
                              'p' = signif.values,
                              "group.size" = group.size,
                              "group.pd" = group.pd)
  return(correlation.results)
  
}

#load total expression data
total.expression.data <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/dpl_tpm_counts_kallisto.csv')
#load metadata
metadata <- read.csv('/home/gabe/Desktop/mtstp/data/experiment_metadata/mtstp_analysis_metadata.tsv', sep='\t')
#remove infected data
total.expression.data <- total.expression.data[substr(total.expression.data$X, 8, 8) != 'i', ]
metadata <- metadata[substr(metadata$sample.id, 8, 8) != 'i', ]
#get developmental stage only
dev.metadata <- data.frame(metadata[, c('developmental.stage', 'sample.id')])

#get medians for each stage
total.data <- merge(total.expression.data, dev.metadata, by.x = "X", by.y = "sample.id")

expr.values <- total.data[, -which(names(total.data) == "developmental.stage")]
median.values <- aggregate(expr.values, by = list(developmental.stage = total.data$developmental.stage), FUN = median)
median.values <- median.values[, !names(median.values) %in% "X"]
median.values <- t(median.values)
#reformat a little
colnames(median.values) <- as.character(unlist(median.values[1, ]))
median.values <- median.values[-1, ]

gene.ids <- data.frame(rownames(median.values))
median.values <- cbind(gene.ids, median.values)
rownames(median.values) <- NULL
names(median.values)[names(median.values) == "rownames.median.values."] <- "genes"
median.values <- median.values[, c("genes", "third-instar", "fifth-instar", "early-pupa", "late-pupa", "adult")]

median.values$`third-instar` <- as.numeric(median.values$`third-instar`)
median.values$`fifth-instar` <- as.numeric(median.values$`fifth-instar`)
median.values$`early-pupa` <- as.numeric(median.values$`early-pupa`)
median.values$`late-pupa` <- as.numeric(median.values$`late-pupa`)
median.values$`adult` <- as.numeric(median.values$`adult`)

#convert gene ids to protein ids
median.values$genes <- sapply(median.values$genes, get.prot.id)

#run
gene.clusters.file <- "/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/psiblast_id-30_e-neg-10_cov-1_sequence_clusters.json"
gene.clusters.data <- rjson::fromJSON(file=gene.clusters.file)
phylogenies.path <- '/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/inferred_phylogenies_close'
file.extension <- '_alignment.fasta.treefile'
distance.correlations <- compute.phylogeny.expression.correlations(gene.clusters.data, median.values, phylogenies.path, file.extension)
distance.correlations.df <- data.frame(distance.correlations)
write.csv(distance.correlations.df, "/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/phylogenetic_vs_expression_distance.csv")

#test that distribution is different from 0
shapiro.test(distance.correlations.df$r)
t.test(distance.correlations.df$r, mu = 0)
wilcox.test(distance.correlations.df$r, mu = 0, alternative = "two.sided")



