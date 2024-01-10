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

#A function to calculate confidence intervals from a vector
get.ci <- function(data) {
  #get mean, standard deviation, and n
  sample.mean <- mean(data)
  sample.sd <- sd(data)
  n <- length(data)
  #calculate standard error
  se <- sample.sd / sqrt(n)
  error.margin <- qt(0.975, df=n-1)*se
  #calculate ci
  lower.i <- sample.mean - error.margin
  upper.i <- sample.mean + error.margin
  confidence.interval <- c(lower.i, upper.i)
  return(confidence.interval)
}

#A function to simulate a null model of how likely a random subset of genes are to
#show a greater than or equal to amount of expression pattern diversity than the input subset
# input[total.expression.data]: A data frame holding the total expression data
# input[subset.size]: Number of rows to subsample from the "total.expression.data" data frame
# input[subset.diversity]: Diversity of observed subset
# input[iterations]: An integer specifying how many times to randomly subsample the "total.expression.data" data frame
#Note: function mean centers and standardizes expressiond data
#Note: rows with only 0 values should be removed prior to use
null.expression.diversity.test <- function(total.expression.data, subset.size, observed.diversity, iterations) {
  #initialize counter for the number of times the random subsample showed equal to or greater than diversity
  s.ge.o <- 0
  for (i in 1:iterations) {
    #randomly subsample data frame
    random.subset <- total.expression.data[sample(nrow(total.expression.data), subset.size), ]    
    #adjust columns for clustering
    row.names(random.subset) <- random.subset$genes
    random.subset <- random.subset[, -1]
    #perform mean centering and standardization
    random.subset.standardized <- mean.center.standardize(random.subset)
    #perform clustering
    clust <- hclust(dist(random.subset.standardized))
    tree <- ape::as.phylo(clust)
    #calculate diversity by summing branch lengths
    branch.lengths <- tree$edge.length
    sum.of.lengths <- sum(branch.lengths)
    #evaluate the random subsample against the observed diversity
    if (sum.of.lengths >= observed.diversity) {
      s.ge.o <- s.ge.o + 1
    }
  }
  #get the probability the random subsample is greater than or equal to observed
  p.s.ge.o <- s.ge.o / iterations
  return(p.s.ge.o)
}

#A function to calculate the Tau index for a row in a data frame
# input[row]: A row or vector of values
# output[tau]: Tau index
compute.tau <- function(row) {
  #normalize to maximal component value
  row <- row / max(row)
  sum.values <- sum(1 - row)
  tau <- sum.values / (length(row) - 1)
  return(tau)
}

#A function to compare gene family phylogenetic diversity to their expression pattern diversity
# input[gene.clusters] : Path to JSON file, where each key is the name of a gene cluster 
#                        and the values are a list of the ids in the corresponding 
#                        cluster. File must be readable by rjson.
# input[expression.data] : A data frame with sequence ids in a column called "genes."
#                          The subsequent columns will be used for clustering.
#                          Values must be numeric.
# output[expression.pattern.diversity.metrics] : A data frame that contains the computed
#                                                expression pattern diversity metrics.
compute.expression.pattern.diversity <- function(gene.clusters, expression.data) {
  #load sequence clusters
  sequence.clusters <- gene.clusters
  #sequence.clusters <- rjson::fromJSON(file=gene.clusters)
  cluster.names <- names(sequence.clusters)
  
  #initialize vectors to store data
  expression.pattern.diversity <- c()
  standardized.expression.pattern.diversity <- c()
  centered.expression.pattern.diversity <- c()
  p.random <- c()
  mean.tau <- c()
  upper_ci_tau <- c()
  lower_ci_tau <- c()
  
  #iterate through clusters
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
    #perform mean centering
    cluster.values.centered <- mean.center(cluster.values)
    #perform mean centering and standardization
    cluster.values.standardized <- mean.center.standardize(cluster.values)
    
    #if there are enough genes to cluster, do so
    if (nrow(cluster.values) > 1) {
      #perform clustering on non-standardized values
      clust <- hclust(dist(cluster.values))
      tree <- ape::as.phylo(clust)
      #calculate diversity by summing branch lengths
      branch.lengths <- tree$edge.length
      sum.of.lengths <- sum(branch.lengths)
      expression.pattern.diversity <- c(expression.pattern.diversity, sum.of.lengths)
      
      #perform clustering on mean centered values
      centered.clust <- hclust(dist(cluster.values.centered))
      centered.tree <- ape::as.phylo(centered.clust)
      #calculate diversity by summing branch lengths
      centered.branch.lengths <- centered.tree$edge.length
      centered.sum.of.lengths <- sum(centered.branch.lengths)
      centered.expression.pattern.diversity <- c(centered.expression.pattern.diversity, centered.sum.of.lengths)
      
      #perform clustering on standardized values
      standard.clust <- hclust(dist(cluster.values.standardized))
      standard.tree <- ape::as.phylo(standard.clust)
      #calculate diversity by summing branch lengths
      standard.branch.lengths <- standard.tree$edge.length
      standard.sum.of.lengths <- sum(standard.branch.lengths)
      standardized.expression.pattern.diversity <- c(standardized.expression.pattern.diversity, standard.sum.of.lengths)
      
      #run random null model
      #remove 0s from expression data
      #non.0.expression.data <- expression.data[rowSums(expression.data != 0) > 0, ]
      #print(paste("Calculating null probability for cluster", id, sep=" "))
      #p.random.diversity <- null.expression.diversity.test(total.expression.data=non.0.expression.data, 
      #subset.size=nrow(cluster.values),
      #observed.diversity=standardized.expression.pattern.diversity,
      #iterations=1000)
      p.random.diversity <- 1
      p.random <- c(p.random, p.random.diversity)
      
      #calculate tau index for each gene
      tau.values <- apply(cluster.values, 1, compute.tau)
      mean.tau.value <- mean(tau.values)
      mean.tau <- c( mean.tau, mean.tau.value)
      #get confidence intervals
      ci <- get.ci(tau.values)
      upper_95_t_ci <- ci[2]
      upper_ci_tau <- c(upper_ci_tau, upper_95_t_ci)
      lower_95_t_ci <- ci[1]
      lower_ci_tau <- c(lower_ci_tau, lower_95_t_ci)
    }
    #if not, add NAs
    else {
      expression.pattern.diversity <- c(expression.pattern.diversity, "NA")
      standardized.expression.pattern.diversity <- c(standardized.expression.pattern.diversity, "NA")
      centered.expression.pattern.diversity <- c(centered.expression.pattern.diversity, "NA")
      p.random <- c(p.random, "NA")
      mean.tau <- c(mean.tau, "NA")
      upper_ci_tau <- c(upper_ci_tau, "NA")
      lower_ci_tau <- c(lower_ci_tau, "NA")
    }
  }
  #assemble and return output
  expression.pattern.diversity.metrics <- data.frame(cluster=cluster.names,
                                                     expression.pattern.diversity = expression.pattern.diversity,
                                                     centered.expression.pattern.diversity = centered.expression.pattern.diversity,
                                                     standardized.expression.pattern.diversity = standardized.expression.pattern.diversity,
                                                     p.random = p.random,
                                                     mean.tau = mean.tau,
                                                     upper_tau_ci = upper_ci_tau,
                                                     lower_tau_ci = lower_ci_tau)
  return(expression.pattern.diversity.metrics)
}

#A function to iteratively calculate phylogenetic diversity of phylogenetic trees
# input[phylogenies.dir.path]: The path to a directory containing phylogenetic trees in newick format
# input[file.extension]: The file extension for phylogeny files
# input[expression.pattern.diversity.data]: A data frame generated by the "compute.expression.pattern.diversity" function.
#                                           The values in "cluster" + file.extension must correspond to the file names
#                                           in the phylogenies.dir.path.
compute.gene.cluster.phylogenetic.diversity <- function(phylogenies.dir.path, file.extension, expression.pattern.diversity.data) {
  #initialize vectors to store output
  gene.cluster.phylogenetic.diversity <- c()
  n.genes <- c()
  #get list of files
  phylogeny.files <- list.files(phylogenies.dir.path)
  
  #iterate through each gene family in the expression diversity data
  for (cluster in expression.pattern.diversity.data$cluster) {
    cluster.gene.tree.file <- paste(cluster, file.extension, sep="")
    #check if the cluster has a gene tree associated
    if (cluster.gene.tree.file %in% phylogeny.files) {
      #assemble path to file
      cluster.gene.tree.path <- paste(phylogenies.path, "/", cluster.gene.tree.file, sep="")
      gene.tree <- ape::read.tree(cluster.gene.tree.path)
      #calculate branch lengths
      branch.lengths <- gene.tree$edge.length
      #calculate phylogenetic diversity
      sum.of.lengths <- sum(branch.lengths)
      gene.cluster.phylogenetic.diversity <- c(gene.cluster.phylogenetic.diversity, sum.of.lengths)
      #add gene family size
      n.genes <- c(n.genes, length(gene.tree$tip.label))
    }
    else {
      gene.cluster.phylogenetic.diversity <- c(gene.cluster.phylogenetic.diversity, 'NA')
      n.genes <- c(n.genes, 'NA')
    }
  }
  #preserve orgision data frame
  phylogenetic.diversity.data <- expression.pattern.diversity.data
  #add results to data
  phylogenetic.diversity.data$gene.cluster.phylogenetic.diversity <- gene.cluster.phylogenetic.diversity
  phylogenetic.diversity.data$n.genes <- n.genes
  #return results
  return(phylogenetic.diversity.data)
}

#A function to randomly generate gene clusters
random.gene.clusters <- function(gene.clusters, total.expression.data) {
  #initialize list for random gene assignments
  random.clusters <- list()
  #load gene clusters
  sequence.clusters <- gene.clusters
  #get list of available gene
  gene.ids <- c(total.expression.data$genes)
  #get cluster ids
  cluster.ids <- names(sequence.clusters)
  #iterate through each cluster
  for (cluster.id in cluster.ids) {
    cluster <- sequence.clusters[[cluster.id]]
    n <- length(cluster)
    #randomly sample n gene ids
    random.ids <- sample(gene.ids, n)
    #assign random ids to cluster in random.cluster list
    random.clusters[[cluster.id]] <- random.ids
  }
  return(random.clusters)
}

#A function to test if the correlation between gene cluster diversity and expression diversity is random
#phylo.data <- compute.gene.cluster.phylogenetic.diversity(phylogenies.path, file.extension, gene.expression.diversity)
null.test <- function(total.expression.data, gene.family.data, clusters.file, iterations) {
  
  #initialize vectors to store local observed values
  local.slopes <- list()
  local.intercepts <- list()
  local.r2 <- list()
  
  #initialize vectors to store null values
  global.null.slope.values <- c()
  global.null.intercept.values <- c()
  global.null.r2.values <- c()
  local.null.slope.values <- list()
  local.null.intercept.values <- list()
  local.null.r2.values <- list()
  
  #convert NAs
  gene.family.data[gene.family.data == "NA"] <- NA
  #gene.family.data <- subset(gene.family.data, !apply(gene.family.data == "NA", 1, any))
  #calculate correlation between standardized expression pattern diversity and gene cluster phylogenetic diversity
  gene.family.data$gene.cluster.phylogenetic.diversity <- as.numeric(gene.family.data$gene.cluster.phylogenetic.diversity)
  gene.family.data$standardized.expression.pattern.diversity <- as.numeric(gene.family.data$standardized.expression.pattern.diversity)
  global.linear.model <- lm(standardized.expression.pattern.diversity~gene.cluster.phylogenetic.diversity, data = gene.family.data, na.action = na.omit)
  global.linear.model.summary <- summary(global.linear.model)
  global.observed.slope <- global.linear.model.summary$coefficients[2]
  global.observed.intercept <- global.linear.model.summary$coefficients[1]
  global.observed.r2 <- global.linear.model.summary$r.squared
  
  #calculate local correlations
  #initialize vector to store gene family sizes that have enough points
  gene.family.sizes.n5 <- c()
  #get unique gene family sizes
  gene.family.sizes <- unique(gene.family.data$n.genes)
  #gene.family.sizes <- gene.family.sizes[gene.family.sizes != NA]
  gene.family.sizes <- gene.family.sizes[!is.na(gene.family.sizes)]
  for (gene.family.size in gene.family.sizes) {
    #get data
    local.data <- gene.family.data[gene.family.data$n.genes == gene.family.size, ]
    #remove NAs
    local.data <- na.omit(local.data)
    #check if there are at least five points
    if (nrow(local.data) >= 5) {
      #add to list
      gene.family.sizes.n5 <- c(gene.family.sizes.n5, gene.family.size)
      #make linear model
      local.linear.model <- lm(standardized.expression.pattern.diversity~gene.cluster.phylogenetic.diversity, data = local.data)
      #print(local.linear.model)
      local.linear.model.summary <- summary(local.linear.model)
      local.observed.slope <- local.linear.model.summary$coefficients[2]
      local.observed.intercept <- local.linear.model.summary$coefficients[1]
      local.observed.r2 <- local.linear.model.summary$r.squared
      
      #add to local values
      local.slopes[[gene.family.size]] <- local.observed.slope
      local.r2[[gene.family.size]] <- local.observed.r2
      
      #add gene family sizes to null list
      local.null.slope.values[[gene.family.size]] <- c()
      local.null.intercept.values[[gene.family.size]] <- c()
      local.null.r2.values[[gene.family.size]] <- c()
    }
  }
  
  for (i in 1:iterations) {
    print(paste("Working on iteration: ", i, sep=""))
    #generate random gene clusters
    random.clusters <- random.gene.clusters(clusters.file, total.expression.data)
    
    #global analysis
    
    #calculate expression pattern diversity
    global.expression.pattern.diversity <- compute.expression.pattern.diversity(random.clusters, total.expression.data)
    #return(list("gene.family.data" = gene.family.data, "global.expression.pattern.diversity" = global.expression.pattern.diversity))
    #filter data
    #unique_clusters <- unique(gene.family.data$cluster)
    #global.expression.pattern.diversity <- global.expression.pattern.diversity[global.expression.pattern.diversity$cluster %in% unique_clusters, ]
    #combine data
    
    global.expression.pattern.diversity <- cbind(global.expression.pattern.diversity, data.frame(gene.family.data$gene.cluster.phylogenetic.diversity))
    global.expression.pattern.diversity <- cbind(global.expression.pattern.diversity, data.frame(gene.family.data$n.genes))
    #fit linear model to standardized gene
    global.expression.pattern.diversity$gene.family.data.gene.cluster.phylogenetic.diversity <- as.numeric(global.expression.pattern.diversity$gene.family.data.gene.cluster.phylogenetic.diversity)
    global.expression.pattern.diversity$standardized.expression.pattern.diversity <- as.numeric(global.expression.pattern.diversity$standardized.expression.pattern.diversity)
    random.linear.model <- lm(standardized.expression.pattern.diversity~gene.family.data.gene.cluster.phylogenetic.diversity, data = global.expression.pattern.diversity)
    #get values
    random.linear.model.summary <- summary(random.linear.model)
    random.global.slope <- random.linear.model.summary$coefficients[2]
    random.global.intercept <- random.linear.model.summary$coefficients[1]
    random.global.r2 <- random.linear.model.summary$r.squared
    
    #add values to output
    global.null.slope.values <- c(global.null.slope.values, random.global.slope)
    global.null.r2.values <- c(global.null.r2.values, random.global.r2)
    global.null.intercept.values <- c(global.null.intercept.values, random.global.intercept)
    
    #local analysis
    
    #iterate through each gene family size that has appropriate sample size
    for (gene.family.size in gene.family.sizes.n5) {
      #subset data
      random.local.data <- global.expression.pattern.diversity[global.expression.pattern.diversity$gene.family.data.n.genes == gene.family.size, ]
      #fit linear model
      random.local.linear.model <- lm(standardized.expression.pattern.diversity~gene.family.data.gene.cluster.phylogenetic.diversity, data = random.local.data)
      #get values
      random.local.linear.model.summary <- summary(random.local.linear.model)
      random.local.slope <- random.local.linear.model.summary$coefficients[2]
      random.local.intercept <- random.local.linear.model.summary$coefficients[1]
      random.local.r2 <- random.local.linear.model.summary$r.squared
      
      #add values to null lists
      local.null.slope.values[[gene.family.size]] <- c(local.null.slope.values[[gene.family.size]], random.local.slope)
      local.null.r2.values[[gene.family.size]] <- c(local.null.r2.values[[gene.family.size]], random.local.r2)
      local.null.intercept.values[[gene.family.size]] <- c(local.null.intercept.values[[gene.family.size]], random.local.intercept)
    }
  }
  #assemble output
  null.model.results <- list("global.slope" = global.observed.slope,
                             "global.intercept" = global.observed.intercept,
                             "global.null.slope.distribution" = global.null.slope.values,
                             "global.null.intercept.distribution" = global.null.intercept.values,
                             "global.r2" = global.observed.r2,
                             "global.null.r2.distribution" = global.null.r2.values,
                             "local.slopes" = local.slopes,
                             "local.r2" = local.r2,
                             "local.null.slope.distributions" = local.null.slope.values,
                             "local.null.intercept.distributions" = local.null.intercept.values,
                             "local.null.r2.distributions" = local.null.r2.values)
  #return(null.model.results)
  return(null.model.results)
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

write.csv(median.values, '/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/total_median_expression.csv', row.names=FALSE)

#convert gene ids to protein ids
median.values$genes <- sapply(median.values$genes, get.prot.id)

#Running
gene.clusters.file <- "/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/psiblast_id-20_e-neg-5_cov-0.7_sequence_clusters.json"
gene.clusters.data <- rjson::fromJSON(file=gene.clusters.file)
gene.expression.diversity <- compute.expression.pattern.diversity(gene.clusters.data, median.values)

phylogenies.path <- '/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/inferred_phylogenies_distant'
file.extension <- '_alignment.fasta.treefile'
phylo.data <- compute.gene.cluster.phylogenetic.diversity(phylogenies.path, file.extension, gene.expression.diversity)

write.csv(phylo.data, "/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/gene_family_diversity_vs_expression_diversity_global_data_distant.csv", row.names=FALSE)

#null testing
#remove rows with only 0s from total data frame
expression.data.no0 <- median.values[rowSums(median.values != 0) > 0, ]



#run null test
null.data <- null.test(total.expression.data = expression.data.no0, 
                        gene.family.data = phylo.data, 
                        clusters.file = gene.clusters.data, 
                        iterations = 1000)

#write results
json.data <- rjson::toJSON(null.data)
write(json.data, file='/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/null_model_results_distant.json')


