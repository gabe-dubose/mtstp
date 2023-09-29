library("edgeR")

#load metadata
sample_metadata <- read.csv('data/mtstp_analysis_metadata.tsv', sep="\t")

#load raw counts data
raw_counts_matrix <- (read.csv('data/counts_tables/dpl_raw_counts_kallisto.csv'))

#function to parse and arrange data
parse.data <- function(developmental.stage) {
  #get instar specific metadata
  instar_metadata <- sample_metadata[sample_metadata$developmental.stage == developmental.stage,]
  #get instar specific data
  instar_data <- raw_counts_matrix[raw_counts_matrix$X %in% instar_metadata$sample.id,]
  
  #separate incarnata data
  incarnata_metadata <- instar_metadata[instar_metadata$plant == 'incarnata',]
  incarnata_data <- instar_data[instar_data$X %in% incarnata_metadata$sample.id,]
  
  #separate curassavica data
  curassavica_metadata <- instar_metadata[instar_metadata$plant == 'curassavica',]
  curassavica_data <- instar_data[instar_data$X %in% curassavica_metadata$sample.id,]
  
  #convert incarnata dataframe to correct format
  incarnata_data_t <- t(incarnata_data)
  colnames(incarnata_data_t) <- incarnata_data_t[1,]
  incarnata_data <- incarnata_data_t[-1,]
  incarnata_data <- matrix(as.numeric(incarnata_data), ncol = ncol(incarnata_data))
  colnames(incarnata_data) <- incarnata_data_t[1,]
  rownames(incarnata_data) <- rownames(incarnata_data_t)[-1]
  
  #convert curassavica dataframe to correct format
  curassavica_data_t <- t(curassavica_data)
  colnames(curassavica_data_t) <- curassavica_data_t[1,]
  curassavica_data <- curassavica_data_t[-1,]
  curassavica_data <- matrix(as.numeric(curassavica_data), ncol = ncol(curassavica_data))
  colnames(curassavica_data) <- curassavica_data_t[1,]
  rownames(curassavica_data) <- rownames(curassavica_data_t)[-1]
  
  #define incarnata comparisons
  incarnata.comparisons <- c()
  for (sample in colnames(incarnata_data)) {
    if (grepl("ii", sample) == TRUE) {
      incarnata.comparisons <- c(incarnata.comparisons, "infected")
    }
    if (grepl("iu", sample) == TRUE) {
      incarnata.comparisons <- c(incarnata.comparisons, "uninfected")
    }
  }
  
  #define curassavica comparisons
  curassavica.comparisons <- c()
  for (sample in colnames(curassavica_data)) {
    if (grepl("ci", sample) == TRUE) {
      curassavica.comparisons <- c(curassavica.comparisons, "infected")
    }
    if (grepl("cu", sample) == TRUE) {
      curassavica.comparisons <- c(curassavica.comparisons, "uninfected")
    }
  }
  
  #assemble output
  formatted_data <- list("incarnata.counts.matrix" = incarnata_data, "incarnata.comparisons" = incarnata.comparisons, 
                      "curassavica.counts.matrix" = curassavica_data, "curassavica.comparisons" = curassavica.comparisons)
  
  return(formatted_data)
}

#function to get and write edgeR data
get.edgeR.dges <- function(developmental.stage, outfile.handle, data) {
  
  #perform DGE analysis for curassavica
  #define dgelist object for curassavica
  curassavica.dge.obj <- DGEList(counts = data$curassavica.counts.matrix, group = factor(data$curassavica.comparisons))
  #filter low quantity genes
  curassavica.dge.obj.keep <- filterByExpr(y = curassavica.dge.obj)
  curassavica.dge.obj.filtered <- curassavica.dge.obj[curassavica.dge.obj.keep, , keep.lib.sizes=FALSE]
  #perform normalization
  curassavica.dge.obj.filtered.normalized <- calcNormFactors(object = curassavica.dge.obj.filtered)
  #estimate dispersion
  curassavica.dge.dispersion <- estimateDisp(y = curassavica.dge.obj.filtered.normalized)
  #test for differential expression
  curassavica.dge <- exactTest(object = curassavica.dge.dispersion)
  #perform fdr correction
  curassavica.dge.fdr <- topTags(object = curassavica.dge, n = "Inf")
  #write file
  write.csv(curassavica.dge.fdr, paste(outfile.handle, "curassavica", sep="_"))
  
  #perform DGE analysis for incarnata
  #define dgelist object for incarnata
  incarnata.dge.obj <- DGEList(counts = data$incarnata.counts.matrix, group = factor(data$incarnata.comparisons))
  #filter low quantity genes
  incarnata.dge.obj.keep <- filterByExpr(y = incarnata.dge.obj)
  incarnata.dge.obj.filtered <- incarnata.dge.obj[incarnata.dge.obj.keep, , keep.lib.sizes=FALSE]
  #perform normalization
  incarnata.dge.obj.filtered.normalized <- calcNormFactors(object = incarnata.dge.obj.filtered)
  #estimate dispersion
  incarnata.dge.dispersion <- estimateDisp(y = incarnata.dge.obj.filtered.normalized)
  #test for differential expression
  incarnata.dge <- exactTest(object = incarnata.dge.dispersion)
  #perform fdr correction
  incarnata.dge.fdr <- topTags(object = incarnata.dge, n = "Inf")
  #write file
  write.csv(incarnata.dge.fdr, paste(outfile.handle, "incarnata", sep="_"))
}

#get third instar dges
third.instar.data <- parse.data("third-instar")
get.edgeR.dges('third-instar', '/home/gabe/Desktop/mtstp/analysis/data/edgeR_data/third_instar_inf_vs_uninf', third.instar.data)

#get fifth instar dges
fifth.instar.data <- parse.data("fifth-instar")
get.edgeR.dges('fifth-instar', '/home/gabe/Desktop/mtstp/analysis/data/edgeR_data/fifth_instar_inf_vs_uninf', fifth.instar.data)

#get early pupa dges
early.pupa.data <- parse.data("early-pupa")
get.edgeR.dges('early-pupa', '/home/gabe/Desktop/mtstp/analysis/data/edgeR_data/early_pupa_inf_vs_uninf', early.pupa.data)

#get late pupa dges
late.pupa.data <- parse.data("late-pupa")
get.edgeR.dges('late-pupa', '/home/gabe/Desktop/mtstp/analysis/data/edgeR_data/late_pupa_inf_vs_uninf', late.pupa.data)

#get adult dges
adult.data <- parse.data("adult")
get.edgeR.dges('adult', '/home/gabe/Desktop/mtstp/analysis/data/edgeR_data/adult_inf_vs_uninf', adult.data)
