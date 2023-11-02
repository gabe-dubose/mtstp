#!/usr/bin/env Rscript

library("edgeR")

#load metadata
sample_metadata <- read.csv('data/mtstp_analysis_metadata.tsv', sep="\t")

#load raw counts data
raw_counts_matrix <- (read.csv('data/counts_tables/dpl_raw_counts_kallisto.csv'))

uninfected_metadata <- sample_metadata[sample_metadata$infection.status != 'infected',]
uninfected_data <- raw_counts_matrix[raw_counts_matrix$X %in% uninfected_metadata$sample.id,]

#function to parse and arrange data
parse.data <- function(from, to, data, metadata, plant, comp1, comp2) {
  
  #separate plant data
  plant_metadata <- metadata[metadata$plant == plant,]
  data <- data[data$X %in% plant_metadata$sample.id,]
  
  #separate data
  categories = c(from, to)
  stages_metadata <- plant_metadata[plant_metadata$developmental.stage %in% categories,]
  stages_data <- data[data$X %in% stages_metadata$sample.id,]
  
  #convert from dataframe to correct format
  stages_data_t <- t(stages_data)
  colnames(stages_data_t) <- stages_data_t[1,]
  stages_data <- stages_data_t[-1,]
  stages_data <- matrix(as.numeric(stages_data), ncol = ncol(stages_data))
  colnames(stages_data) <- stages_data_t[1,]
  rownames(stages_data) <- rownames(stages_data_t)[-1]
  
  #define comparisons
  comparisons <- c()
  for (sample in colnames(stages_data)) {
    if (grepl(comp1, sample) == TRUE) {
      comparisons <- c(comparisons, from)
    }
    else if (grepl(comp2, sample) == TRUE) {
      comparisons <- c(comparisons, to)
    }
  }
  
  #assemble output
  formatted_data <- list("matrix" = stages_data, "comparisons" = comparisons)
  
  return(formatted_data)
}

#function to get and write edgeR data
get.edgeR.dges <- function(data, outfile.handle, plant) {
  
  #perform DGE analysis
  dge.obj <- DGEList(counts = data$matrix, group = factor(data$comparisons))
  #filter low quantity genes
  dge.obj.keep <- filterByExpr(y = dge.obj)
  dge.obj.filtered <- dge.obj[dge.obj.keep, , keep.lib.sizes=FALSE]
  #perform normalization
  dge.obj.filtered.normalized <- calcNormFactors(object = dge.obj.filtered)
  #estimate dispersion
  dge.dispersion <- estimateDisp(y = dge.obj.filtered.normalized)
  #test for differential expression
  dge <- exactTest(object = dge.dispersion)
  #perform fdr correction
  dge.fdr <- topTags(object = dge, n = "Inf")
  #write file
  write.csv(dge.fdr, paste(outfile.handle, plant, sep="_"))
  
}

#3-5
third.to.fifth.inc.data <- parse.data('third-instar', 'fifth-instar', uninfected_data, sample_metadata, 'incarnata', "3", '5')
get.edgeR.dges(third.to.fifth.inc.data, 'data/edgeR/across_development/third_to_fifth_inc', 'incarnata')
third.to.fifth.cur.data <- parse.data('third-instar', 'fifth-instar', uninfected_data, sample_metadata, 'curassavica', "3", '5')
get.edgeR.dges(third.to.fifth.cur.data, 'data/edgeR/across_development/third_to_fifth_cur', 'curassavica')

#5-E
fifth.to.early.inc.data <- parse.data('fifth-instar', 'early-pupa', uninfected_data, sample_metadata, 'incarnata', "5", 'E')
get.edgeR.dges(fifth.to.early.inc.data, 'data/edgeR/across_development/fifth_to_early_inc', 'incarnata')
fifth.to.early.cur.data <- parse.data('fifth-instar', 'early-pupa', uninfected_data, sample_metadata, 'curassavica', "5", 'E')
get.edgeR.dges(fifth.to.early.cur.data, 'data/edgeR/across_development/fifth_to_early_cur', 'curassavica')

#E-L
early.to.late.inc.data <- parse.data('early-pupa', 'late-pupa', uninfected_data, sample_metadata, 'incarnata', "E", 'L')
get.edgeR.dges(early.to.late.inc.data, 'data/edgeR/across_development/early_to_late_inc', 'incarnata')
early.to.late.cur.data <- parse.data('early-pupa', 'late-pupa', uninfected_data, sample_metadata, 'curassavica', "E", 'L')
get.edgeR.dges(early.to.late.cur.data, 'data/edgeR/across_development/early_to_late_cur', 'curassavica')

#L-A
late.to.adult.inc.data <- parse.data('late-pupa', 'adult', uninfected_data, sample_metadata, 'incarnata', "L", 'A')
get.edgeR.dges(late.to.adult.inc.data, 'data/edgeR/across_development/late_to_adult_inc', 'incarnata')
late.to.adult.cur.data <- parse.data('late-pupa', 'adult', uninfected_data, sample_metadata, 'curassavica', "L", 'A')
get.edgeR.dges(late.to.adult.cur.data, 'data/edgeR/across_development/late_to_adult_cur', 'curassavica')
