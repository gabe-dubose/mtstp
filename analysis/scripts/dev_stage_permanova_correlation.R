library(vegan)
library(tidyverse)
#load data
data <- read.csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/dpl_tpm_counts_kallisto.csv')
colnames(data)[1] <- "sample.id"

#load metadata
metadata <- read.csv('/home/gabe/Desktop/mtstp/data/experiment_metadata/mtstp_analysis_metadata.tsv', sep='\t')
metadata <- subset(metadata, infection.status != "infected")
colnames(metadata)
#remove infected

test <- metadata %>%
  left_join(data, by = c("sample.id"))

trt <- test %>%
  select(c("sample.id","developmental.stage", "plant"))

data2 <- test %>%
  select(-c("sample.id", "plant", "infection.status", "developmental.stage", "lineage"))


#calculat distance matrix using pearson distances
data3 <- t(data2)
pearson.correlation.matrix <- cor(data3, method="pearson")
pearson.correlation.matrix <- (1 - pearson.correlation.matrix)/2
pearson.correlation.matrix <- as.dist(pearson.correlation.matrix)

ad.test <- adonis2(pearson.correlation.matrix ~ developmental.stage+plant, data = trt)
ad.test

#This code was retrieved from https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- eval(x1[[2]], environment(x1), globalenv())
  environment(x1) <- environment()
  # extract factors on right hand side of formula
  rhs <- x1[[3]]
  # create model.frame matrix
  x1[[2]] <- NULL
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)
  
  # create unique pairwise combination of factors
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors
  for(elem in 1:ncol(co)){
    
    #reduce model elements
    if(inherits(eval(lhs),'dist')){
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))
    }else{
      xnew <- as.formula(paste('xred' ,
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names
  class(res) <- c("pwadstrata", "list")
  return(res)
}


### Method summary
summary.pwadstrata = function(object, ...) {
  cat("Result of pairwise.adonis2:\n")
  cat("\n")
  print(object[1], ...)
  cat("\n")
  
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

#run pairwise permanova
ad.pairwise <- pairwise.adonis2(pearson.correlation.matrix ~ developmental.stage+plant, data = trt)
ad.pairwise$`third-instar_vs_fifth-instar`
ad.pairwise$`fifth-instar_vs_early-pupa`
ad.pairwise$`early-pupa_vs_late-pupa`
ad.pairwise$`adult_vs_late-pupa`

ad.pairwise$`third-instar_vs_adult`
ad.pairwise$`fifth-instar_vs_adult`
ad.pairwise$`third-instar_vs_early-pupa`
