library("ROntoTools")
library("graph")
library("KEGGREST")

dpl_pathways <- keggPathwayGraphs("dpl", updateCache = TRUE, verbose = FALSE)

res <- keggList("dpl")
