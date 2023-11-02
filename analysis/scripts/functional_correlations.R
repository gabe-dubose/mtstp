#!/usr/bin/env Rscript

metabolism <- read.csv('data/function_tables/metabolism_nonredundant.csv')
organismal.systems <- read.csv('data/function_tables/organismal_systems_nonredundant.csv')
genetic.information.processing <- read.csv('data/function_tables/genetic_information_processing_nonredundant.csv')
cellular.processes <- read.csv('data/function_tables/cellular_processes_nonredundant.csv')

drop <- c("plant", "developmental_stage")
metabolism = metabolism[,!(names(metabolism) %in% drop)]
organismal.systems = organismal.systems[,!(names(organismal.systems) %in% drop)]
genetic.information.processing = genetic.information.processing[,!(names(genetic.information.processing) %in% drop)]

data <- merge(metabolism, organismal.systems, by="sample.id")
data <- merge(data, genetic.information.processing, by="sample.id")
data <- merge(data, cellular.processes, by="sample.id")

df <- data[,!(names(data) %in% drop)]
df <- df[,-1]
cor.data <- cor(df)
cor.data <- data.frame(cor.data)

write.csv(cor.data, 'data/function_tables/functional_correlations.csv')
