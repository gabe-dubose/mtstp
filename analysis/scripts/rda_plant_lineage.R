#!/usr/bin/env Rscript

#load distance matrix
distances <- read.csv('data/distance_matricies/euclidian.csv')
#rename x to sampe.id
colnames(distances)[1] <- 'sample.id'
metadata <- read.csv('data/mtstp_analysis_metadata.tsv', sep='\t')
#remove infected data
infected.ids <- metadata[metadata$infection.status == 'infected',]$sample.id
distances <- distances[,!names(distances) %in% infected.ids]
distances <- distances[!distances$sample.id %in% infected.ids,]

#add metadata categories
genetic.background <- c()
plant.species <- c()
developmental.stages <- c()

for (sample in distances$sample.id) {
  plant <- metadata[metadata$sample.id == sample,]$plant
  lineage <- metadata[metadata$sample.id == sample,]$lineage
  stage <- metadata[metadata$sample.id == sample,]$developmental.stage
  
  genetic.background <- c(genetic.background, lineage)
  plant.species <- c(plant.species, plant)
  developmental.stages <- c(developmental.stages, stage)
}

#add metadata to distance
distances <- cbind(distances, genetic.background)
distances <- cbind(distances, plant.species)
distances <- cbind(distances, developmental.stages)

#remove infected data
infected.ids <- metadata[metadata$infection.status == 'infected',]$sample.id
distances <- distances[,!names(distances) %in% infected.ids]
distances <- distances[!distances$sample.id %in% infected.ids,]
metadata <- metadata[metadata$infection.status != 'infected',]

#partition data by developmental stage
non.third.ids <- metadata[metadata$developmental.stage != 'third-instar',]$sample.id
third.instar.data <- distances[,!names(distances) %in% non.third.ids]
third.instar.data <- third.instar.data[!third.instar.data$sample.id %in% non.third.ids,]

non.fifth.ids <- metadata[metadata$developmental.stage != 'fifth-instar',]$sample.id
fifth.instar.data <- distances[,!names(distances) %in% non.fifth.ids]
fifth.instar.data <- fifth.instar.data[!fifth.instar.data$sample.id %in% non.fifth.ids,]

non.early.ids <- metadata[metadata$developmental.stage != 'early-pupa',]$sample.id
early.pupa.data <- distances[,!names(distances) %in% non.early.ids]
early.pupa.data <- early.pupa.data[!early.pupa.data$sample.id %in% non.early.ids,]

non.late.ids <- metadata[metadata$developmental.stage != 'late-pupa',]$sample.id
late.pupa.data <- distances[,!names(distances) %in% non.late.ids]
late.pupa.data <- late.pupa.data[!late.pupa.data$sample.id %in% non.late.ids,]

non.adult.ids <- metadata[metadata$developmental.stage != 'adult',]$sample.id
adult.data <- distances[,!names(distances) %in% non.adult.ids]
adult.data <- adult.data[!adult.data$sample.id %in% non.adult.ids,]

#initialize dataframe to store data
rda.data <- data.frame()


rda.data.columns <- c("developmental.stage", 
                      "explained.by.genetic.background", 
                      "explained.by.plant",
                      "p.genetic.background",
                      "p.plant")

#run dbRDA for third instar
#partition actual distance matrix
third.distances <- third.instar.data[2:11]
#make model
third.dbRDA.model <- vegan::rda(formula = third.distances~genetic.background+plant.species, data=third.instar.data)
#get model summary
third.dbRDA.model.summary <- summary(third.dbRDA.model)

#perform anova and identify significant 
third.dbRDA.anova <- anova(third.dbRDA.model, by="terms", permu=999)

#partition variance
third.variance.partitioned <- vegan::varpart(third.distances, ~genetic.background, ~plant.species, data=third.instar.data)
third.variance.partitioned.summary <- summary(third.variance.partitioned)
#get variance explained 
third.explained <- third.variance.partitioned.summary$uniqpart

#run dbRDA for fifth instar
#partition actual distance matrix
fifth.distances <- fifth.instar.data[2:11]
#make model
fifth.dbRDA.model <- vegan::rda(formula = fifth.distances~genetic.background+plant.species, data=fifth.instar.data)
#get model summary
fifth.dbRDA.model.summary <- summary(fifth.dbRDA.model)

#perform anova and identify significant 
fifth.dbRDA.anova <- anova(third.dbRDA.model, by="terms", permu=999)

#partition variance
fifth.variance.partitioned <- vegan::varpart(fifth.distances, ~genetic.background, ~plant.species, data=fifth.instar.data)
fifth.variance.partitioned.summary <- summary(fifth.variance.partitioned)
#get variance explained 
fifth.explained <- fifth.variance.partitioned.summary$uniqpart

#run dbRDA for early pupa
#partition actual distance matrix
early.distances <- early.pupa.data[2:10]
#make model
early.dbRDA.model <- vegan::rda(formula = early.distances~genetic.background+plant.species, data=early.pupa.data)
#get model summary
early.dbRDA.model.summary <- summary(early.dbRDA.model)

#perform anova and identify significant 
early.dbRDA.anova <- anova(early.dbRDA.model, by="terms", permu=999)

#partition variance
early.variance.partitioned <- vegan::varpart(early.distances, ~genetic.background, ~plant.species, data=early.pupa.data)
early.variance.partitioned.summary <- summary(early.variance.partitioned)
#get variance explained 
early.explained <- early.variance.partitioned.summary$uniqpart

#run dbRDA for late pupa
#partition actual distance matrix
late.distances <- late.pupa.data[2:11]
#make model
late.dbRDA.model <- vegan::rda(formula = late.distances~genetic.background+plant.species, data=late.pupa.data)
#get model summary
late.dbRDA.model.summary <- summary(late.dbRDA.model)

#perform anova and identify significant 
late.dbRDA.anova <- anova(late.dbRDA.model, by="terms", permu=999)

#partition variance
late.variance.partitioned <- vegan::varpart(late.distances, ~genetic.background, ~plant.species, data=late.pupa.data)
late.variance.partitioned.summary <- summary(late.variance.partitioned)
#get variance explained 
late.explained <- late.variance.partitioned.summary$uniqpart

#run dbRDA for adult
#partition actual distance matrix
adult.distances <- adult.data[2:11]
#make model
adult.dbRDA.model <- vegan::rda(formula = adult.distances~genetic.background+plant.species, data=adult.data)
#get model summary
adult.dbRDA.model.summary <- summary(adult.dbRDA.model)

#perform anova and identify significant 
adult.dbRDA.anova <- anova(adult.dbRDA.model, by="terms", permu=999)

#partition variance
adult.variance.partitioned <- vegan::varpart(adult.distances, ~genetic.background, ~plant.species, data=adult.data)
adult.variance.partitioned.summary <- summary(adult.variance.partitioned)
#get variance explained 
adult.explained <- adult.variance.partitioned.summary$uniqpart

#assemble data
third.rda.results <- c("third-instar", third.explained["X1"], third.explained["X2"], third.dbRDA.anova$`Pr(>F)`[1], third.dbRDA.anova$`Pr(>F)`[2])
fifth.rda.results <- c("fifth-instar", fifth.explained["X1"], fifth.explained["X2"], fifth.dbRDA.anova$`Pr(>F)`[1], fifth.dbRDA.anova$`Pr(>F)`[2])
early.rda.results <- c("early-pupa", early.explained["X1"], early.explained["X2"], early.dbRDA.anova$`Pr(>F)`[1], early.dbRDA.anova$`Pr(>F)`[2])
late.rda.results <- c("late-pupa", late.explained["X1"], late.explained["X2"], late.dbRDA.anova$`Pr(>F)`[1], late.dbRDA.anova$`Pr(>F)`[2])
adult.rda.results <- c("adult", adult.explained["X1"], adult.explained["X2"], adult.dbRDA.anova$`Pr(>F)`[1], adult.dbRDA.anova$`Pr(>F)`[2])

rda.data <- rbind(rda.data, third.rda.results)
rda.data <- rbind(rda.data, fifth.rda.results)
rda.data <- rbind(rda.data, early.rda.results)
rda.data <- rbind(rda.data, late.rda.results)
rda.data <- rbind(rda.data, adult.rda.results)

rda.data.columns <- c("developmental.stage", 
                      "explained.by.genetic.background", 
                      "explained.by.plant",
                      "p.genetic.background",
                      "p.plant")
colnames(rda.data) <- rda.data.columns

write.csv(rda.data, 'data/variance_paritioned.csv', row.names=FALSE)