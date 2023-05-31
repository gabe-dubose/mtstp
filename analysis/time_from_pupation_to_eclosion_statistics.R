#load metadata
data <- read.csv("../metadata/mtstp_sample_metadata.tsv", sep="\t")

#make a quick boxplot for visualization
boxplot(time.from.pupation.to.eclosion.days~infection.status*plant, data=data)

#run an ANOVA test
anova.test <- aov(time.from.pupation.to.eclosion.days~infection.status*plant, data=data)
summary(anova.test)

#use Tukey test for post-hoc analysis
tukey.test <- TukeyHSD(anova.test)
tukey.test

