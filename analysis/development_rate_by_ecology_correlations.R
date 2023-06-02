#load metadata
data <- read.csv("../metadata/mtstp_sample_metadata.tsv", sep="\t")

#Correlations between infection+plant with time from hatching to pupation
#a quick boxplot for visualization
boxplot(days.hatch.to.pupation~infection.status*plant.diet, data=data)
#This looks kind of weird. Might be better to visualize with histogram
#run an ANOVA to look for significant differences
days.hatch.to.pupation.anova <- aov(days.hatch.to.pupation~infection.status*plant.diet, data=data)
summary(days.hatch.to.pupation.anova)
#looks like there are differences by plant diet (p=1.16e-08), but not by infection status (p=0.110) or plant by infection (p=0.858)
#and now for a post-hoc test
days.hatch.to.pupation.tukey <- TukeyHSD(days.hatch.to.pupation.anova)
days.hatch.to.pupation.tukey

#Correlations between infection*plant with time from pupation to eclosion
#a quick boxplot for visualization
boxplot(days.pupation.to.eclosion~infection.status*plant.diet, data=data)
#doesn't seem like many differences
#but I guess I'll look with an ANOVA
days.pupation.to.eclosion.anova <- aov(days.pupation.to.eclosion~infection.status*plant.diet, data=data)
summary(days.pupation.to.eclosion.anova)
#looks like there might be difference by plant by diet (p=0.0186), but not infection status (p=0.3502) or plant diet (p=0.8223)
#and now for a post-hoc test
days.pupation.to.eclosion.tukey <- TukeyHSD(days.pupation.to.eclosion.anova)
days.pupation.to.eclosion.tukey
#looks like there might be a difference between infect and uninfected on incarnata if there was more power (p=0.0884416), but nothing else was even close

#Correlations between infection*plant with time from hatching to eclosion
#a quick boxplot for visualization
boxplot(days.hatch.to.eclosion~infection.status*plant.diet, data=data)
#still looks a little funny. I guess there just isn't much variation in this specific data
#run an ANOVA to look for significant differences
days.hatch.to.eclosion.anova <- aov(days.hatch.to.eclosion~infection.status*plant.diet, data=data)
summary(days.hatch.to.eclosion.anova)
#looks like plant diet and infection status by plant diet might be significant with more statistical power (p=0.0685 and p=0.0778, respectively)
#since nothing was actually significant, I'm not doing a post-hoc analysis
