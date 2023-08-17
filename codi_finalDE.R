###############################################################
#### DIFFERENCIAL EXPRESION ANALYSIS WITH edgeR & limma-trend
## TFM KAREN JUNE 2023

### Tots els comentaris que he posat poden ser una guia per tu per saber el que s'ha fet
#Jo ho editaria una mica, deixant alguns comentaris on s'expliqui que fa el codi
#Al TFM posa la informació més rellevant, es a dir quins paràmetres i paquet de R s'han fet servir per filtrar, normalitzar...
#Quins paràmetres i pq s'han fet servir per el DE...
#És a dir, aquí només deixa la info referent a que fa cada linea de codi i al TFM dona una explicació de perquè has fet una cosa o l'altre
## Una altra cosa! Alguna passa d'aquest codi l'he treta de chatGPT, per al codi no passa res però canvia tots els comentaris i posa'ls amb les teves paraules
#No sigui cosa que tenguin algun programa que ho detecti i et diguin alguna cosa!

library(edgeR)
library(limma)

rm(list=ls())

#To visualize a DGEList object just type the name, do not use str() or summary() far less informative

#Normalization info:
#When you first create a DGEList object, all the normalization factors are initially set to 1 by default. 
#This is equivalent to no normalization. Once you use calcNormFactors, the normalization factors will be set appropriately.

### CAREFUL! MANY DOWNSTREAM ANALYSIS RELY ON ORDER OF SAMPLES. MAKE SURE THAT ALL THIS FILES HAVE THE SAMPLES SORTED IN THE SAME ORDER!
#PREFERABLY IN ALPHABETICAL ORDER OR SORTED FROM SMALLER TO LARGEST.

#No deixis els meus paths a l'arxiu que posis al TFM...
counts <- read.table("/Users/mbarcelo/Dropbox/KarenTFM/Codi_R/RSEM_counts.txt", header = TRUE, row.names = 1)
summary(counts)

#Aquests dos arxius que hi ha a continuació es podrien condensar amb un que tengui tota la info... 
#Ho pots canviar si vols, pero pensa que si ho canvies també ho has de canviar a la resta de script on es facin servir les variables que fan referència a aquests arxius
info_samples <- read.table("/Users/mbarcelo/Dropbox/KarenTFM/Codi_R/samples.txt", header = TRUE)
summary(info_samples)

samples <- read.table("/Users/mbarcelo/Dropbox/KarenTFM/Codi_R/samples_edgeR.txt", header = TRUE)
samples$Sex <- as.factor(samples$Sex)

summary(samples)

### Use edgeR to create a DGEList, do the filtering and normalization.

y <- DGEList(counts = counts, group = NULL, lib.size = info_samples$mapped, samples = samples ) 

### Plot log of counts to see how the distribution looks like
logCPM <- cpm(y, log= TRUE)
plotDensities(logCPM, legend = FALSE,  main = "Before Filtering")
#Genes before filter: 26690


######## FILTER
## Low read counts. 
#Use function filterByExpr to go over count matrix and decide what are low counts
keep <- filterByExpr(y)
table(keep)
#Create a new object y that will have only the genes that are informative and do not have very low counts
y <- y[keep, , keep.lib.sizes=FALSE]
#Genes in matrix after filter: 16278

### Plot log of counts to see how the distribution looks like after filter
logCPM <- cpm(y, log= TRUE)
plotDensities(logCPM, legend = FALSE,  main = "After Filtering")

#####Plot to see effects of normalization:
#Library sizes per sample

par(mfrow=c(1,2))
boxplot(logCPM)
title(main="Unnormalised data", ylab="Log-cpm")

#As suspected the libraries are fairly similar, it would not be a big deal if we did not normalize...

##### Calculate the normalization Factors
y <- calcNormFactors(y)
y
plotMDS(y)

##Plot again with new normalization factors:
logCPM <- cpm(y, log= TRUE)
boxplot(logCPM)
title(main="Normalised data", ylab="Log-cpm")

##Apply normalization factors to the data and get the normalized counts:

normalized_counts <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)

# logCPM <- cpm(dge, log=TRUE, prior.count=3)
# normalized_counts

#Problem, limma works with count matrices.. need to create a matrix:

normalized_matrix <- as.matrix(normalized_counts)

summary(normalized_matrix)
## Now we need to define our experimental design to calculate dif. expr

#First we define our variables, activity onset (AO) and sex (Sex)
AO <- samples$A_all
Sex <- samples$Sex


#With this design (~ 0 + AO*Sex) we specify that we are not interested in the intercept.
# When you use 0 + in the formula, it indicates that the design matrix should not include an intercept term. 
# Instead, it will include separate coefficients for each level of the categorical variable, representing the differences 
# from the overall mean expression level. This can be useful when you want to compare the expression levels 
# between different levels of the categorical variable without assuming a baseline or reference level.
# When we remove the intercept (model above) we are forcing the regression line to cross the y axis at 0... 
# The fitting will be way worse! What it was written above it makes sense for categorical variables, not recommended for continuous...
# In general, we suggest the inclusion of an intercept term for modelling explanatory variables that are 
# covariates (continuous) since it provides a more flexible fit to the data points. A model without an intercept 
# term would only be recommended in cases where there is a strong biological reason why a zero covariate should be 
# associated with a zero expression value, and such contexts are rare in gene expression modelling.


#We add the terms Activity Onset and Sex as an interaction, in this way we will control for the
#differences in AO related to sex.

design <- model.matrix(~ AO*Sex)


#Check with voom plot if filtering was good enough 
#We can ingore this for now.... No need to mention in the TFM
v <- voom(y, design, plot=TRUE)


#Fitting the linear model using the normalized counts and the experimental design:

fit <- lmFit(normalized_matrix, design)
#Ensuring correction with eBayes usling limma-trend
fit <- eBayes(fit, trend=TRUE)
fit$coefficients
#This p-value will be already corrected using the eBayes function. 
fit$p.value
str(fit)
#However, if we check the distribution of the coefficients and p-values we can see that we are finding
#as many statistically significant genes as expected by random chance. This is because we are doing many
#tests (one for each gene! 16K). 
hist(fit$coefficients[,2])
hist(fit$p.value[,2]) #, xlim = c(0, 0.05))


#We need to do a multiple test adjustment of p-values.
#In the literature the BH correction is found to be the most robust.

#I finally understood what the term coef was! 
#In this case you are suposed to use the coeficients of the variable you want to adjust the p-value!
#coef = 2 corresponds to the coefficients of the second variable AO (the first one is the intercept)
#coef = 3 corresponds to the coefficients of the third variable Sex
#coef = 4 corresponds to the coefficients of the fourth variable AO*Sex

#Results with adjusted p-values for AO
resultsAO <- topTable(fit, coef = 2, number = 10, adjust.method = "BH")

#Results with adjusted p-values for Sex
resultsSex <- topTable(fit, coef = 3, number = 10, adjust.method = "BH")

#Results with adjusted p-values for interaction AO*Sex
resultsAO_Sex <- topTable(fit, coef = 4, number = 10, adjust.method = "BH")

file_path <- "/Users/mbarcelo/Dropbox/KarenTFM/Codi_R/AO_DE.txt"
# Save the table to a text file
write.table(resultsAO, file = file_path, sep = "\t", quote = FALSE, row.names = TRUE)

file_path <- "/Users/mbarcelo/Dropbox/KarenTFM/Codi_R/Sex_DE.txt"
# Save the table to a text file
write.table(resultsSex, file = file_path, sep = "\t", quote = FALSE, row.names = TRUE)

file_path <- "/Users/mbarcelo/Dropbox/KarenTFM/Codi_R/AO_Sex_DE.txt"
# Save the table to a text file
write.table(resultsAO_Sex, file = file_path, sep = "\t", quote = FALSE, row.names = TRUE)

