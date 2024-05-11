
library(matrixStats)
library(limma)
library(gplots)

dataset_a = readLines("GDS3709.soft")
cleaned_file = dataset_a[!grepl("^[!^#]", dataset_a)]
writeLines(cleaned_file,"GDS3709.soft")

mT = read.table("GDS3709.soft",sep = "\t",header = T,na.strings = "null")



sampls = colnames(mT)[-c(1,2)]

gender = ifelse(grepl("female",sampls),"female","male")


smoking = ifelse(grepl("never smoker",sampls),"never_smoker","smoker")


gender = factor(gender, levels = c("male", "female"))
smoking = factor(smoking, levels = c("smoker", "never_smoker"))


design = model.matrix(~ 0 + gender*smoking)

fit_a = lmFit(mT[-c(1, 2)], design, intercept=TRUE)


contrasts = makeContrasts(....,levels = design)

fit_b = contrasts.fit(fit_a,contrasts_main)

fit_b = eBayes(fit_b)

tb = topTable(fit_b,n=Inf)

