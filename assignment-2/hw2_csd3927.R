#################################
#         BIO476 - R            #
#         Assignment-2          #
#         @csd3927              #
#         Giorgos Geramoutsos   #
#################################


library(matrixStats)
library(limma)
library(gplots)

# PART A

gds_ds = readLines("GDS3709.soft")
gds_cl = gds_ds[!grepl("^[!^#]", gds_ds)]
writeLines(gds_cl,"GDS3709.soft")

myDt = read.table("GDS3709.soft", sep="\t", header=T, na.strings="null")
head(myDt)


gender = factor(rep(c('f', 'm'), each = 40))
gender = gender[2:80]
smoking = factor(rep(rep(c('s', 'ns'), each = 20), 2))
smoking = smoking[2:80]

expr = myDt[,3:ncol(myDt)]
head(expr)

# Q1
design_1 = model.matrix(~ 0 + gender*smoking)
colnames(design_1) = c("female", "male", "smoking", "male_smoking")

fit_1 = lmFit(expr, design_1, intercept=T)
fit_1
rm 

# 1.1 
contrasts_11 = makeContrasts(female-male, levels=design_1)

fit_11 = contrasts.fit(fit_1,contrasts_11)
fit_11 = eBayes(fit_11)
fit_11

tb_11 = topTable(fit_11, n=Inf)
head(tb_11, n=10)


# 1.2 
contrasts_12 = makeContrasts(smoking, levels=design_1)

fit_12 = contrasts.fit(fit_1, contrasts_12)
fit_12

fit_12 = eBayes(fit_12)
fit_12

tb_12 = topTable(fit_12, n=Inf)
head(tb_12, n=10)



# Q2
design_2 = model.matrix(~ 0 + gender+smoking)

colnames(design_2) = c("female", "male", "smoking")

fit_2 = lmFit(expr, design_2, intercept=T)
fit_2


# 2.1
contrasts_21 = makeContrasts(female-male, levels=design_2)
contrasts_21

fit_21 = contrasts.fit(fit_2, contrasts_21)
fit_21

fit_21 = eBayes(fit_21)
fit_21

tb_21 = topTable(fit_21, n=Inf)
head(tb_21, n=10)

# 2.2 
contrasts_22 = makeContrasts("smoking", levels=design_2)
contrasts_22

fit_22 = contrasts.fit(fit_2, contrasts_22)
fit_22

fit_22 = eBayes(fit_22)
fit_22

tb_22 = topTable(fit_22, n=Inf)
head(tb_22, n=10)


# Q3
# diferent approach
library(GEOquery)

gds <- getGEO("GDS3709", GSEMatrix=TRUE)
eset <- GDS2eSet(gds, do.log2=TRUE)


phen = pData(eset)
phen


gender = factor(phen$gender, levels = c("female","male"))
gender
smoking = factor(phen$agent, levels = c("cigarette smoke", "control"))
smoking


# 3.1
smoker_ind = which(phen$agent == "cigarette smoke")
non_smoker_in = which(phen$agent == "control")

ttest_smoking = apply(expr, 1, function(x) t.test(x[smoker_ind], x[non_smoker_in])$p.value)
significant_genes_smoking = which(ttest_smoking < 0.05)
# significant_genes_smoking


for (i in head(significant_genes_smoking, 20)) {
  result <- paste("Gene position:", i, "p-value:", ttest_smoking[i])
  cat(result, "\n")
}



# 3.2
female_ind = which(phen$gender == "female")
male_ind = which(phen$gender == "male")


ttest_gender = apply(expr, 1, function(x) t.test(x[female_ind], x[male_ind])$p.value)

significant_genes_gender = which(ttest_gender < 0.05)

for (i in head(significant_genes_gender, 20)) {
  result <- paste("Gene position:", i, "p-value:", ttest_smoking[i])
  cat(result, "\n")
}



# Q4 


# PART B

# Processing the file 
processFile = function(filepath){
  
  con = file(filepath, "r")
  on.exit(close(con))
  seqs = list()
  
  while(TRUE){
    
    line = readLines(con, n=1)
    
    if(length(line)==0){
      break
    }
    
    isNewSeq = length(grep(">", line, ignore.case=TRUE, perl=TRUE)) > 0
    
    if(isNewSeq){
      motiv = ""
      name = gsub(">([^_]*_[^_ ]*).*", x=line, replacement="\\1", perl=TRUE)
    }else{
      motiv = paste(motiv, line, sep="")
      seqs[[name]] = motiv
    }
  }
  return(seqs) 
}



# 1 
# all patterns list
patterns_a = list()

# exporting the patterns for question 1
getPatternsA = function(ptr_string, base_length=6){
  
  start_pos = 1
  end_pos = nchar(ptr_string) - base_length + 1
  
  v = strsplit(ptr_string, split="")[[1]]
  
  for(i in start_pos:end_pos){
    
    motiv = paste(v[i:(i+base_length-1)], collapse="")
    
    if(!is.null(patterns_a[[motiv]])){
      patterns_a[[motiv]] = patterns_a[[motiv]] + 1
      
    } else{
      patterns_a[[motiv]] = 1
    }
  }
  
  return(patterns_a)
}

human_file_a = processFile("test_file.fa")
# human_file_a = processFile("upstream1000_human.fa")
# chimp_file_a = processFile("upstream1000_chimpanzee.fa")
# mouse_file_a = processFile("upstream1000_mouse.fa")
# rat_file_a = processFile("upstream1000_rat.fa")



cat("Demo File: -test_file.fa-\n\n")
cat("Names and sequences:\n\n")

for(i in seq_along(human_file_a)){
  cat(names(human_file_a)[i], human_file_a[[i]], "\n")
  cat("\n")
}

cat("Total counts for each pattern:\n\n")

for(i in names(human_file_a)){
  print(getPatternsA(human_file_a[[i]]))
}



# 2
#unique patterns list
patterns_b = list()

# exporting all patterns
getPatternsB = function(ptr_string, base_length=6){

  start_pos = 1
  end_pos = nchar(ptr_string) - base_length + 1

  v = strsplit(ptr_string, split="")[[1]]

  # patterns_a = list()
  
  for(i in start_pos:end_pos){

    motiv = paste(v[i:(i+base_length-1)], collapse="")

    if(!is.null(patterns_b[[motiv]])){
      patterns_b[[motiv]] =  1 

    } else{
      patterns_b[[motiv]] = 1
    }
  }

  return(patterns_b)
}


human_file_b = processFile("test_file.fa")
# human_file_b = processFile("upstream1000_human.fa")
# chimp_file_b = processFile("upstream1000_chimpanzee.fa")
# mouse_file_b = processFile("upstream1000_mouse.fa")
# rat_file_b = processFile("upstream1000_rat.fa")

cat("Demo File: -test_file.fa-\n\n")
cat("Names and sequences:\n\n")

for(i in seq_along(human_file_b)){
  cat(names(human_file_b)[i], human_file_b[[i]], "\n")
  cat("\n")
}

cat("Total counts for each pattern:\n\n")

for(i in names(human_file_a)){
  print(getPatternsB(human_file_b[[i]]))
}


















