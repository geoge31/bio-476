####################3333
# Giorgos Geramoutsos
# csd3927
# midterm bio476



library(limma)
library(gplots)
library(matrixStats)


# EXERCISE 1 

# 1. Φορτώστε το dataset στην R.

data_set = readLines("GDS4894.soft")
data_set_cleaned = data_set[!grepl("[!^#]",data_set)]

# write new file 
writeLines(data_set_cleaned,"GDS4894.soft")

# read in table
newTb = read.table("GDS4894.soft",sep = "\t", header=T, na.strings="null")

head(newTb)

# 2. Πάρτε μόνο τα 2000 πρώτα γονίδια και φτιάξτε το dataset D.

dataset_D = newTb[1:2000, 3:ncol(newTb)]
dataset_D

# 3. Κανονικοποιήστε αφαιρώντας τον μέσο όρο και διαιρώντας με τυπική απόκλιση (στο D).

cat("Pinakas: \n")
dataset_D

cat("\n\nAnastrofos pinakas : \n")
standarized_D = t(dataset_D)
standarized_D

rowmeans_a = rowMeans(standarized_D)
cat("Mesoi oroi grammwn: \n", rowmeans_a, "\n\n")

rowsds_a = rowSds(standarized_D)
cat("Typikes apokleiseis grammwn: \n",rowsds_a, "\n\n")

cat("Telikos Pinakas\n")
final_t = t((standarized_D-rowmeans_a)/rowsds_a)
final_t


# 4. Κάντε ένα boxplot για τα πρώτα 100 γονίδια μονο.

bxpl = final_t[1:100,]
boxplot(bxpl)


# 5. Κάντε ένα t-test μεταξύ trained vs untrained (και τα δύο φύλα μαζί) (από το D)

genderm = factor(rep(c("m"),8))
genderm
genderf = factor(rep(c("f"),6))
genderf

gender = unlist(list(genderm,genderf))
gender

gender = as.factor(gender)
gender



trained = c(1:4, 9:11)
untrained = c(5:8, 12:14)

tb_b  = log(newTb[1:2000, 3:ncol(newTb)])
                  
tbst = t((t(tb_b) - colMeans(tb_b))/colSds(as.matrix(tb_b)))


pvalues = c()


for(i in 1:nrow(tbst))
{
  p_value = t.test(tbst[i, trained], tbst[i, untrained], alternative="greater")$p.value
  pvalues = c(pvalues, p_value)
}

print(pvalues)

# pvalues.adjust(pvalues, method=pvalues.adjust.methods, n=length(pvalues))



# EXERCISE 2 

# 1. Γράψτε μια συνάρτηση που θα επιστρέφει την τέταρτη μεγαλύτερη τιμή από κάποιο vector. Αν
# το μήκος του vector ειναι μικρότερο από 4 θα επιστρέφει NULL.


find_frth_Num = function(vec_a){
  
  if(length(vec_a)<4){
    cat("The length of vector is less than 4\n")
    return(NULL)
  }
  
  big_nums = head(sort(vec_a, decreasing = TRUE), n = 4)
  min_val = min(big_nums)
  all_results = append(all_results, min_val)
  cat("The  4th greatest value of vector is: ", min_val, "\n")

}

my_vec = c(1,80,3,44,120,6,21,11)
my_nvec = c(1,2,3)

cat("calling function with vector : ", my_vec, " ")
find_frth_Num(my_vec)

cat("calling function with vector : ", my_nvec, " ")  
find_frth_Num(my_nvec)




# 2. Γράψτε μια συνάρτηση που θα εφαρμόζει την παραπάνω συνάρτηση (1) σε κάθε γραμμή
# ενός πίνακα Α. Ο πίνακας Α εχει διαστάσεις 20x10 και περιλαμβάνει τιμές από κανονική
# κατανομή με μέση τιμή 100 και τυπική απόκλιση 10.

mtx_a = matrix(rnorm(25,100,10), nrow=5, byrow=T)
mtx_a


for(i in 1:nrow(mtx_a)){
  
  new_vec = as.vector(mtx_a[i,])
  print(new_vec)
  # cat("The 4th greatest number of line is : ")
  find_frth_Num(new_vec)
}


# 3. Γράψτε μια συνάρτηση που θα διάβαζε σωστά ένα αρχέιο με τα εξής περιεχόμενα

dm = read.table("demo.txt", sep='\t', header = T)
dm
