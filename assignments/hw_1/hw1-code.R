# assignment 1 code


install.packages("matrixStats")
library(matrixStats)

# Data set 1 

ds_a = readLines("GDS6063.soft")
print(head(ds_a,n=5))

ds_a_cleaned = ds_a[!grepl("[!^#]",ds_a)]
print(head(ds_a_cleaned,n=5))
writeLines(ds_a_cleaned,'GDS6063_cleaned.soft')


ta = read.table("GDS6063_cleaned.soft",sep = "\t",header = T,na.strings = "null")

print(head(ta,n=5))

ta = ta[,3:ncol(ta)]
ta = na.omit(ta)


boxplot(ta)


mo_a = colMeans(as.matrix(ta))
sds_a = colSds(as.matrix(ta))

t_a = t((t(ta)- mo_a)/sds_a)

boxplot(t_a)


heatmap(as.matrix(t_a))



# Data set 2 

ds_b = readLines("GDS6100.soft")
print(head(ds_b,n=5))

ds_b_cleaned = ds_b[!grepl("[!^#]",ds_b)]
print(head(ds_b_cleaned,n=5))
writeLines(ds_b_cleaned,'GDS6100_cleaned.soft')

tb = read.table("GDS6100_cleaned.soft",sep = "\t",header = T,na.strings = "null")

print(head(tb,n=5))

tb = tb[,3:ncol(tb)]
tb = na.omit(tb)


boxplot(tb)


mo_b = colMeans(as.matrix(tb))
sds_b = colSds(as.matrix(tb))

t_b = t((t(tb)- mo_b)/sds_b)

boxplot(t_b)

heatmap(as.matrix(t_b))




