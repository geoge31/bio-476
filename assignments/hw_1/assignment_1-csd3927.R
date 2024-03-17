
### Assignment 1 - BIO476

### Giorgos Geramoutsos - csd3927




# Task 1 

ds1 = readLines('GDS6063.soft')

print(head(ds1,n=10))

ds1_cleaned = ds1[!grepl('[!^#]',ds1)]

print(head(ds1_cleaned,n=10))

writeLines(ds1_cleaned,'GDS6063_cleaned.soft')




ds2 = readLines('GDS6100.soft')
print(head(ds2,n=10))
ds2_cleaned = ds2[!grepl('[!^#]',ds2)]
print(head(ds2_cleaned,n=10))
writeLines(ds2_cleaned,'GDS6100_cleaned.soft')