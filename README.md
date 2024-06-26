# BIO-476
Genetics


### Libraries

-   limma

-   matrixStats

-   gplots

-   GEOquery

```{r}
library(limma)
library(gplots)
library(matrixStats)
library(GEOquery)
```

------------------------------------------------------------------------

### Functions

```{r}
doubleNum = function(num){
  dn = 2*num
  return(dn)
}

mn = doubleNum(15)
mn
```

------------------------------------------------------------------------

# Κατανομες

-   **Ομοιομορφη Κατανομη**

    -   v = runif(Χ, Υ, Ζ)

        -   δημιουργια βεκτορ Χ τιμων απο Υ εως Ζ

        -   v = runif(100,10,20)

-   **Κανονικη Κατανομη**

    -   r = rnorm(N, Ε, S)

        -   δημιουργια βεκτορ με Ν τιμες με μεση τιμη E και τυπικη αποκλιση S

        -   r = rnorm(200, 10, 6)

-   **Εκθετικη κατανομη**

    -   r = rexp(N, rate=n)

        -   δημιουργια βεκτορ Ν τιμων και ρυθμο n

        -   r = rexp(1000, rate=10)

------------------------------------------------------------------------

# 1. Vectors

```{r}
w = 10 
w2 = c(1,2,3,4,5,6,1,7,8,2,1)
w2

```

## 1.1 Vector arithmetic operations

```{r}
v = c(1,2,3)
u = c(4,5,6)

a = v + u 
b = v * u
c = v %*% u
d = log(v)

cat("Operations: \n","v + u =(",a,")\n","v * u =(",b,")\n","v %*% u =",c,"\n","log(v) =",d,"\n")


```

```{r}
# vector με ολους τους ζυγους απο 2 εως 100
l = 2*(1:50)

# vector με ολους τους ζυγους απο 20 εως 100
l2 = 2*(10:50)

cat("l = (",l,")\n\n","l2 = (",l2, ")\n")


```

------------------------------------------------------------------------

## 1.2 Accessing Vector

```{r}
v = runif(100, 10, 50)
v

# 1st element
v[1]

# change 1st elements
v[1] = 21
v[1]


# δημιουργια βεκτορ απο βεκτορ
w = v[1:10]
w

# δημιουργια βεκτορ απο βεκτορ παιρνοντας μονο τις ζυγες θεσεις
w = v[2*(1:50)]
w


# Ολα τα στοιχεια εκτος απο το 10 
w = v[1:10]
w[-10]

# ολα τα στοιχεια εκτος απο το 3 το 5 και το 9 
w[c(-3,-5,-9)]

```

```{r}
a = c(1,2,3)
b = c("str1","b",3)
a
b
```

\~ Parsing στοιχεια βεκτορ μεσω TRUE, FALSE

```{r}
# d = runif(100, 7, 10)
d = c(1,3,5)
d
d[c(T,T,F)]

```

-   **Κανονικη κατανομη (rnorm)**

    -   `u = rnorm(100, 10, 1)`

        -   κανονική κατανομή με 100 τιμων, μεση τιμη 10, τυπικη αποκλιση 1

```{r}
z = rnorm(10,5,1)
z
```

------------------------------------------------------------------------

## 1.3 Vector Functions

-   `mean(vec)`

    -   μεση τιμη τιμων του βεκτορ

```{r}
mean(z)
```

-   `which(vec)`

    -   επιστρεφει τις τιμες του vector που ειναι TRUE

```{r}
u2 = c(T,T,F,T,F,F,F,T,T,T)
which(u2)
```

-   `which.max(vec)`

    -   επιστρεφει την θεση με μεγαλυτερη τιμη

-   `which.min(vec)`

    -   επιστρεφει την θεση με την μικροτερη τιμη

```{r}
cat("indexes of max and min\n","max =",which.max(z),"\n","min=", which.min(z))

```

-   `max(vec), min(vec)`

    -   επιστρεφει τις max & min τιμες ενος βεκτορ

```{r}
cat("max value =", max(z), "\n","min value =", min(z))
```

-   `length(vec)`

    -   μηκος του βεκτορ

```{r}
cat("length of vector : ",length(z))
```

-   `sum(vec)`

    -   αθροισμα των στοιχειων του βεκτορ

```{r}
cat("sum of vector's values :",sum(z))
```

-   `median(vec)`

    -   παιρνω την μεσαια τιμη του βεκτορ

```{r}
q1 = c(1,3,5,6,8,9) # zygos arithmos --> mesh timh twn mesaiwn
cat("median q1 = ",median(q1),"\n")
 
q2 = c(1,3,5,6,8) # monos arithmos --> mesaia timh
cat("median q2 = ",median(q1))
```

-   `sd(vec)`

    -   τυπικη αποκλιση του βεκτορ

```{r}
cat("sd = ",sd(z))
```

-   `var(vec)`

    -   διασπορα

```{r}
cat("var = ", var(z))
```

-   `sort(vec)`

    -   αυξουσα σειρα τιμων

```{r}
cat("vector = (", z, ")\n\n", "sort(vector) = ", sort(z), "\n")
```

-   `order(vec)`

    -   επιστρεφει τις θεσεις της αυξουσας σειρας των τιμων του βεκτορ

```{r}
cat("Indexes of sort :",order(z))
```

------------------------------------------------------------------------

## 1.4 Απεικονιση Vectors

-   **Plots**

    -   `plot(vector)`

Δυο διαφορετικα clusters & χρωματισμος

```{r}
greeks = rnorm(500,10,6)
sweedish = rnorm(500,20,6)

mx = c(greeks,sweedish)

mycol = rep(c("navyblue","yellow"),each=500) # ftiaxnoume to xrwma 

plot(mx, xlim=c(0,1000), ylim=c(0,40), col=mycol, xlab="1000 samples", ylab="values")
```

-   **Boxplots**

    -   `boxplot(vec)`

```{r}
boxplot(list("greeks"=greeks,"sweedish"=sweedish), col=c("blue","yellow"))
```

```{r}
myExp = rexp(500, rate=10)
boxplot(myExp)
```

-   **Density**

    -   `plot(density(vec))`

```{r}

r1 = runif(100,10,20)
r2 = rnorm(100,10,6)
r3 = rexp(100, rate=10)
mn3 = mean(r3)
# mn3

plot(density(r1), main = "R - uniform")
plot(density(r2), main = "R - normal")
plot(density(r3), main = "R - exponential", xlab="Time")


# abline(v=mn3, col="blue")
```

-   **Histogram**

    -   `hist(vec)`

```{r}

hist(r1, main="R uniform", col="lightblue")
hist(r2, main="R normal", col="orange")
hist(r3, main="R exponential", col="lightgreen")

# dev.off()
```

Create pdf

```{r}
pdf("r2.pdf")
hist(r2, main="R uniform", col="lightblue", xlab="Values of normal distribution with mean around 10")
dev.off()

```

------------------------------------------------------------------------

## Examples in vectors

**\~** Πως θα παρω ολες τις τιμες του βεκτορ που ειναι **\> 9**

```{r}

u = rnorm(100,10,1)

u>9 # θα επιστρεψει true false στις θεσεις που οι τιμες του βεκτορ ειναι μεγαλυτερες απο 9

a[a>9] # θα επιστρεψει τις τιμες του βεκτορ που ειναι μεγαλυτερες απο 9
```

------------------------------------------------------------------------

\~ Θελω να παρω ολα τα ατομα που εχουν ηλικια **\> 25**

```{r}

ilikies = c(19,27,31,21)
onomata = c("giorgos", "maria", "kostas", "foteini")

which(ilikies>25) # επιστρεφει τις θεσεις του βεκτορ που ειναι T 

onomata[ilikies>25] # επιστρεφει τα ονοματα που εχουν τιμη μεγαλυτερη απο 25 


```

------------------------------------------------------------------------

\~ Παιρνουμε 100 ατομα με μεση τιμη βαρους 60

```{r}
weight = rnorm(100, 60, 5)
height = 2 * weight + rnorm(100, 0, 10) # συσχετιση υψους και βαρους
# οσο αυξανεται το υψος μεγαλωνει και το βαρος 

plot(height, weight, xlab="Height", ylab="Weight", main="Correlation between heigth and weight")
```

------------------------------------------------------------------------

# 2. Read & Clean Dataset

```{r}
myDataSet = readLines("file.soft")
myDCleaned = myDataSet[!grepl("[!^#]",myDataSet)]

# write new file 
writeLines(myDCleaned,"newFile.soft")
```

------------------------------------------------------------------------

# 3. Matrices

## 3.1 Δημιουργια Πινακα

-   `byrow = T`

    -   προσθετει τις τιμες του πινακα κατα γραμμες και οχι κατα στηλες οπως κανει by default η R

-   `dim(matrix)`

    -   επιστρεφει τις διαστασεις του πινακα (rows X cols)

```{r}
# δημιουργια πινακα 
mtrx = matrix(1:100, nrow=10, ncol=10)
cat("Matrix:\n")
mtrx


# δημιουργια πινκα με τυχαιες τιμες κανονικης κατανομης
mtrx2 = matrix(rnorm(100,10,10), nrow=10)
mtrx2

dim(mtrx2)
```

\~ Δημιουργια πινακα με τιμες απο μια κανονικη κατανομη με μεση τιμη 100 και τυπικη αποκλειση 10. Ο πινακας πρεπει να ειναι 20χ10

```{r}
mtrx3 = matrix(data=rnorm(200,100,10), nrow=20)
mtrx3
```

\~ Δημιουργια πινακα απο μια ομοιομορφη κατανομη με ελαχιστο 10 και μεγιστο 100

```{r}
mtrx4 = matrix(data=runif(200,10,100), nrow=20)
mtrx4

cat("\n\nThe min value of matrix is: ", min(mtrx4),"\n")
```

## 3.2 Προσβαση στα στοιχεια του πινακα

```{r}
mtrx5 = matrix(1:25, nrow=5, ncol=5)
mtrx5

cat("\nMatrix[1,4] = ", mtrx5[1,4], "\n")
cat("\nMatrix[1:2,4] = ",mtrx5[1:4, 4], "\n")
cat("\nMatrix[5,3:5] = ",mtrx5[5, 3:5], "\n")

```

-   H R διαβαζει τους πινακες σαν vectors κατα στηλες

```{r}
mtrx5[mtrx5>13]
```

```{r}
which(mtrx5>20)
```

## 3.3 Πραξεις μεταξυ πινακων και βεκτορ

```{r}

newMatrix = matrix(1:100, nrow=20, byrow=T) 
newMatrix

newMatrix + 1
```

```{r}
# m1 = matrix(runif(100,0,10), nrow=20)
# m2 = matrix(rnorm(100,10,10), nrow=20)

m1 = matrix(1:100, nrow=10)
m2 = matrix(101:200, nrow=10)

m3 = m1 + m2 
m4 = m1 * m2
m5 = m1 ^ m2
cat("matrix 1 :\n\n",m1,"\n\nmatrix 2 :\n\n",m2,"\n\nmatrix 3 (matrix 1 + matrix 2) :\n\n",m3,"\n\n matrix 4 (matrix 1 * matrix 2) : \n\n",m4, "\n\nmatrix 5 (matrix 1 ^ matrix 2) :\n\n", m5)
# m2
# m3
```

## 3.4 Συνενωση πινακων

-   ειτε κατα γραμμες ειτε κατα στηλες

-   `cbind(matrix_1,matrix_2)`

```{r}

tab1 = matrix(1:10, nrow=5, byrow=T)
tab2 = matrix(1:10, nrow=5, byrow=T)

tb = cbind(tab1,tab2)
tb
```

-   `colSums(matrix)`

    -   επιστρεφει το αθροισμα των τιμων της καθε στηλης του πινακα

-   `rowSums(matrix)`

    -   επιστρεφει το αθροισμα των τιμων της καθε γραμμης του πινακα

```{r}
cat("Sum of each column of matrix tb is: ", colSums(tb), "\n")
cat("Sum of each of row matrix tb is: ", rowSums(tb), "\n")
```

-   `colMeans(matrix)`

    -   επιστρεφει τους μεσους ορους της καθε στηλης του πινακα

-   `rowMeans(matrix)`

    -   επιστρεφει τους μεσους ορους της καθε γραμμης του πινακα

```{r}
cat("Mean of each column of matrix tb is: ", colMeans(tb), "\n")
cat("Mean of each row of matrix tb is: ", rowMeans(tb), "\n")
```

-   `colSds(matrix)`

    -   επιστρεφει την τυπικη αποκλιση καθε στηλης του πινακα

-   `rowSds(matrix)`

    -   επιστρεφει την τυπικη αποκλιση καθε γραμμης του πινακα

```{r}
tab3 = matrix(rnorm(200,10,5), nrow=20)
tab3
colSds(tab3)
rowSds(tab3)
```

## 3.5 Kανονικοποιηση Πινακα

1.  Βρισκω τον μεσο ορο καθε στηλης

2.  Βρισκω την τυπικη αποκλιση καθε στηλης

3.  Για καθε στοιχειο της στηλης αφαιρω τον μεσο ορο και διαιρω με την τυπικη αποκλιση

    3.1 Για να αφαιρεσω απο καθε στοιχειο τον μεσο ορο χρειαζεται να αναστρεψω τον πινακα

```{r}

mx_a = matrix(rnorm(20,10,3), nrow=4)
cat("pinakas: \n")
mx_a
cat("\n")

# μεσος ορος καθε στηλης
colmeans_a =  colMeans(mx_a)
cat("mesoi oroi sthlhs : ", colmeans_a, "\n\n")

# τυπικη αποκλιση καθε στηλης
colsds_a = colSds(mx_a)
cat("typikes apokleiseis: ", colsds_a, "\n\n")

# για καθε στοιχειο αφαιρω με τον μ.ο. και διαιρω με την τ.α.
  # βρισκοντας τον αναστροφο πινακα 
cat("Anastrofos tou pinaka : \n\n")
mx_b = t(mx_a)
mx_b
cat("\n")

# μεσος ορος καθε γραμμης 
rowmeans_a = rowMeans(mx_b)
cat("mesoi oroi grammhs : ", rowmeans_a, "\n\n")

# τυπικες αποκλισεις γραμμης 
rowsds_a = rowSds(mx_b)
cat("typikes apokliseis grammhs : ", rowsds_a, "\n\n")

mx_b = mx_b - rowmeans_a
cat("exontas afairesei tous m.o apo tis grammes\n")
mx_b 
cat("\n\n")


mx_b = mx_b / rowsds_a
cat("exontas diairesei me tis typikes apokleiseis\n")
mx_b
cat("\n\n")


cat("epanafora toy pinaka sthn arxikh tou morfi\n")
mx_b = t(mx_b)
mx_b
```
