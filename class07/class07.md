Class 7: More on R functions
================
Yi Fu
4/23/2019

Note: the [original](https://tinyurl.com/rescale-R) *both\_na*,
*gene\_intersect* functions

## 1\. *both\_na* funtion

**Start with a simple example of the larger problem I am trying to
solve**

``` r
x = c(1 ,2, NA, 3, NA)
y = c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
which (c (F, F, T, F, T))
```

    ## [1] 3 5

``` r
which (is.na (c (1, 2, NA, 4)))
```

    ## [1] 3

**Now, I can make this into a function.**

### Version 1:

Check how many NA entries in both input vectors.

``` r
both_na1 <- function (x, y) {
  sum (is.na (x) & is.na (y))
}
```

Examples:

``` r
both_na1 (c(1 ,2, NA, 3, NA), c(NA, 3, NA, 2, NA))
```

    ## [1] 2

``` r
both_na1 (c(NA, NA, NA), c(1, NA, NA))
```

    ## [1] 2

``` r
both_na1 (c(NA, NA, NA), c(1, NA, NA, NA, 1, NA)) # both_na(c(NA, NA, NA,NA, NA, NA), c(1, NA, NA, NA, 1, NA))
```

    ## [1] 4

### Version 2:

Version 1 + check the length of two vectors.

``` r
both_na2 <- function(x, y) {
  if(length(x) != length(y)) {
    print("Input x and y should be vectors of the same length")
  }
  sum(is.na(x) & is.na(y))
}
```

Example:

``` r
both_na2 (c(NA, NA, NA), c(1, NA, NA))
```

    ## [1] 2

``` r
both_na2 (c(NA, NA, NA), c(1, NA, NA, NA, 1, NA))
```

    ## [1] "Input x and y should be vectors of the same length"

    ## [1] 4

### Version 3:

Version 2 + check NA position(s).

``` r
both_na3 <- function(x, y) {
  if(length(x) != length(y)) {
    print("Input x and y should be vectors of the same length")
  }
  na.in.both = is.na(x) & is.na(y)
  na.number = sum(na.in.both)
  na.which = which(na.in.both)

  message("Found ", na.number, " NA's at position(s):", paste(na.which, collapse=", ")) 
  
  return(list(number=na.number, which=na.which))
}
```

Examples:

``` r
both_na3 (c(NA ,2, NA, 3, NA), c(NA, 3, NA, NA, 4))
```

    ## Found 2 NA's at position(s):1, 3

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 1 3

## 2\. *gene\_intersect* function

**Start with a simple example of the larger problem I am trying to
solve**

``` r
df1 = data.frame(IDs=c("gene1", "gene2", "gene3"),
                 exp=c(2,1,1),
                 stringsAsFactors=FALSE)
df2 = data.frame(IDs=c("gene2", "gene4", "gene3", "gene5"),
                 exp=c(-2, NA, 1, 2),
                 stringsAsFactors=FALSE)
x = df1$IDs
y = df2$IDs
```

``` r
intersect(x, y)
```

    ## [1] "gene2" "gene3"

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
y %in% x
```

    ## [1]  TRUE FALSE  TRUE FALSE

``` r
y[y %in% x]
```

    ## [1] "gene2" "gene3"

``` r
cbind( x[ x %in% y ], y[ y %in% x ] )
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

**We can use RStudio shortcut ‘CODE \> EXTRACT FUNCTION’ to turn our
snippet into a working function**

### Version 1:

Find the intersection of two gene lists.

``` r
gene_intersect1 <- function(x, y) {
  cbind( x[ x %in% y ], y[ y %in% x ] )
}
```

Example:

``` r
gene_intersect1 (df1$IDs, df2$IDs)
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

**It seems to contain duplicated information. Let’s improve that.**

### Version 2:

Version 1 + get values for each gene.

``` r
gene_intersect2 <- function(df1, df2) { 
   cbind( df1[ df1$IDs %in% df2$IDs, ], 
          df2[ df2$IDs %in% df1$IDs, "exp"] )
}
```

Example:

``` r
gene_intersect2 (df1, df2)
```

    ##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
    ## 2 gene2   1                               -2
    ## 3 gene3   1                                1

### Version 3:

Version 2 + specify the name of gene column.

``` r
gene_intersect3 <- function(df1, df2, gene.colname="IDs") { 
   cbind( df1[ df1[,gene.colname] %in% df2[,gene.colname], ], 
          exp2=df2[ df2[,gene.colname] %in% df1[,gene.colname], "exp"] )
}
```

Example:

``` r
gene_intersect3 (df1, df2)
```

    ##     IDs exp exp2
    ## 2 gene2   1   -2
    ## 3 gene3   1    1

### Version 4:

Neat version of Version 3.

``` r
gene_intersect4 <- function(df1, df2, gene.colname="IDs") { 

  df1.name <- df1[,gene.colname]
  df2.name <- df2[,gene.colname]

  df1.inds <- df1.name %in% df2.name
  df2.inds <- df2.name %in% df1.name

   cbind( df1[ df1.inds, ], 
          exp2=df2[ df2.inds, "exp"] )
}
```

Example:

``` r
gene_intersect4 (df1, df2)
```

    ##     IDs exp exp2
    ## 2 gene2   1   -2
    ## 3 gene3   1    1

## 3\. *score* function

dropped the lowest score

``` r
score <- function(x) {
  (sum (x) - min (x))/ (length (x) - 1)
}
```

Examples:

``` r
score (c (100, 100, 100, 100, 100, 100, 100, 90))
```

    ## [1] 100

``` r
score (c (100, 90, 90, 90, 90, 90, 97, 80))
```

    ## [1] 92.42857
