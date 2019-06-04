Class 4: Bioinformatics data analysis with R
================
Yi Fu
4/11/2019

## 1\. Simple Calculations

Example of Addition:

``` r
5 + 3
```

    ## [1] 8

Example of Subtraction:

``` r
5 - 3
```

    ## [1] 2

Example of Multiplication:

``` r
5 * 3
```

    ## [1] 15

Example of Division:

``` r
5 / 3
```

    ## [1] 1.666667

## 2\. Object Assignment

Examples of object Assignment:

``` r
x = 3 * 4
x <- 3 * 4
```

## 3\. Function Calling

Examples of calling *seq()* function:

``` r
seq(1, 10)
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10

``` r
seq(1, 10, by=2)
```

    ## [1] 1 3 5 7 9

Example of calling *date()* function:

``` r
date()
```

    ## [1] "Tue May 28 13:26:57 2019"

## 4\. Getting help

When we don’t remember what a function does and what are the parameters
to take in, we can use *help()* or *?* to check.

``` r
help(seq)
?seq
```

We can also take a look at demos of the function.

``` r
example(seq)
```

    ## 
    ## seq> seq(0, 1, length.out = 11)
    ##  [1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    ## 
    ## seq> seq(stats::rnorm(20)) # effectively 'along'
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    ## 
    ## seq> seq(1, 9, by = 2)     # matches 'end'
    ## [1] 1 3 5 7 9
    ## 
    ## seq> seq(1, 9, by = pi)    # stays below 'end'
    ## [1] 1.000000 4.141593 7.283185
    ## 
    ## seq> seq(1, 6, by = 3)
    ## [1] 1 4
    ## 
    ## seq> seq(1.575, 5.125, by = 0.05)
    ##  [1] 1.575 1.625 1.675 1.725 1.775 1.825 1.875 1.925 1.975 2.025 2.075
    ## [12] 2.125 2.175 2.225 2.275 2.325 2.375 2.425 2.475 2.525 2.575 2.625
    ## [23] 2.675 2.725 2.775 2.825 2.875 2.925 2.975 3.025 3.075 3.125 3.175
    ## [34] 3.225 3.275 3.325 3.375 3.425 3.475 3.525 3.575 3.625 3.675 3.725
    ## [45] 3.775 3.825 3.875 3.925 3.975 4.025 4.075 4.125 4.175 4.225 4.275
    ## [56] 4.325 4.375 4.425 4.475 4.525 4.575 4.625 4.675 4.725 4.775 4.825
    ## [67] 4.875 4.925 4.975 5.025 5.075 5.125
    ## 
    ## seq> seq(17) # same as 1:17, or even better seq_len(17)
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17

## 5\. R Vectors

Let’s play with vectors in R, with the following two vectors *x* and
*y*.

``` r
x = c(56, 95.3, 0.4)
length(x)
```

    ## [1] 3

``` r
y = c(3.2, 1.1, 0.2)
length(y)
```

    ## [1] 3

Here are the vector arithmetics:

``` r
x + y
```

    ## [1] 59.2 96.4  0.6

``` r
x - y
```

    ## [1] 52.8 94.2  0.2

``` r
x * y
```

    ## [1] 179.20 104.83   0.08

``` r
x / y
```

    ## [1] 17.50000 86.63636  2.00000

``` r
sqrt(x)
```

    ## [1] 7.4833148 9.7621719 0.6324555

``` r
round(sqrt(x), 3)
```

    ## [1] 7.483 9.762 0.632

``` r
log(x) / 2 + 1
```

    ## [1] 3.0126758 3.2785149 0.5418546

Get vector values:

``` r
x[1]
```

    ## [1] 56

``` r
x[4]
```

    ## [1] NA

Set vector values:

``` r
x[3] = 0.5
x[4] = 1
```

Now, x is a vector of length 4.

``` r
x
```

    ## [1] 56.0 95.3  0.5  1.0

## RStudio Shortcuts

Use **Alt+Shift+K** to see shortcut quick reference. More shortcuts can
be found
[here](https://support.rstudio.com/hc/en-us/articles/200711853-Keyboard-Shortcuts).

### Some Useful Shortcuts

``` bash
Command+Shift+K    # knit document
Command+Option+I   # Insert chunk (Sweave and Knitr)
Command+Enter      # Run current line/selection
Command+Option+R   # Run the current document
Command+Option+B   # Run from document beginning to current line
Command+Shift+C    # comment and uncomment
```
