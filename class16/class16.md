Class 16: Essential UNIX for bioinformatics
================
Yi Fu
5/23/2019

``` r
data=read.csv("mm-second.x.zebrafish.tsv",sep="\t")
colnames(data)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(data)
```

    ##        qseqid         sseqid pident length mismatch gapopen qstart qend
    ## 1 YP_220551.1    NP_059332.1 44.509    346      188       3      1  344
    ## 2 YP_220551.1    NP_059341.1 24.540    163      112       3    112  263
    ## 3 YP_220551.1    NP_059340.1 26.804     97       65       2     98  188
    ## 4 YP_220552.1    NP_059333.1 88.132    514       61       0      1  514
    ## 5 YP_220552.1 XP_021326074.1 31.818     66       32       2    427  482
    ## 6 YP_220552.1 XP_005162943.1 31.818     66       32       2    427  482
    ##   sstart send   evalue bitscore
    ## 1      1  344 8.62e-92    279.0
    ## 2    231  393 5.15e-06     49.7
    ## 3    200  296 1.00e-01     35.8
    ## 4      1  514 0.00e+00    877.0
    ## 5     16   78 6.70e+00     29.3
    ## 6     48  110 7.50e+00     29.6

``` r
hist(data$bitscore,breaks=30)
```

![](class16_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
plot(data$pident  * (data$qend - data$qstart), data$bitscore)
```

![](class16_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
if (!require("ggplot2")) {
  install.packages("ggplot2")
}
```

    ## Loading required package: ggplot2

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

``` r
library(ggplot2)
ggplot(data, aes(pident, bitscore)) + geom_point(alpha=0.1) 
```

![](class16_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggplot(data, aes((data$pident * (data$qend - data$qstart)), bitscore)) + geom_point(alpha=0.1) + geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](class16_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
