Class 17:
================
Yi Fu
5/28/2019

## Import data

> **File \> Import \> Network from File…** and choose
> **galFiltered.sif**

> **File \> Import \> Table from File…** and choose **galExpData.csv**

``` r
library(RCy3)
library(igraph)
library(RColorBrewer)
library(ggraph)
```

Let’s check if we can talk to cytoscape

``` r
cytoscapePing()
```

    ## [1] "You are connected to Cytoscape!"

``` r
g = makeSimpleIgraph()
createNetworkFromIgraph(g, "myGraph")
```

    ## Loading data...
    ## Applying default style...
    ## Applying preferred layout...

    ## networkSUID 
    ##        1073

We can simply plot the graph

``` r
plot(g)
```

![](class17_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Or We can include this Cytoscape rendered network image in our report.

``` r
fig <- exportImage(filename="demo", type="png", height=350)
```

    ## Warning: This file already exists. A Cytoscape popup 
    ##                 will be generated to confirm overwrite.

``` r
knitr::include_graphics("./demo.png")
```

![](./demo.png)<!-- -->

Cytoscape provides a number of canned visual
    styles.

``` r
getVisualStyleNames()
```

    ##  [1] "BioPAX_SIF"           "Solid"                "Curved"              
    ##  [4] "size_rank"            "default"              "Ripple"              
    ##  [7] "size_rank_0"          "BioPAX"               "Nested Network Style"
    ## [10] "Sample2"              "Gradient1"            "Sample3"             
    ## [13] "Directed"             "Sample1"              "default black"       
    ## [16] "Big Labels"           "Universe"             "Marquee"             
    ## [19] "Minimal"

``` r
setVisualStyle("Marquee")
```

    ##                 message 
    ## "Visual Style applied."

We can again include this Cytoscape rendered network image in our
report.

``` r
fig <- exportImage(filename="demo_marquee", type="png", height=350)
knitr::include_graphics("./demo_marquee.png")
```

![](./demo_marquee.png)<!-- -->

``` r
prok_vir_cor <- read.delim("./data/virus_prok_cor_abundant.tsv", stringsAsFactors = FALSE)
head(prok_vir_cor)
```

    ##       Var1          Var2    weight
    ## 1  ph_1061 AACY020068177 0.8555342
    ## 2  ph_1258 AACY020207233 0.8055750
    ## 3  ph_3164 AACY020207233 0.8122517
    ## 4  ph_1033 AACY020255495 0.8487498
    ## 5 ph_10996 AACY020255495 0.8734617
    ## 6 ph_11038 AACY020255495 0.8740782

``` r
g <- graph.data.frame(prok_vir_cor, directed = FALSE)
plot(g)
```

![](class17_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Lets turn off text labels

``` r
plot(g, vertex.label=NA)
```

![](class17_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Make the vertex much smaller

``` r
plot(g, vertex.size=3, vertex.label=NA)
```

![](class17_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ggraph(g, layout = 'auto') +
  geom_edge_link(alpha = 0.25) +
  geom_node_point(color="steelblue") +
  theme_graph()
```

    ## Using `nicely` as default layout

![](class17_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
