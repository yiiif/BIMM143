Class 12: Bioinformatics in drug discovery and design
================
Yi Fu
5/9/2019

First, let’s check if “bio3d” package is installed. And then, load the
package.

``` r
library(bio3d)
```

``` r
file.name <- get.pdb("1hsg.pdb")
```

    ## Warning in get.pdb("1hsg.pdb"): ids should be standard 4 character PDB-IDs:
    ## trying first 4 characters...

    ## Warning in get.pdb("1hsg.pdb"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
```

``` r
prot <- atom.select(hiv, "protein",value=T)
lig <- atom.select(hiv, "ligand",value=T)
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

``` bash
#~/Downloads/autodock_vina_1_1_2_mac/bin/vina --config config.txt --log log.txt
```

## process our docking results

``` r
res <- read.pdb("all.pdbqt", multi=T)
write.pdb(res, "results.pdb")
```

## calculate the RMSD (root mean square distance) between each of the docking results and the known crystal structure using the bio3d package.

``` r
# res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

## Normal Mode analysis for flexibility prediction

``` r
pdb <- read.pdb("1HEL")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.017 seconds.
    ##  Diagonalizing Hessian...    Done in 0.092 seconds.

``` r
plot(modes, sse=pdb)
```

![](class12_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
mktrj(modes, mode=7, file="nma_7.pdb")
```
