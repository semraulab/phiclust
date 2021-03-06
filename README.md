phiclust: A clusterability measure for scRNA-seq data
==================================================
[![DOI](https://zenodo.org/badge/438647623.svg)](https://zenodo.org/badge/latestdoi/438647623)

A package to calculate clusterability in clusters of scRNA-seq data. By
applying this measure you can assess if meaningful variability is left
in clusters. The corresponding paper can be found here [Phiclust: a clusterability measure for single-cell transcriptomics reveals phenotypic subpopulations](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02590-x)

Installation
------------

``` r
# install.packages("devtools")
devtools::install_github("semraulab/phiclust")
```

Example
-------

This how you can use the functions in this R package. The most important
one is phiclust, which calculates the clusterability per cluster.

``` r
library(phiclust)
library(splatter)
library(ggplot2)

#Load sample data simulated with splatter
data("splatO")

expr <- counts(splatO)
expr <- expr[rowSums(expr)>0,]

#Normalize and log-transform the data
expr.norm <- t(t(expr)/colSums(expr))*10000
expr.norm.log <- log(expr.norm + 1)

#Create toy example of a data set
test.cluster <- as.character(splatO$Group)
test.cluster[test.cluster == "Group3"] <- "Group2"
test.cluster[test.cluster == "Group4"] <- "Group2"

#Main funcion that calculates the clusterability
out <- phiclust(expr = expr.norm.log, clusters = test.cluster,
                   exclude = data.frame(clsm = log(colSums(expr) + 1)))
```

For ways to evaluate the results of this clusterability measure check
out the vignette [Guide\_to\_phiclust](https://github.com/semraulab/phiclust/blob/main/vignettes/Guide_to_phiclust.md) and for a real scRNA-seq example you can have a look at our analysis of a fetal kidney [Analysis_kidney](https://github.com/semraulab/phiclust/blob/main/vignettes/Analysis_kidney.md).
