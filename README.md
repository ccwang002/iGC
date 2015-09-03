## iGC
iGC is an integrated analysis package of gene expression and copy number alteration. It can identify differentially expressed gene driven by copy number alterations (CNA) from samples with both gene expression and CNA data.


### Installation
One can directly install iGC from Bioconductor 3.2+,

```{r}
source("http://www.bioconductor.org/biocLite.R")
biocLite('iGC')
```

Otherwise, including the development version, one can install iGC by [devtools].

```{r}
install.packages("devtools")
devtools::install_github("ccwang002/iGC", build_vignettes = TRUE)
```

> Note: Some requisites are needed on some platforms. 
> 
> For example, [Rtools] is required on Windows; Xcode is required on OSX; A compiler and various development libraries are needed on Linux (details vary across different distros).


### Documentation
One can view the vignettes and man pages by the following commands,

```{r}
browseVignettes("iGC")
help(iGC)
```

### Citation
If you use iGC package, please consider citing our publication.

> Yi-Pin Lai, Liang-Bo Wang, Liang-Chuan Lai, Mong-Hsun Tsai, Tzu-Pin Lu, and Eric Y. Chuang. iGC–an integrated analysis package of Gene expression and Copy number alteration, *(pending publication)*, (2015).


[devtools]: https://github.com/hadley/devtools
[Rtools]: https://cran.r-project.org/bin/windows/Rtools/
