R Notebook
================

  - [INSTALLATION DES “PACKAGE” POUR L’ANALYSE DES
    DONNÉES](#installation-des-package-pour-lanalyse-des-données)
  - [Installation Dada2](#installation-dada2)
  - [Installation phyloseq](#installation-phyloseq)
  - [Installation Sudo](#installation-sudo)
  - [Installation gridExtra](#installation-gridextra)
  - [Installation Cran\_packages](#installation-cran_packages)
  - [Installation Deseq2](#installation-deseq2)
  - [Installation Rmarkdown](#installation-rmarkdown)
  - [Installation knitr](#installation-knitr)

# INSTALLATION DES “PACKAGE” POUR L’ANALYSE DES DONNÉES

# Installation Dada2

``` r
library(ggplot2)
```

``` r
BiocManager::install("dada2")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'dada2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'cli', 'GenomeInfoDb', 'pillar'

# Installation phyloseq

``` r
library(phyloseq)
```

``` r
BiocManager::install("phyloseq")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phyloseq'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'cli', 'GenomeInfoDb', 'pillar'

# Installation Sudo

``` bash
sudo apt-get install -y libglpk-dev
```

    ## sudo: unable to resolve host 9143e498ed98: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## libglpk-dev is already the newest version (4.65-2).
    ## 0 upgraded, 0 newly installed, 0 to remove and 28 not upgraded.

# Installation gridExtra

``` r
install.packages("gridExtra")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

# Installation Cran\_packages

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
.inst <- .github_packages %in% installed.packages()
if (any(!.inst))
  devtools::install_github(.github_packages[!.inst])
```

    ## Skipping install of 'phyloseqGraphTest' from a github remote, the SHA1 (3fb6c274) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
BiocManager::install(".bioc_packages")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) '.bioc_packages'

    ## Warning: package '.bioc_packages' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'cli', 'GenomeInfoDb', 'pillar'

``` r
#if(any(!.inst)){source("http://bioconductor.org/biocLite.R")
#biocLite(.bioc_packages[!.inst])}
```

# Installation Deseq2

``` r
BiocManager::install("DESeq2")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'DESeq2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'cli', 'GenomeInfoDb', 'pillar'

# Installation Rmarkdown

``` r
install.packages("rmarkdown")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

# Installation knitr

``` r
install.packages("knitr")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)
