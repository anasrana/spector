# Installation

Currently you can only install spector from GitHub, there are two methods to do this:

## Installing using Bioconductor

The first step is to check you have the _Bioconductor_ release 3.4 or higher. To do that run:


```r
source("https://bioconductor.org/biocLite.R")
biocLite()
```

This will tell you the version of _Bioconductor_ and ask you to update if necessary. If you are having trouble updating please see the instructions to [troubleshoot](#troubleshooting).

The next step is to install the package with all dependencies and vignettes.

```r
library(devtools)
biocLite("anasrana/spector", 
         dependencies = TRUE,  build_vignettes = TRUE, ref = "v0.8")
```

## Manually install dependencies

The other option is to use `devtools` to install the package with [CRAN](https://cran.r-project.org) dependencies and to install _Bioconductor_ dependencies manually.

#### Bioconductor dependencies

_Versions required in brackets_

- `GenomicAlignments (>= 1.6.3)`
- `GenomicRanges`
- `IRanges (>= 2.7.0)`
- `Rsamtools (>= 1.25.0)`

#### Installing package

```r
library(devtools)

install_github("anasrana/spector", 
               dependencies = TRUE, build_vignettes = TRUE, ref = "v0.8") 
```



# Troubleshooting

_Below is information compiled from the Bioconductor homepage, more information can be found [here](http://www.bioconductor.org/install/#troubleshoot-bioconductor-packages)_

If some of the _Bioconductor_ packages won't update or you are seeing a message like this:

```
BiocInstaller version 3.2 is too old for R version 3.3
```

you have to take the following steps:

- Quit R session
- Start a new session from the command line using `R --vanilla`
- Run this code in R

```R
 remove.packages("BiocInstaller")
 source("https://bioconductor.org/biocLite.R")
 biocValid()
 ```

 If all steps are completed successfully go back to the [installation](#installation)
