-----------------------------------------------------------

# miRBaseConverter
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/taoshengxu/miRBaseConverter?branch=master&svg=true)](https://ci.appveyor.com/project/taoshengxu/miRBaseConverter)
[![bioc](http://www.bioconductor.org/shields/downloads/miRBaseConverter.svg)](http://bioconductor.org/packages/stats/bioc/miRBaseConverter.html)
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/miRBaseConverter.svg)](http://bioconductor.org/packages/miRBaseConverter/)
[![bioc](http://bioconductor.org/shields/availability/3.6/miRBaseConverter.svg)](http://bioconductor.org/packages/miRBaseConverter/)
[![bioc](http://www.bioconductor.org/shields/build/devel/bioc/miRBaseConverter.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/miRBaseConverter.html)


`miRBaseConverter` is an R/Bioconductor package for converting and retrieving the definition of miRNAs ( Names, Accessions, Sequences, Families and others) in different miRBase versions ( From miRBase version 6 to version 22 [ The latest version ] ). A tiny built-in database is embedded in the `miRBaseConverter` R package for retrieving miRNA information efficiently.

------------------------------------------------------

## Installing miRBaseConverter

```{r,eval=FALSE,warning=FALSE,message=FALSE}
devtools::install_github("taoshengxu/miRBaseConverter")
```

-----------------------------------------------------------------

## Manual
The tutorial and examples of `miRBaseConverter` pacakge can be found [here](https://bioconductor.org/packages/devel/bioc/vignettes/miRBaseConverter/inst/doc/miRBaseConverter-vignette.html).

## Shiny application
The online application with interactive interface of this package can be accessed in
[http://nugget.unisa.edu.au:3838/miRBaseConverter](http://nugget.unisa.edu.au:3838/miRBaseConverter)  or 

[https://taoshengxu.shinyapps.io/miRBaseConverter/](https://taoshengxu.shinyapps.io/miRBaseConverter/).



------------------------------------------------------------------------

## Citation

Please cite the following article when using `miRBaseConverter`:

[![doi](https://img.shields.io/badge/doi.org/10.1186/s12859-018-2531-5-green.svg?style=flat)]( https://doi.org/10.1186/s12859-018-2531-5) [![citation](https://img.shields.io/badge/cited%20by-10-green.svg?style=flat)](https://doi.org/10.1093/bioinformatics/btx378) 

Taosheng Xu, Ning Su, Lin Liu, Junpeng Zhang, Hongqiang Wang, Weijia Zhang, Jie Gui, Kui Yu, Jiuyong Li and Thuc Duy Le. miRBaseConverter: An R/Bioconductor Package for Converting and Retrieving miRNA Name, Accession, Sequence and Family Information in Different Versions of miRBase. BMC Bioinformatics 19 (Suppl 19) :514.
doi: https://doi.org/10.1186/s12859-018-2531-5
