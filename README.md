# Overview

This repository holds presentation materials used for seminars about FacileData.

The FacileData ecosystem is a set of R packages built to make the analysis and exploration of large, consortia-scale genomic data more fruitful for computational and non-computational biologists, alike. It also provides a mechanisms that enables better collaboration between these two parties by enabling them to pass analyses off to each other "in flight" so that they can work together on a common goal.

We plan to open source these tools in September, 2017.

## FacileDataSet

The `FacileDataSet` package provides a novel container for the out-of-memory storage of multi-assay, multi-dataset, genomics data (think of the type of data generated for the [TCGA][tcga]).

It provides fast and efficient ways to query and retrieve these data, as well as a [dplyr][dplyr]-like syntax to more easily manipulate them within the context of an exploratory data analysis. These data can also be easily coerced into more standard [bioconductor][bioconductor] assay containers, like an `ExpressionSet`, `DGEList`, etc. to better integrate with the rest of the bioconductor ecosystem.

## FacileExplorer

The `FacileExplorer` package provides a set of independent [shiny modules][modules], each of which interacts with a `FacileDataSet` in order to produce a specific mode of interaction and/or analysis of arbitrary subsets of the data. For instance you can invoke an interactive scatter plot, PCA, boxplot, differential expression analysis, etc.

These modules can be used in two modes:

1. Modules can be invoked as [shiny gadgets][gadget] by an analyst while she is in the middle of an active (code-driven) analysis; or
2. Modules can be combined together into a stand-alone shiny application (also included in this package) that gives non-informaticians the ability to independently explore and analyze these data using a GUI.

# Talks

## 2017-plotcon

This was the first public debut of FacileData, given on May 3, 2017 at [PLOTCON][plotcon]. We outlined the need for tools that enable better collaboration between biologists and computational biologists, and showcased how the `FacileDataSet` and `FacileExplorer` serves this purpose.

The [video of this talk][plotconvideo] has been [posted to YouTube][plotconvideo].

[tcga]: https://cancergenome.nih.gov/
[gadget]: https://shiny.rstudio.com/articles/gadgets.html
[modules]: https://shiny.rstudio.com/articles/modules.html
[bioconductor]: https://bioconductor.org
[plotcon]: https://plotcon.plot.ly/
[plotconvideo]: https://www.youtube.com/watch?v=-qepBa5vYxU
