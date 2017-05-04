library(FacileExplorer)
library(FacileTCGADataSet)
library(ggplot2)
theme_set(theme_bw())

tcga <- FacileTCGADataSet()

samples <- tcga %>%
  filter_samples(indication == "LUSC")

features <- tcga %>%
  filter_features(name %in% c("CCNB2", "KRT17")) %>%
  mutate(assay="rnaseq") %>%
  collect

## I am independently exploring CCNB2 vs KRT17 in tumors and normals
dat <- samples %>%
  with_sample_covariates('sample_type') %>%
  with_assay_data(features, normalized=TRUE, spread="name")

ggplot(dat, aes(CCNB2, KRT17, color=sample_type)) +
  geom_point(size=2.5) +
  theme(legend.position = "bottom")

## Collaborator informs me to retrieve custom covariate that
## was just added to the system.
## The with_* functions decorate results with more data from the
## FacileDataSet
dat <- with_sample_covariates(dat, 'plotcon2017')
View(dat)

ggplot(dat, aes(CCNB2, KRT17, color=plotcon2017)) +
  geom_point(size=2.5) +
  theme(legend.position="bottom")

## Invoke the FacileExplorer scatterplot module as a shiny gadget so that we
## can explore this interactively within RStudio, as well as capture custom
## brushings we apply.
live <- fscatter(tcga, features[1,], features[2,], samples=samples)

## We can access custom annotations I brushed onto these data within the
## `live$annotation` data.frame
View(live$annotation)

## Save it so that collaborator can look
save_custom_sample_covariates(tcga, 'check_it', live$annotation)


## Don't like dplyr-ish interface? Can convert to "standard" bioconductor
## data structures. This creates a DGEList with the expression data and all
## sample covariate information we have for the `samples` specified.
library(edgeR)
y <- as.DGEList(samples)

