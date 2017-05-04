library(FacileExplorer)
library(FacileTCGADataSet)
library(ggplot2)
theme_set(theme_bw())

tcga <- FacileTCGADataSet()

samples <- tcga %>%
  filter_samples(indication == "LUSC")
head(samples)
##    dataset                    sample_id
## 1|    LUSC TCGA-18-3406-01A-01R-0980-07
## 2|    LUSC TCGA-18-3407-01A-01R-0980-07
## 3|    LUSC TCGA-18-3408-01A-01R-0980-07

features <- tcga %>%
  filter_features(name %in% c("CCNB2", "KRT17")) %>%
  mutate(assay="rnaseq")
select(features, feature_type, feature_id, name, assay)
# feature_type feature_id  name  assay
#        <chr>      <chr> <chr>  <chr>
#       entrez       9133 CCNB2 rnaseq
#       entrez       3872 KRT17 rnaseq

## I am independently exploring CCNB2 vs KRT17 in tumors and normals
dat <- samples %>%
  with_sample_covariates('sample_type') %>%
  with_assay_data(features, normalized=TRUE, spread="name")
head(dat)
# dataset                    sample_id    CCNB2     KRT17 sample_type
#   <chr>                        <chr>    <dbl>     <dbl>      <fctr>
#    LUSC TCGA-18-3406-01A-01R-0980-07 5.091566 13.758303       tumor
#    LUSC TCGA-18-3407-01A-01R-0980-07 4.774941 12.631085       tumor
#    LUSC TCGA-18-3408-01A-01R-0980-07 5.239986  9.227850       tumor
#    LUSC TCGA-18-3409-01A-01R-0980-07 4.645760 12.166467       tumor

ggplot(dat, aes(CCNB2, KRT17, color=sample_type)) +
  geom_point(size=2.5) +
  theme(legend.position = "bottom")

## Collaborator informs me to retrieve custom covariate that
## was just added to the system.
## The with_* functions decorate results with more data from the
## FacileDataSet
dat <- with_sample_covariates(dat, 'plotcon2017')
head(dat)
# dataset                    sample_id    CCNB2     KRT17 sample_type plotcon2017
#   <chr>                        <chr>    <dbl>     <dbl>      <fctr>       <chr>
#    LUSC TCGA-18-3406-01A-01R-0980-07 5.091566 13.758303       tumor  KRT17_high
#    LUSC TCGA-18-3407-01A-01R-0980-07 4.774941 12.631085       tumor  KRT17_high
#    LUSC TCGA-18-3408-01A-01R-0980-07 5.239986  9.227850       tumor  KRT17_high
#    LUSC TCGA-18-3409-01A-01R-0980-07 4.645760 12.166467       tumor   unlabeled
#    LUSC TCGA-18-3410-01A-01R-0980-07 7.177209  7.319666       tumor   KRT17_low


ggplot(dat, aes(CCNB2, KRT17, color=plotcon2017)) +
  geom_point(size=2.5) +
  theme(legend.position = "bottom")

## Let me REALLY check this out and add my own thoughts
features <- collect(features)
live <- fscatter(tcga, features[1,], features[2,], samples=samples)

## We can access custom annotations I brushed onto these data within the
## `live$annotation` data.frame
head(live$annotation)
# dataset                    sample_id     variable     value            type
#   <chr>                        <chr>        <chr>    <fctr>           <chr>
#    LUSC TCGA-22-1017-01A-01R-0692-07 active_brush low_mygrp user_annotation
#    LUSC TCGA-22-4594-01A-01R-1201-07 active_brush low_mygrp user_annotation
#    LUSC TCGA-22-5481-01A-31R-1949-07 active_brush low_mygrp user_annotation
#    LUSC TCGA-37-4130-01A-01R-1100-07 active_brush low_mygrp user_annotation

## Save it so that collaborator can look
save_custom_sample_covariates(tcga, 'check_it', live$annotation)


## Don't like dplyr-ish interface? Can convert to "standard" bioconductor
## data structures
library(edgeR)
y <- as.DGEList(samples)

