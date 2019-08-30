library(here)
if (!requireNamespace("multiGSEA")) {
  devtools::install_github("lianos/multiGSA@develop")
}
library(multiGSEA)
library(FacileData)
library(FacileAnalysis)
library(dplyr)
library(plotly)
theme_set(theme_bw())

# Setup some cpu-local parameters ----------------------------------------------
fdata.dir <- "~/workspace/data/FacileData/dockerlink"
dat.dir <- here("2019-08-29_ssf-compbio", "data")


# Code with Me -----------------------------------------------------------------
tcga <- FacileDataSet(file.path(fdata.dir, "FacileTcgaDataSet"))

# What data do we have in here?
assay_info(tcga)

(xs <- samples(tcga))

# The facile_frame `xs` looks like a "normal" tibble, but has an internal
# connection to the remote data store that it was pulled out of.
#
# We can use the remote connection to send the facile_frame back to the
# remote data store and decorate with more data found there.

# Perhaps we forgot to pull out some sample-level covariates:
(xs <- with_sample_covariates(xs, "indication"))

# There are facile_frame S3 methods defined for dplyr's munging verbs
dat <- xs %>%
  filter(indication %in% c("BLCA", "BRCA", "LUSC", "UCEC")) %>%
  with_sample_covariates(c("sex", "sample_type"))
head(dat)

# HER2 CNV vs Expression -------------------------------------------------------
# Let's explore the relationship between HER2 amplification and its expression
# across the entire data set over the subset of data we are working with.
(her2 <- features(tcga) %>% filter(name == "ERBB2"))

# Pull out gistic2 cnv scores for samples and plot

adat <- with_assay_data(dat, her2, assay_name = "cnv_score") # ~1s for 2067 samples

ggplot(adat, aes(ERBB2, color = indication)) +
  geom_density() +
  xlab("GISTIC2 score") +
  ggtitle("GISTIC2 CNV scores per sample")

# Add expression data to `adat`
adat <- with_assay_data(adat, her2, normalized = TRUE,
                        assay_name = "rsem_gene") # ~ 1s
head(adat)

ggplot(adat, aes(x = ERBB2, y = rsem_gene_ERBB2, color = sex)) +
  geom_point() +
  facet_wrap(~ indication)

# FacileAnalyses ---------------------------------------------------------------
library(edgeR)
(crc <- filter_samples(tcga, indication == "CRC"))

# Bioc-way .....................................................................
y.crc <- as.DGEList(crc) # 60498 genes, 430 samples
head(y.crc$samples)

des <- model.matrix(~ sample_type + sex, y.crc$samples)
y <- y.crc[filterByExpr(y.crc, des),,keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

vm <- voom(y, des, plot = TRUE)
fit <- lmFit(vm, vm$design) %>% eBayes()
(tt.vm <- topTable(fit, coef = "sample_typetumor", n = Inf))


# facile way ...................................................................
# FacileAnalysis::fdge() -- Facile Differential Gene Expression

?fdge

dge.tumor <- crc %>%
  flm_def(covariate = "sample_type",
          numer = "tumor", denom = "normal",
          batch = "sex") %>%
  fdge(method = "voom", filter = tt.vm$feature_id)

# Lets interact with our results)
(tt.facile <- tidy(dge.tumor))


cmp <- inner_join(
  transmute(tt.vm, feature_id, symbol, logFC, t, pval = P.Value),
  select(tt.facile, feature_id, logFC, t, pval),
  by = "feature_id")
plot(cmp$t.x, cmp$t.y, xlab = "bioc", ylab = "facile", main = "t-statistics")
abline(0, 1, col = "red")

# Deep interaction with shine
shine(dge.tumor)

# JS interactivity for reporting with viz() [used by shine]
viz(dge.tumor)

# Where does new data come from? ...............................................
# Guinney 2015 published new consensus molecular subtyping for CRC
crc.info <- read.csv(file.path(dat.dir, "guinney-2015_crc-cms.csv"),
                     stringsAsFactors = FALSE)
head(crc.info)

# Our FacileTtestAnalysisResult has brings with it all the information it needs
# to reconstruct its results.
#
# sample_type and sex are there because they were used to define the linear
# model
head(samples(dge.tumor))

# We can pull these samples out, merge them with new data and analyze subtypes
xsamples <- samples(dge.tumor) %>%
  filter(sample_type == "tumor") %>%
  transmute(dataset, sample_id, xref = substr(sample_id, 1, 12)) %>%
  inner_join(crc.info, by = c("xref" = "sample")) %>%
  select(-xref)

# Pretend you didn't know how to code, but you still wanted to analyze these
# newly annotated data
cms.anova <- fdgeGadget(xsamples)

cms.2v1 <- xsamples %>%
  flm_def("cms_label", "CMS2", "CMS1") %>%
  fdge(method = "voom", filter = features(cms.anova))
cms.4v2 <- xsamples %>%
  flm_def("cms_label", "CMS4", "CMS2") %>%
  fdge(method = "voom", filter = features(cms.anova))
cmp2 <- compare(cms.4v2, cms.2v1)

# gsea -------------------------------------------------------------------------
library(multiGSEA)
gdb <- getMSigGeneSetDb("h", "human", "ensembl")
gsea <- ffsea(cms.2v1, gdb, methods = c("cameraPR", "enrichtest"))
shine(gsea)

gsea2 <- ffseaGadget(cms.4v2, gdb)

crc.pca <- fpca(xsamples)
crc.pca.gsea <- ffsea(crc.pca, gdb, dim = 1)
