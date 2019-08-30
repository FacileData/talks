# Parses the `clinical_molecular_public_all.txt` metadata from Guinney et al.
# to retrieve new subtyping scheme and other covariates.
#
# Paper: https://www.ncbi.nlm.nih.gov/pubmed/26457759
# Data File: https://www.synapse.org/#!Synapse:syn4978510
library(here)

orig.fn <- file.path(
  "~/workspace/data/publications/guinney_natmed_2015_crc_cms",
  "synapse/clinical_molecular_public_all.txt")

out.dir <- here("2019-08-29_ssf-compbio", "data")
stopifnot(dir.exists(out.dir))

library(dplyr)
library(readr)
dat.all <- read_tsv(orig.fn)
dat <- dat.all %>%
  filter(grepl("^TCGA", sample, ignore.case = TRUE)) %>%
  transmute(sample,
            cms_label,
            stage = paste("stage", stage, sep = "_"),
            msi = msi,
            kras_mut = paste0("kras_", ifelse(kras_mut == 0, "WT", "mutant")),
            bram_mut = paste0("braf_", ifelse(braf_mut == 0, "WT", "mutant")))
write.csv(dat, file.path(out.dir, "guinney-2015_crc-cms.csv"), row.names = FALSE)

