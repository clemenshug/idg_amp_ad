if (!require(devtools))
  install.packages("devtools")

# https://github.com/labsyspharm/DRIAD
devtools::install_github("labsyspharm/DRIAD")

library(tidyverse)
library(DRIAD)
library(here)

dir_data <- here("data")
dir.create(dir_data, showWarnings = FALSE)

# Wrangle AMP-AD data ----------------------------------------------------------
###############################################################################T

fn_rosmap <- wrangleROSMAP(dir_data)
fn_msbb <- wrangleMSBB(dir_data)

amp_ad_raw <- tribble(
  ~dataset, ~fn,
  "rosmap", fn_rosmap,
  "msbb10", fn_msbb[[1]],
  "msbb22", fn_msbb[[2]],
  "msbb36", fn_msbb[[3]],
  "msbb44", fn_msbb[[4]]
) %>%
  mutate(data = map(fn, read_tsv))

# Assemble count and metadata --------------------------------------------------
###############################################################################T

braak_map <- set_names(
  c("A", "A", "A", "B", "B", "C", "C"),
  as.character(0:6)
)

meta <- amp_ad_raw %>%
  select(dataset, data) %>%
  mutate(
    data = map(
      data,
      select, ID, Braak
    ) %>%
      map(mutate_at, vars(ID, Braak), as.character)
  ) %>%
  unnest(data) %>%
  mutate(ID = paste(dataset, ID, sep = "_")) %>%
  inner_join(enframe(braak_map, "Braak", "Stage"), by = "Braak")

counts <- amp_ad_raw %>%
  select(dataset, data) %>%
  mutate(
    data = map(
      data,
      ~select(.x, ID, which(colnames(.x) == "A1BG"):ncol(.x)) %>%
        mutate_at(vars(ID), as.character)
    )
  ) %>%
  unnest(data) %>%
  mutate(ID = paste(dataset, ID, sep = "_"))

write_csv(
  meta,
  file.path(dir_data, "meta_wangled.csv.gz")
)

write_csv(
  counts,
  file.path(dir_data, "counts_wangled.csv.gz")
)

# Differential expression ------------------------------------------------------
###############################################################################T

library(DESeq2)
library(biomaRt)
library(rtracklayer)

download.file(
  "ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz",
  file.path(dir_data, "Homo_sapiens.GRCh38.99.gtf.gz")
)

gene_info <- GTFFile(file.path(dir_data, "Homo_sapiens.GRCh38.99.gtf.gz")) %>%
  import() %>%
  as_tibble()

protein_coding <- gene_info %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(gene_name) %>%
  unique()

meta_deseq <- meta %>%
  filter(ID %in% counts$ID) %>%
  arrange(ID) %>%
  column_to_rownames("ID")

counts_deseq <- counts %>%
  filter(ID %in% rownames(meta_deseq)) %>%
  arrange(ID) %>%
  dplyr::select(ID, one_of(protein_coding[protein_coding %in% colnames(.)]))
%>%
  column_to_rownames("ID") %>%
  as.matrix() %>%
  t()

deseq_dataset <- DESeqDataSetFromMatrix(
  counts_deseq, meta_deseq, design = ~ dataset * Stage
)
