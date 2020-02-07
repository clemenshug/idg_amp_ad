
library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("data"))

dir_data <- here("data")
dir.create(dir_data, showWarnings = FALSE)

# Wrangle AMP-AD data ----------------------------------------------------------
###############################################################################T

amp_ad_meta_raw <- tribble(
  ~dataset, ~fn,
  "rosmap", file.path(dir_data, "rosmap_meta.csv"),
  "msbb", file.path(dir_data, "msbb_meta.csv")
) %>%
  mutate(
    data = map(fn, read_csv)
  )

amp_ad_counts_raw <- tribble(
  ~dataset, ~synapse_id,
  "rosmap", "syn8691134",
  "msbb", "syn8691099"
) %>%
  mutate(
    data = map(synapse_id, syn) %>%
      map(read_tsv)
  )

# Assemble count and metadata --------------------------------------------------
###############################################################################T

braak_map <- set_names(
  c("A", "A", "A", "B", "B", "C", "C"),
  as.character(0:6)
)

meta <- amp_ad_meta_raw %>%
  select(dataset, data) %>%
  mutate(
    data = map(
      data,
      select, rnaseq_id, braak, batch
    ) %>%
      map(mutate_at, vars(braak), as.character)
  ) %>%
  unnest(data) %>%
  inner_join(enframe(braak_map, "braak", "stage"), by = "braak") %>%
  distinct()

counts <- amp_ad_counts_raw %>%
  pull(data) %>%
  reduce(inner_join, by = "feature") %>%
  filter(str_starts(feature, "ENSG"))

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

# library(biomaRt)
library(rtracklayer)

# Data was counted using Gencode v24
download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.basic.annotation.gtf.gz",
  file.path(dir_data, "gencode.v24.chr_patch_hapl_scaff.basic.annotation.gtf.gz")
)

gene_info <- GTFFile(file.path(dir_data, "gencode.v24.chr_patch_hapl_scaff.basic.annotation.gtf.gz")) %>%
  import() %>%
  as_tibble()

protein_coding <- gene_info %>%
  filter(gene_type == "protein_coding") %>%
  pull(gene_id) %>%
  unique()

meta_deseq <- meta %>%
  filter(rnaseq_id %in% colnames(counts)) %>%
  arrange(rnaseq_id) %>%
  column_to_rownames("rnaseq_id")

counts_deseq <- counts %>%
  filter(feature %in% protein_coding) %>%
  dplyr::select(feature, one_of(rownames(meta_deseq))) %>%
  column_to_rownames("feature") %>%
  as.matrix()

library(DESeq2)
unloadNamespace("synapser")
unloadNamespace("PythonEmbedInR")

deseq_dataset <- DESeqDataSetFromMatrix(
  counts_deseq, meta_deseq, design = ~ batch + stage
)

deseq_de <- DESeq(deseq_dataset)

write_rds(
  deseq_de,
  file.path(dir_data, "deseq.rds"),
  compress = "gz"
)

resultsNames(deseq_de)

res_c_vs_a <- results(deseq_de, name = "stage_C_vs_A")
res_c_vs_a_shrunken <- lfcShrink(deseq_de, coef = "stage_C_vs_A", type = "apeglm")

res_c_vs_a_df <- res_c_vs_a_shrunken %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  left_join(
    gene_info %>%
      distinct(gene_id, gene_name),
    by = "gene_id"
  ) %>%
  left_join(
    res_c_vs_a %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      select(gene_id, log2FoldChange_MLE = log2FoldChange, lfcSE_MLE = lfcSE),
    by = "gene_id"
  ) %>%
  select(gene_id, gene_name, everything()) %>%
  arrange(padj)

write_csv(
  res_c_vs_a_df,
  here("amp_ad_stage_C_vs_A.csv")
)

