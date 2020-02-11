
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
      select, rnaseq_id, braak, batch, age_death
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
  file.path(dir_data, "meta_wrangled.csv.gz")
)

write_csv(
  counts,
  file.path(dir_data, "counts_wrangled.csv.gz")
)

# Differential expression ------------------------------------------------------
###############################################################################T

# library(biomaRt)
library(rtracklayer)

# Data was counted using Gencode v24
if (!file.exists(file.path(dir_data, "gencode.v24.chr_patch_hapl_scaff.basic.annotation.gtf.gz")))
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
  mutate(
    age_death = age_death %>%
      recode(`90+` = "90") %>%
      as.double(),
    age_death_scaled = scale(age_death)[, 1]
  ) %>%
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
  counts_deseq, meta_deseq, design = ~ age_death_scaled + batch + stage
)

deseq_de <- DESeq(deseq_dataset)

write_rds(
  deseq_de,
  file.path(dir_data, "deseq_age.rds"),
  compress = "gz"
)

extract_result <- function(de, name) {
  res <- results(de, name = name, parallel = FALSE)
  shrunken <- lfcShrink(de, coef = name, res = res, type = "apeglm", parallel = FALSE)
  shrunken %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    left_join(
      gene_info %>%
        distinct(gene_id, gene_name),
      by = "gene_id"
    ) %>%
    left_join(
      res %>%
        as.data.frame() %>%
        rownames_to_column("gene_id") %>%
        select(gene_id, log2FoldChange_MLE = log2FoldChange, lfcSE_MLE = lfcSE),
      by = "gene_id"
    ) %>%
    select(gene_id, gene_name, everything()) %>%
    arrange(padj)
}

resultsNames(deseq_de)

res_c_vs_a <- extract_result(deseq_de, name = "stage_C_vs_A")
res_age <- extract_result(deseq_de, name = "age_death_scaled")

write_csv(
  res_c_vs_a,
  here("amp_ad_stage_C_vs_A_age_adjusted.csv")
)

bind_rows(
  "c_vs_a" = res_c_vs_a,
  "age_death" = res_age,
  .id = "comparison"
) %>%
  ggplot(aes(baseMean, log2FoldChange, color = padj < 0.05)) +
    geom_point() +
    scale_x_log10() +
    facet_wrap(vars(comparison))

volcano_plot <- bind_rows(
  "c_vs_a" = res_c_vs_a,
  "age_death" = res_age,
  .id = "comparison"
) %>%
  ggplot(aes(log2FoldChange, -log10(padj))) +
  geom_point() +
  facet_wrap(vars(comparison))

ggsave(
  file.path(dir_plots, "volcano_age_c_vs_a.pdf"),
  volcano_plot
)

# Age matching -----------------------------------------------------------------
###############################################################################T

library(MatchIt)

meta_matching <- meta_deseq %>%
  rownames_to_column("sample") %>%
  mutate_at(vars(stage, batch), as.factor) %>%
  as_tibble() %>%
  filter(stage %in% c("A", "B")) %>%
  mutate(stage = if_else(stage == "A", 1L, 0L)) %>%
  select(stage, batch, age_death)

meta_matched <- matchit(
  stage ~ age_death,
  data = meta_matching,
  method = "nearest"
)

summary(meta_matched)
## Ugh, somehow makes difference in ages even worse...
