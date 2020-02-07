library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("data"))

dir_data <- here("data")
dir.create(dir_data, showWarnings = FALSE)

# Wrangle MSBB data ------------------------------------------------------------
###############################################################################T

rnaseq_meta <- syn("syn6100548") %>%
  read_csv()

clinical_meta <- syn("syn6101474") %>%
  read_csv()

combined_meta <- rnaseq_meta %>%
  inner_join(
    clinical_meta %>%
      select(individualIdentifier, braak = bbscore, age_death = AOD),
    by = c("individualIdentifier.inferred" = "individualIdentifier")
  ) %>%
  mutate(
    rnaseq_id = sampleIdentifier,
    batch = paste0("MSBB_", BrodmannArea)
  )

write_csv(
  combined_meta,
  file.path(dir_data, "msbb_meta.csv")
)
