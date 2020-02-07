
library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("data"))

dir_data <- here("data")
dir.create(dir_data, showWarnings = FALSE)

# Wrangle ROSMAP data ----------------------------------------------------------
###############################################################################T

biospecimen_meta <- syn("syn21323366") %>%
  read_csv()

clinical_meta <- syn("syn3191087") %>%
  read_csv()

rnaseq_meta <- syn("syn21088596") %>%
  read_csv()

combined_meta <- rnaseq_meta %>%
  inner_join(
    biospecimen_meta %>%
      select(specimenID, individualID),
    by = "specimenID"
  ) %>%
  drop_na(individualID) %>%
  inner_join(
    clinical_meta %>%
      select(individualID, Study, braak = braaksc, age_death),
    by = "individualID"
  ) %>%
  mutate(
    rnaseq_id = specimenID,
    batch = Study
  )

write_csv(
  combined_meta,
  file.path(dir_data, "rosmap_meta.csv")
)
