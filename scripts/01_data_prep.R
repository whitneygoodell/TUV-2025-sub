# title: "01_data_prep"
# scripts/prep/01_data_prep.R

# This script pulls submersible annotated data from the PS Drive and cleans and preps it
# for data analysis

# ==== SETUP =============================================
library(tidyverse)
library(readxl)
library(PristineSeasR2)

# ==== CONFIGURATION =============================================
ps_paths <- ps_science_paths()
exp_path <- file.path(ps_paths$expeditions, "TUV-2025")

# where data comes from and where it goes
raw_dir       <- file.path(exp_path, "data/primary/raw/sub")
processed_dir <- file.path(exp_path, "data/primary/processed/sub")
annotations_dir <- file.path(raw_dir, "annotations")


# ==== LOAD DATA =========
excel_file <- file.path(annotations_dir, "TUV_DOEX0112_sub_Data_Summary_v3.xlsx")
tator_data <- read_excel(excel_file, sheet = "Tator_Data_Summary")%>%
  rename(
    dive_id = deployment
  )
metadata   <- read_excel(excel_file, sheet = "Field_Log_Metadata") %>%
  rename(
    dive_id = `Site ID`,
    filename = Filename,
    avg_depth = `Average Depth (m)`,
    region = Region,
    island = Subregion
  )

# ==== MERGE AND CLEAN =========

# Merge metadata with annotation data
master_data <- tator_data %>%
  left_join(
    metadata %>%
      select(
        dive_id,
        region,
        island
        ),
    by = "dive_id"
  )

# ====  MANUAL EXCLUSIONS ========

## In dropcam data, this section is used for cleaning out taxa exclusions
# such as for taxa observed at the surfaces,
# of Deployment exclusions (For freq. of occ. calculations), flagging Partial Deployments.

# Don't yet know what kinds of exclusions the sub workflow may need, but leaving this section here in case.

# Exclusions are applied to the master_data so that they are not incorporated in downstream analyses.

# ======= EXPORT PROCESSED DATA =================================

# Just in case the processed folder doesn't exist on the Drive yet, this creates it:
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# Save clean dataset using processed directory defined in configuration
output_file <- file.path(processed_dir, "clean_master_data.rds")
write_rds(master_data, output_file)

message(paste("Data prep complete! Clean data saved to:", output_file))