# title: "03_exploratory_analysis"
# TUV-2025-sub/scripts/03_exploratory_analysis.R

# This script is to explore the sub annotation data
# and create initial figures and tables for the TUV report
# This hasn't quite yet been smoothed out to not be redundant with previous scripts

# ========== SETUP ==========================
library(tidyverse)
library(gt)

# ========== CONFIGURATION ==========================
ps_paths <- ps_science_paths()
exp_path <- file.path(ps_paths$expeditions, "TUV-2025")

# Load phylogenetic key locally from within the project, where the script is located
source("scripts/phylo_keys.R") 

# where data comes from and where it goes
raw_dir       <- file.path(exp_path, "data/primary/raw/sub")
processed_dir <- file.path(exp_path, "data/primary/processed/sub")
output_dir <- file.path(exp_path, "data/primary/output/sub")
annotations_dir <- file.path(raw_dir, "annotations")
out_path <- output_dir

# Load cleanly processed dataset (from 01_data_prep)
clean_taxa <- read_rds(file.path(processed_dir, "clean_master_data.rds"))

# Load transect video metadata
excel_file <- file.path(raw_dir, "TUV_2025_sub_fieldbook.xlsx")
vid_meta   <- read_excel(excel_file, sheet = "transects") %>%
  rename(
    dive_id = `ps_site_id`
  )

# ========= SAMPLING EFFORT BY DEPTH (HISTOGRAM) =========================================

p_effort <- vid_meta %>%
  filter(!is.na(start_depth)) %>%
  ggplot(aes(x = start_depth)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "black", alpha = 0.8) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = scales::breaks_width(1) # Integer ticks
  ) +
  labs(
    x = "Transect Start Depth (m)",
    y = "Number of Transects"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid = element_blank(), # No grid
    axis.line = element_line(color = "black"), # Solid baseline
    axis.ticks = element_line(color = "black", linewidth = 0.8), # Physical ticks restored
    axis.ticks.length = unit(0.25, "cm")
  )
print(p_effort)
# Save Histogram
ggsave(file.path(out_path, "sampling_effort_histogram.png"), plot = p_effort, width = 9, height = 6)


# ========= TAXA INVENTORY  ==========================================

## -------- Taxa inventory prep ---------
total_dives <- n_distinct(clean_taxa$dive_id)

taxa_inventory <- clean_taxa %>%
  group_by(phylum, class, order, family, scientificName, taxonRank) %>%
  summarize(
    count_obs = n(),
    min_depth = round(min(`depth(m)`, na.rm = TRUE), 0),
    max_depth = round(max(`depth(m)`, na.rm = TRUE), 0),
    taxa_deployments = n_distinct(dive_id), 
    .groups = "drop"
  ) %>%
  mutate(
    Depth_m = if_else(min_depth == max_depth, 
                      as.character(min_depth), 
                      paste0(min_depth, " - ", max_depth)),
    # round freq. of occ. to nearest whole #
    Freq_occ = round((taxa_deployments / total_dives) * 100, 0), 
    
    # add italics and " sp." conditionally
    # If the taxonRank is genus, species (or subspecies), wrap the name in markdown asterisks
    scientificName = case_when(
      # If Genus: italicize the name and add " sp." in regular font
      tolower(taxonRank) == "genus" ~ paste0("*", scientificName, "* sp."),
      # If Species/Subspecies/Morpho: just italicize the whole thing
      tolower(taxonRank) %in% c("species", "subspecies", "morphospecies") ~ paste0("*", scientificName, "*"),
      # Otherwise (Phylum, Family, etc.): leave it completely alone
      TRUE ~ scientificName
    ),
    
    # This taps the phylo_keys order, but tacks any unknown taxa onto the end 
    phylum = factor(phylum, levels = unique(c(phylo_phylum, na.omit(phylum)))),
    class  = factor(class,  levels = unique(c(phylo_class,  na.omit(class)))),
    order  = factor(order,  levels = unique(c(phylo_order,  na.omit(order)))),
    family = factor(family, levels = unique(c(phylo_family, na.omit(family))))
  ) %>%
  arrange(phylum, class, order, family, desc(count_obs)) %>%
  select(phylum, class, order, family, scientificName, Depth_m, Freq_occ, count_obs)


## ------- Full taxa inventory table (for appendix) ----------------

taxa_table_gt <- taxa_inventory %>%
  gt() %>%
  # replace NAs with blank space for the table
  sub_missing(
    columns = everything(),
    missing_text = ""
  ) %>%
  # figure labeling
   cols_label(
    phylum = "Phylum",
    class = "Class",
    order = "Order",
    family = "Family",
    scientificName = "Lowest Taxonomic ID",
    Depth_m = "Depth Range (m)",
    Freq_occ = "Freq. of Occ. (%)",
    count_obs = "Obs. Count"
  ) %>%
  # Alignments
  cols_align(align = "left", columns = c(phylum, class, order, family, scientificName)) %>%
  cols_align(align = "center", columns = c(Depth_m, Freq_occ, count_obs)) %>%
  # Strict CSS left-align to beat markdown/html defaults
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = c(phylum, class, order, family, scientificName))
  ) %>%
  # Widths & Styling
  cols_width(
    Depth_m ~ px(140),  
    Freq_occ ~ px(110),
    count_obs ~ px(90),  
    scientificName ~ px(220) 
  ) %>%
  fmt_markdown(columns = scientificName) %>%
  tab_style(
    style = cell_text(whitespace = "nowrap"),
    locations = cells_body(columns = Depth_m)
  ) %>%
  opt_row_striping() %>%
  tab_options(data_row.padding = px(5), table.font.size = px(12))

# Save Full Table
gtsave(taxa_table_gt, file.path(out_path, "Taxa_Inventory_Table.html"))


## --------- Top 15 tables (for main text) -------------

# 1. Get Top 15 Fish (Chordata)
top15_fish <- taxa_inventory %>%
  filter(phylum == "Chordata") %>%
  arrange(desc(count_obs)) %>%
  slice_head(n = 15) %>%
  mutate(fauna_group = "Fish") # Add a grouping column

# 2. Get Top 15 Invertebrates (Everything except Chordata)
top15_inverts <- taxa_inventory %>%
  filter(phylum != "Chordata") %>%
  arrange(desc(count_obs)) %>%
  slice_head(n = 15) %>%
  mutate(fauna_group = "Invertebrates")

# 3. Combine them back together
top15_combined <- bind_rows(top15_fish, top15_inverts)

# 4. Build the gt() table using the grouping column
top15_table_gt <- top15_combined %>%
  # gt() will automatically create beautiful divider rows based on this column
  gt(groupname_col = "fauna_group") %>%
  
  # Replace NA with blank space
  sub_missing(
    columns = everything(),
    missing_text = ""
  ) %>%
  
  cols_label(
    phylum = "Phylum", class = "Class", order = "Order", family = "Family",
    scientificName = "Lowest Taxonomic ID", Depth_m = "Depth Range (m)",
    Freq_occ = "Freq. of Occ. (%)", count_obs = "Obs. Count"
  ) %>%
  cols_align(align = "left", columns = c(phylum, class, order, family, scientificName)) %>%
  cols_align(align = "center", columns = c(Depth_m, Freq_occ, count_obs)) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = c(phylum, class, order, family, scientificName))
  ) %>%
  cols_width(
    Depth_m ~ px(140), Freq_occ ~ px(110),
    count_obs ~ px(90), scientificName ~ px(220) 
  ) %>%
  fmt_markdown(columns = scientificName) %>%
  tab_style(
    style = cell_text(whitespace = "nowrap"),
    locations = cells_body(columns = Depth_m)
  ) %>%
  
  # Add a light gray background to the group header rows
  tab_style(
    style = cell_fill(color = "gray90"),
    locations = cells_row_groups()
  ) %>%
  
  opt_row_striping() %>%
  tab_options(data_row.padding = px(5), table.font.size = px(12))

# Save the combined Top 15 Table
gtsave(top15_table_gt, file.path(out_path, "Top15_Fish_Inverts_Table.html"))


# ============= HABITAT ASSOCIATIONS (STACKED BAR CHART) =====================
# This displays a bar plot of the proportion of observations (in transect data), by phyla,
# observed upon which substrate types.

# Custom Colors
habitat_colors <- c(
  "boulder" = "gray",
  "cobble" = "gold",
  "in water column" = "lightblue",
  "rock" = "#800020",    # Burgundy
  "sand" = "#E2C290"     # Sand
)

top_phyla <- clean_taxa %>% count(phylum, sort = TRUE) %>% head(5) %>% pull(phylum)
top_upon <- clean_taxa %>% count(upon, sort = TRUE) %>% head(5) %>% pull(upon)

p_habitat <- clean_taxa %>%
  filter(phylum %in% top_phyla, upon %in% top_upon) %>%
  count(phylum, upon) %>%
  group_by(phylum) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = fct_infreq(phylum), y = pct, fill = upon)) +
  # width = 0.65 shrinks the bars to create a wider gap
  geom_col(position = "stack", color = "black", linewidth = 0.2, width = 0.65) + 
  scale_y_continuous(
    labels = function(x) x * 100, 
    expand = expansion(mult = c(0, 0.05)) # first 0 removes bottom gap, 0.05 gives gap at the top
  ) +
  scale_fill_manual(values = habitat_colors) +
  labs(
    x = "Phylum",
    y = "Relative Proportion of Observations (%)",
    fill = "Observed Upon"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    # Add margins to axis titles to create space
    # margin(t, r, b, l)
    axis.title.x = element_text(margin = margin(t = 15)), # Adds 15pts of space to Top of x-axis title
    axis.title.y = element_text(margin = margin(r = 15)), # Adds 15pts of space to Right of y-axis title
    # size = 10 shrinks the x-axis text to prevent squeezing
    axis.text.x = element_text(color = "black", size = 10), 
    axis.text.y = element_text(color = "black"),
    # Gridline cleanup
    panel.grid.major.x = element_blank(), # Removes main vertical gridlines
    panel.grid.minor.x = element_blank(), # Removes secondary vertical gridlines
    panel.grid.minor.y = element_blank()  # Removes the intermediate horizontal gridlines
  )

print(p_habitat)
# Save Habitat Plot
ggsave(file.path(out_path, "habitat_associations.png"), plot = p_habitat, width = 8, height = 5)

