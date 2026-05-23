library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(scales)

#### 1. IMPORT ####

path <- "Data/"

ps <- readRDS(paste0(path, "phyloseq_no_negatives_ITS.rds"))

DATABASE <- read_excel(paste0(path, "Fungal_trait.xlsx"))
DATABASE <- as.data.frame(DATABASE)

#### 2. TAX ####

tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)
tax$Genus <- gsub("^g__", "", as.character(tax$Genus))

DATABASE$Genus <- gsub(
  "^g__", "",
  DATABASE[[grep("GENUS|Genus", colnames(DATABASE), value = TRUE)[1]]]
)

lifestyle_col <- grep(
  "lifestyle",
  colnames(DATABASE),
  value = TRUE,
  ignore.case = TRUE
)[1]

merge_ft <- merge(tax, DATABASE, by = "Genus", all.x = TRUE)

#### 3. ECM ASVs ####

ecm_asv <- merge_ft$ASV[
  grepl("ecto", merge_ft[[lifestyle_col]], ignore.case = TRUE) &
    grepl("mycorrh", merge_ft[[lifestyle_col]], ignore.case = TRUE)
]

#### 4. OTU TABLE ####

otu <- as(otu_table(ps), "matrix")

if(taxa_are_rows(ps)){
  otu <- t(otu)
}

#### 5. ECM ONLY ####

ASV_ecm <- otu[, colnames(otu) %in% ecm_asv, drop = FALSE]

#### 6. RELATIVE ABUNDANCE WITHIN ECM COMMUNITY ####

rel <- ASV_ecm / rowSums(ASV_ecm)
rel[is.na(rel)] <- 0

#### 7. GENUS AGGREGATION ####

genus_vec <- merge_ft$Genus[match(colnames(rel), merge_ft$ASV)]

rel_genus <- rowsum(t(rel), group = genus_vec)
rel_genus <- t(rel_genus)

#### 8. METADATA ####

meta <- data.frame(sample_data(ps))
meta$ID <- rownames(meta)

meta$Landuse <- factor(
  meta$Landuse,
  levels = c("Forest", "UW", "Park", "Lawn", "Street")
)

#### 9. DATAFRAME ####

df <- as.data.frame(rel_genus)
df$ID <- rownames(df)

df <- merge(df, meta, by = "ID")

#### 10. TOP GENERA ####

top_genus <- colMeans(rel_genus) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(8)

#### 11. LONG FORMAT ####

df_long <- df %>%
  pivot_longer(
    cols = colnames(rel_genus),
    names_to = "Genus",
    values_to = "RelAbund"
  ) %>%
  mutate(
    Genus = ifelse(Genus %in% top_genus, Genus, "Other")
  ) %>%
  group_by(ID, Type, Landuse, Neighborhood, Genus) %>%
  summarise(
    RelAbund = sum(RelAbund),
    .groups = "drop"
  )

#### 12. MEAN ABUNDANCE FOR PLOTS ####

df_mean <- df_long %>%
  group_by(Landuse, Genus) %>%
  summarise(
    mean_abund = mean(RelAbund),
    .groups = "drop"
  )

#### 13. ORDER GENERA ####

genus_order <- df_mean %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(mean_abund), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  pull(Genus)

df_mean$Genus <- factor(df_mean$Genus, levels = genus_order)

#### 14. COLORS ####

genus_cols <- hue_pal()(length(levels(df_mean$Genus)))
names(genus_cols) <- levels(df_mean$Genus)

#### 15. STACKED BARPLOT ####

p_bar <- ggplot(df_mean, aes(x = Landuse, y = mean_abund, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_cols) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Dominant ECM genera across land-use gradient",
    y = "Mean relative abundance within ECM community",
    x = NULL
  )

#### 16. HEATMAP ####

p_heatmap <- ggplot(df_mean, aes(x = Landuse, y = Genus, fill = mean_abund)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal(base_size = 14) +
  labs(
    title = "ECM genus distribution across land use",
    x = NULL,
    y = NULL,
    fill = "Mean relative abundance"
  )

p_bar
p_heatmap
