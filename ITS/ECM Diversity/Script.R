library(phyloseq)
library(dplyr)
library(readxl)
library(ggplot2)
library(FSA)
library(rcompanion)

#### 1. IMPORT ####

path <- "Data/"

ps <- readRDS(paste0(path, "phyloseq_no_negatives_ITS.rds"))

DATABASE <- read_excel(paste0(path, "Fungal_trait.xlsx"))
DATABASE <- as.data.frame(DATABASE)

#### 2. TAXONOMY ####

tax <- as.data.frame(tax_table(ps))
tax$Genus <- gsub("^g__", "", as.character(tax$Genus))
tax$ASV <- rownames(tax)

DATABASE$Genus <- gsub(
  "^g__",
  "",
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
  grepl("ecto", merge_ft[[lifestyle_col]], ignore.case = TRUE)
]

#### 4. SETTINGS ####

types <- c("Leaf", "Root", "Soil")

landuse_levels <- c("Forest", "UW", "Park", "Lawn", "Street")

cols <- c(
  "Forest" = "#1b5e20",
  "UW"     = "#66bb6a",
  "Park"   = "#a5d6a7",
  "Lawn"   = "#dce775",
  "Street" = "#616161"
)

#### 5. FUNCTION : ECM SHANNON ####

get_ecm_shannon <- function(type_name){

  ps_sub <- prune_samples(sample_data(ps)$Type == type_name, ps)

  ps_sub <- prune_taxa(taxa_names(ps_sub) %in% ecm_asv, ps_sub)

  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  if(nsamples(ps_sub) < 3) return(NULL)

  set.seed(123)

  ps_sub <- rarefy_even_depth(
    ps_sub,
    sample.size = min(sample_sums(ps_sub)),
    replace = FALSE,
    verbose = FALSE
  )

  alpha <- estimate_richness(ps_sub, measures = "Shannon")
  meta <- data.frame(sample_data(ps_sub))

  d <- cbind(alpha, meta)
  d$ID <- rownames(d)

  d$Landuse <- factor(d$Landuse, levels = landuse_levels)

  return(d)
}

#### 6. FUNCTION : KRUSKAL + DUNN + PLOT ####

plot_ecm_shannon_kruskal <- function(type_name){

  d <- get_ecm_shannon(type_name)

  if(is.null(d)) return(NULL)

  d <- d %>%
    filter(!is.na(Shannon), !is.na(Landuse))

  print(type_name)

  kw <- kruskal.test(Shannon ~ Landuse, data = d)
  print(kw)

  dunn <- dunnTest(
    Shannon ~ Landuse,
    data = d,
    method = "bh"
  )

  print(dunn)

  letters_df <- cldList(
    P.adj ~ Comparison,
    data = dunn$res,
    threshold = 0.05
  )

  colnames(letters_df) <- c("Landuse", ".group")

  ypos <- d %>%
    group_by(Landuse) %>%
    summarise(y = max(Shannon, na.rm = TRUE), .groups = "drop")

  letters_df <- merge(letters_df, ypos, by = "Landuse")

  letters_df$y <- letters_df$y +
    0.1 * diff(range(d$Shannon, na.rm = TRUE))

  letters_df$Landuse <- factor(
    letters_df$Landuse,
    levels = landuse_levels
  )

  p <- ggplot(d, aes(x = Landuse, y = Shannon, fill = Landuse)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.65) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.3, color = "black") +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    geom_text(
      data = letters_df,
      aes(x = Landuse, y = y, label = .group),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = cols) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = paste(type_name, "ECM Shannon diversity"),
      x = NULL,
      y = "Shannon diversity"
    )

  return(p)
}

#### 7. RUN ####

p_leaf_ecm_shannon <- plot_ecm_shannon_kruskal("Leaf")
p_root_ecm_shannon <- plot_ecm_shannon_kruskal("Root")
p_soil_ecm_shannon <- plot_ecm_shannon_kruskal("Soil")

p_leaf_ecm_shannon
p_root_ecm_shannon
p_soil_ecm_shannon
