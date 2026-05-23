library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)

#### 1. IMPORT ####

path <- "Data/"

ps_ITS <- readRDS(paste0(path, "phyloseq_CLR_normalized_ITS.rds"))
dist_ITS <- readRDS(paste0(path, "Aitchison_distance_matrix_ITS.rds"))

data_env <- read.csv(
  paste0(path, "Data_env.csv"),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

#### 2. METADATA ####

meta_ITS <- data.frame(sample_data(ps_ITS))

meta_ITS$Tree_ID <- sub("^[0-9]+_", "", meta_ITS$ID)
meta_ITS$Tree_ID <- sub("_(Soil|Root|Leaf)$", "", meta_ITS$Tree_ID)

#### 3. SOIL ONLY ####

idx_ITS <- meta_ITS$Type == "Soil"

meta_soil_ITS <- meta_ITS[idx_ITS, ]
dist_soil_ITS <- as.dist(as.matrix(dist_ITS)[idx_ITS, idx_ITS])

#### 4. ENV ALIGNMENT ####

env_ITS <- data_env[match(meta_soil_ITS$Tree_ID, data_env$ID), ]

vars <- c(
  "Temperature",
  "Moisture",
  "pH",
  "SOM",
  "NH4_ug.g",
  "NO3_ug.g",
  "OlsenP_ug.g"
)

env_ITS <- env_ITS[, vars, drop = FALSE]

#### 5. REMOVE NA ####

keep_ITS <- complete.cases(env_ITS)

env_ITS <- env_ITS[keep_ITS, ]
meta_soil_ITS <- meta_soil_ITS[keep_ITS, ]

dist_soil_ITS <- as.dist(
  as.matrix(dist_soil_ITS)[keep_ITS, keep_ITS]
)

#### 6. SCALE ENV VARIABLES ####

env_ITS <- as.data.frame(scale(env_ITS))

#### 7. PERMANOVA ENVIRONMENTAL DRIVERS ####

formula_ITS <- as.formula(
  paste("dist_soil_ITS ~", paste(colnames(env_ITS), collapse = " + "))
)

res_ITS <- adonis2(
  formula_ITS,
  data = env_ITS,
  permutations = 999,
  by = "margin"
)

print(res_ITS)

#### 8. BARPLOT R2 ####

df_perm <- as.data.frame(res_ITS)
df_perm$Variable <- rownames(df_perm)
df_perm <- df_perm[df_perm$Variable %in% vars, ]

ggplot(df_perm, aes(x = Variable, y = R2)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Environmental drivers of ITS soil communities",
    y = "R²",
    x = NULL
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#### 9. db-RDA ####

rda_ITS <- capscale(
  dist_soil_ITS ~ .,
  data = env_ITS
)

#### 10. TESTS ####

anova(rda_ITS)
anova(rda_ITS, by = "margin")

#### 11. VIF / COLLINEARITY ####

vif.cca(rda_ITS)

#### 12. PLOT db-RDA ####

plot_rda <- function(rda_model, meta, title){

  sites <- scores(rda_model, display = "sites")
  env_scores <- scores(rda_model, display = "bp")

  df_sites <- as.data.frame(sites)
  df_sites$Landuse <- meta$Landuse

  df_env <- as.data.frame(env_scores)
  df_env$Variable <- rownames(df_env)

  var_exp <- summary(rda_model)$cont$importance[2, 1:2] * 100

  ggplot(df_sites, aes(CAP1, CAP2, color = Landuse)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(aes(group = Landuse), linewidth = 1) +
    geom_segment(
      data = df_env,
      aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.25, "cm"))
    ) +
    geom_text_repel(
      data = df_env,
      aes(CAP1, CAP2, label = Variable),
      inherit.aes = FALSE
    ) +
    labs(
      title = title,
      x = paste0("CAP1 (", round(var_exp[1], 1), "%)"),
      y = paste0("CAP2 (", round(var_exp[2], 1), "%)")
    ) +
    theme_minimal(base_size = 14) +
    coord_equal()
}

plot_rda(rda_ITS, meta_soil_ITS, "db-RDA - ITS Soil") 
