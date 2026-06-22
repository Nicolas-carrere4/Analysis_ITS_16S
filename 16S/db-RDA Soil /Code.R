library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)

#### 1. IMPORT ####

path <- "Data/"

ps_16S <- readRDS(paste0(path, "phyloseq_CLR_normalized_16S.rds"))
dist_16S <- readRDS(paste0(path, "Aitchison_distance_matrix_16S.rds"))

data_env <- read.csv(
  paste0(path, "Data_env.csv"),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

#### 2. METADATA ####

meta_16S <- data.frame(sample_data(ps_16S))
meta_16S$ID <- rownames(meta_16S)

meta_16S$Tree_ID <- sub("^[0-9]+_", "", meta_16S$ID)
meta_16S$Tree_ID <- sub("_(Soil|Root|Leaf)$", "", meta_16S$Tree_ID)

#### 3. SOIL ONLY ####

idx_16S <- meta_16S$Type == "Soil"

meta_soil_16S <- meta_16S[idx_16S, ]

dist_soil_16S <- as.dist(
  as.matrix(dist_16S)[idx_16S, idx_16S]
)

#### 4. ENVIRONMENTAL ALIGNMENT ####

vars <- c(
  "Temperature",
  "Moisture",
  "pH",
  "SOM",
  "NH4_ug.g",
  "NO3_ug.g",
  "OlsenP_ug.g"
)

env_16S <- data_env[
  match(meta_soil_16S$Tree_ID, data_env$ID),
  vars,
  drop = FALSE
]

#### 5. REMOVE NA ####

keep_16S <- complete.cases(env_16S)

env_16S <- env_16S[keep_16S, , drop = FALSE]
meta_soil_16S <- meta_soil_16S[keep_16S, , drop = FALSE]

dist_soil_16S <- as.dist(
  as.matrix(dist_soil_16S)[keep_16S, keep_16S]
)

#### 6. CHECK ALIGNMENT ####

stopifnot(nrow(env_16S) == nrow(meta_soil_16S))
stopifnot(attr(dist_soil_16S, "Size") == nrow(meta_soil_16S))

print(dim(env_16S))
print(table(meta_soil_16S$Landuse))

#### 7. SCALE ENVIRONMENTAL VARIABLES ####

env_16S_scaled <- as.data.frame(scale(env_16S))

#### 8. PERMANOVA ENVIRONMENTAL DRIVERS ####

formula_16S <- as.formula(
  paste("dist_soil_16S ~", paste(colnames(env_16S_scaled), collapse = " + "))
)

res_16S_env <- adonis2(
  formula_16S,
  data = env_16S_scaled,
  permutations = 999,
  by = "margin"
)

print(res_16S_env)

#### 9. BARPLOT R2 ####

df_perm_16S <- as.data.frame(res_16S_env)
df_perm_16S$Variable <- rownames(df_perm_16S)

df_perm_16S <- df_perm_16S[
  df_perm_16S$Variable %in% vars,
]

p_r2_16S <- ggplot(
  df_perm_16S,
  aes(
    x = reorder(Variable, R2),
    y = R2
  )
) +
  geom_bar(
    stat = "identity",
    fill = "darkgreen"
  ) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Environmental drivers of 16S soil bacterial communities",
    x = NULL,
    y = "Marginal R²"
  )

p_r2_16S

#### 10. db-RDA ####

rda_16S <- capscale(
  dist_soil_16S ~ .,
  data = env_16S_scaled
)

#### 11. GLOBAL AND MARGINAL TESTS ####

rda_16S_global_test <- anova(
  rda_16S,
  permutations = 999
)

rda_16S_margin_test <- anova(
  rda_16S,
  by = "margin",
  permutations = 999
)

print(rda_16S_global_test)
print(rda_16S_margin_test)

#### 12. VIF / COLLINEARITY ####

vif_16S <- vif.cca(rda_16S)

print(vif_16S)

#### 13. PLOT db-RDA ####

plot_rda <- function(rda_model, meta, title){

  sites <- scores(
    rda_model,
    display = "sites"
  )

  env_scores <- scores(
    rda_model,
    display = "bp"
  )

  df_sites <- as.data.frame(sites)
  df_sites$Landuse <- factor(
    meta$Landuse,
    levels = c("Forest", "UW", "Park", "Lawn", "Street")
  )

  df_env <- as.data.frame(env_scores)
  df_env$Variable <- rownames(df_env)

  var_exp <- summary(rda_model)$cont$importance[2, 1:2] * 100

  p <- ggplot(
    df_sites,
    aes(
      x = CAP1,
      y = CAP2,
      color = Landuse
    )
  ) +
    geom_point(
      size = 2.3,
      alpha = 0.75
    ) +
    stat_ellipse(
      aes(group = Landuse),
      linewidth = 0.9,
      alpha = 0.8
    ) +
    geom_segment(
      data = df_env,
      aes(
        x = 0,
        y = 0,
        xend = CAP1,
        yend = CAP2
      ),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.25, "cm")),
      linewidth = 0.8
    ) +
    geom_text_repel(
      data = df_env,
      aes(
        x = CAP1,
        y = CAP2,
        label = Variable
      ),
      inherit.aes = FALSE,
      size = 4
    ) +
    scale_color_manual(
      values = c(
        "Forest" = "#1b5e20",
        "UW"     = "#66bb6a",
        "Park"   = "#a5d6a7",
        "Lawn"   = "#dce775",
        "Street" = "#616161"
      )
    ) +
    theme_minimal(base_size = 14) +
    coord_equal() +
    labs(
      title = title,
      x = paste0("CAP1 (", round(var_exp[1], 1), "%)"),
      y = paste0("CAP2 (", round(var_exp[2], 1), "%)"),
      color = "Land use"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(face = "bold")
    )

  return(p)
}

p_rda_16S <- plot_rda(
  rda_model = rda_16S,
  meta = meta_soil_16S,
  title = "db-RDA - 16S soil bacterial communities"
)

p_rda_16S
