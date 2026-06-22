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

meta_16S$Tree_ID <- sub("^[0-9]+_", "", meta_16S$ID)
meta_16S$Tree_ID <- sub("_(Soil|Root|Leaf)$", "", meta_16S$Tree_ID)

meta_16S$Landuse <- factor(
  meta_16S$Landuse,
  levels = c("Forest", "UW", "Park", "Lawn", "Street")
)

#### 3. SOIL ONLY ####

idx_16S <- meta_16S$Type == "Soil"

meta_soil_16S <- meta_16S[idx_16S, ]
dist_soil_16S <- as.dist(
  as.matrix(dist_16S)[idx_16S, idx_16S]
)

#### 4. ENV ALIGNMENT ####

env_16S <- data_env[
  match(meta_soil_16S$Tree_ID, data_env$ID),
]

vars <- c(
  "Temperature",
  "Moisture",
  "pH",
  "SOM",
  "NH4_ug.g",
  "NO3_ug.g",
  "OlsenP_ug.g"
)

env_16S <- env_16S[, vars, drop = FALSE]

#### 5. REMOVE NA ####

keep_16S <- complete.cases(env_16S)

env_16S <- env_16S[keep_16S, ]
meta_soil_16S <- meta_soil_16S[keep_16S, ]

dist_soil_16S <- as.dist(
  as.matrix(dist_soil_16S)[keep_16S, keep_16S]
)

#### 6. SCALE ENV VARIABLES ####

env_16S <- as.data.frame(
  scale(env_16S)
)

#### 7. PERMANOVA ENVIRONMENTAL DRIVERS ####

formula_16S <- as.formula(
  paste(
    "dist_soil_16S ~",
    paste(colnames(env_16S), collapse = " + ")
  )
)

res_16S <- adonis2(
  formula_16S,
  data = env_16S,
  permutations = 999,
  by = "margin"
)

print(res_16S)

#### 8. BARPLOT R2 ####

df_perm <- as.data.frame(res_16S)

df_perm$Variable <- rownames(df_perm)

df_perm <- df_perm[
  df_perm$Variable %in% vars,
]

ggplot(
  df_perm,
  aes(
    x = Variable,
    y = R2
  )
) +
  geom_bar(
    stat = "identity",
    fill = "darkgreen"
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Environmental drivers of 16S soil communities",
    y = "R²",
    x = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

#### 9. db-RDA ####

rda_16S <- capscale(
  dist_soil_16S ~ .,
  data = env_16S
)

#### 10. TESTS ####

anova(rda_16S)

anova(
  rda_16S,
  by = "margin"
)

#### 11. VIF / COLLINEARITY ####

vif.cca(rda_16S)

#### 12. PLOT db-RDA ####

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

  df_sites$Landuse <- meta$Landuse

  df_env <- as.data.frame(env_scores)

  df_env$Variable <- rownames(df_env)

  var_exp <- summary(rda_model)$cont$importance[2, 1:2] * 100

  ggplot(
    df_sites,
    aes(
      x = CAP1,
      y = CAP2,
      color = Landuse
    )
  ) +
    geom_point(
      size = 2,
      alpha = 0.7
    ) +
    stat_ellipse(
      aes(group = Landuse),
      linewidth = 1
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
      arrow = arrow(
        length = unit(0.25, "cm")
      )
    ) +
    geom_text_repel(
      data = df_env,
      aes(
        x = CAP1,
        y = CAP2,
        label = Variable
      ),
      inherit.aes = FALSE
    ) +
    labs(
      title = title,
      x = paste0(
        "CAP1 (",
        round(var_exp[1], 1),
        "%)"
      ),
      y = paste0(
        "CAP2 (",
        round(var_exp[2], 1),
        "%)"
      )
    ) +
    theme_minimal(
      base_size = 14
    ) +
    coord_equal()
}

plot_rda(
  rda_16S,
  meta_soil_16S,
  "db-RDA - 16S Soil"
) 
