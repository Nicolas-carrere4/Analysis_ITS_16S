library(dplyr)
library(tidyr)
library(ggplot2)

#### 1. IMPORT ####

path <- "Data/"

guild_clr <- read.csv(
  paste0(path, "16S_functional_guilds_CLR_values.csv"),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

data_env <- read.csv(
  paste0(path, "Data_env.csv"),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

#### 2. SETTINGS ####

vars <- c(
  "Temperature",
  "Moisture",
  "pH",
  "SOM",
  "NH4_ug.g",
  "NO3_ug.g",
  "OlsenP_ug.g"
)

groups <- c(
  "Nitrification",
  "Denitrification",
  "N_fixation",
  "Assim_nitrite_reduction",
  "Dissim_nitrite_reduction",
  "Assim_nitrate_reduction",
  "Dissim_nitrate_reduction",
  "Cellulolytic",
  "Chitinolytic",
  "Lignolytic",
  "Methanotroph",
  "Copiotroph"
)

#### 3. KEEP SOIL ONLY ####

guild_soil <- guild_clr %>%
  dplyr::filter(Type == "Soil")

#### 4. CREATE TREE_ID ####

guild_soil$Tree_ID <- sub("^[0-9]+_", "", guild_soil$ID)
guild_soil$Tree_ID <- sub("_(Soil|Root|Leaf)$", "", guild_soil$Tree_ID)

#### 5. ALIGN ENVIRONMENTAL DATA ####

env_soil <- data_env[
  match(guild_soil$Tree_ID, data_env$ID),
  ,
  drop = FALSE
]

df <- cbind(
  guild_soil,
  env_soil[, vars, drop = FALSE]
)

#### 6. REMOVE NA ####

df <- df %>%
  dplyr::filter(complete.cases(dplyr::select(., all_of(vars), CLR)))

#### 7. CHECK ####

print(dim(df))
print(table(df$Function))
print(table(df$Landuse))

#### 8. SPEARMAN CORRELATIONS ####

res_cor_16S <- data.frame()

for(g in groups){

  d_g <- df %>%
    dplyr::filter(Function == g)

  if(nrow(d_g) < 10){
    message(paste("Skipped:", g, "- too few samples"))
    next
  }

  for(v in vars){

    test <- cor.test(
      d_g$CLR,
      d_g[[v]],
      method = "spearman",
      exact = FALSE
    )

    res_cor_16S <- rbind(
      res_cor_16S,
      data.frame(
        Function = g,
        Variable = v,
        Rho = unname(test$estimate),
        p.value = test$p.value,
        N = nrow(d_g)
      )
    )
  }
}

#### 9. MULTIPLE TEST CORRECTION ####

res_cor_16S$p.adj <- p.adjust(
  res_cor_16S$p.value,
  method = "BH"
)

#### 10. SIGNIFICANCE LABELS ####

res_cor_16S$Signif <- ""
res_cor_16S$Signif[res_cor_16S$p.adj < 0.05] <- "*"
res_cor_16S$Signif[res_cor_16S$p.adj < 0.01] <- "**"
res_cor_16S$Signif[res_cor_16S$p.adj < 0.001] <- "***"

res_cor_16S$Function <- factor(
  res_cor_16S$Function,
  levels = groups
)

res_cor_16S$Variable <- factor(
  res_cor_16S$Variable,
  levels = vars
)

print(res_cor_16S)

#### 11. HEATMAP ####

p_cor_16S <- ggplot(
  res_cor_16S,
  aes(
    x = Variable,
    y = Function,
    fill = Rho
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = Signif),
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Spearman correlations between 16S functional guild CLR abundance and soil variables",
    x = NULL,
    y = NULL,
    fill = "Spearman rho"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_cor_16S
