library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)

#### 1. IMPORT ####

path <- "Data/"

dist_16S <- readRDS(paste0(path, "Aitchison_distance_matrix_16S.rds"))
ps_16S <- readRDS(paste0(path, "phyloseq_no_negatives_16S.rds"))

data_env <- read.csv(
  paste0(path, "Data_env.csv"),
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

guild_clr <- read.csv(
  paste0(path, "16S_functional_guilds_CLR_values.csv"),
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

#### 2. METADATA ####

meta_16S <- data.frame(sample_data(ps_16S))
meta_16S$ID <- rownames(meta_16S)

meta_16S$Tree_ID <- sub("^[0-9]+_", "", meta_16S$ID)
meta_16S$Tree_ID <- sub("_(Soil|Root|Leaf)$", "", meta_16S$Tree_ID)

#### 3. SOIL ONLY ####

idx_16S <- meta_16S$Type == "Soil"

meta_soil_16S <- meta_16S[idx_16S, , drop = FALSE]

dist_soil_16S <- as.dist(
  as.matrix(dist_16S)[idx_16S, idx_16S]
)

#### 4. MANAGEMENT INTENSITY ####

mgmt <- c(
  Forest = 1,
  UW = 2,
  Park = 3,
  Lawn = 4,
  Street = 5
)

meta_soil_16S$Management <- mgmt[as.character(meta_soil_16S$Landuse)]

#### 5. ENVIRONMENTAL VARIABLES ####

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

#### 6. 16S FUNCTIONAL GUILDS CLR ####

guild_soil <- guild_clr %>%
  dplyr::filter(Type == "Soil") %>%
  dplyr::select(ID, Function, CLR)

guild_wide <- guild_soil %>%
  tidyr::pivot_wider(
    names_from = Function,
    values_from = CLR
  )

guild_wide <- as.data.frame(guild_wide)

rownames(guild_wide) <- guild_wide$ID
guild_wide$ID <- NULL

guild_wide <- guild_wide[
  meta_soil_16S$ID,
  ,
  drop = FALSE
]

#### 7. CLEAN COLUMN NAMES FOR ENVFIT ####

colnames(guild_wide) <- make.names(colnames(guild_wide))

#### 8. SCALE ENVIRONMENTAL VARIABLES ####

env_16S_scaled <- as.data.frame(scale(env_16S))

#### 9. NMDS 3D + ENVFIT FUNCTION ####

make_panel_3d <- function(dist_obj, meta, env, extra = NULL, title){

  fit_mat <- env

  if(!is.null(extra)){
    fit_mat <- cbind(fit_mat, extra)
  }

  keep <- complete.cases(fit_mat)

  fit_mat <- fit_mat[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  dist_obj <- as.dist(
    as.matrix(dist_obj)[keep, keep]
  )

  ord <- metaMDS(
    dist_obj,
    k = 3,
    trymax = 200,
    trace = FALSE,
    autotransform = FALSE
  )

  cat("\n====================================\n")
  cat(title, "\n")
  cat("Stress =", round(ord$stress, 3), "\n")
  cat("====================================\n")

  sites <- as.data.frame(scores(ord, display = "sites"))
  sites$Management <- meta$Management
  sites$Landuse <- factor(
    meta$Landuse,
    levels = c("Forest", "UW", "Park", "Lawn", "Street")
  )

  fit <- envfit(
    ord,
    fit_mat,
    permutations = 999
  )

  print(fit)

  vec <- as.data.frame(scores(fit, display = "vectors"))
  vec$Variable <- rownames(vec)
  vec$p <- fit$vectors$pvals
  vec$r2 <- fit$vectors$r

  vec_sig <- vec[
    vec$p <= 0.05,
    ,
    drop = FALSE
  ]

  if(nrow(vec_sig) > 0){

    max_site <- max(abs(c(sites$NMDS1, sites$NMDS2)), na.rm = TRUE)
    max_vec <- max(abs(c(vec_sig$NMDS1, vec_sig$NMDS2)), na.rm = TRUE)

    mult <- 0.85 * max_site / max_vec

    vec_sig$NMDS1 <- vec_sig$NMDS1 * mult
    vec_sig$NMDS2 <- vec_sig$NMDS2 * mult
  }

  p <- ggplot(
    sites,
    aes(
      x = NMDS1,
      y = NMDS2
    )
  ) +
    geom_point(
      aes(color = Management),
      size = 3,
      alpha = 0.85
    ) +
    scale_color_gradient(
      low = "#2e7d32",
      high = "#d32f2f"
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(title, " | 3D NMDS stress = ", round(ord$stress, 3)),
      color = "Management intensity"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(face = "bold")
    )

  if(nrow(vec_sig) > 0){

    p <- p +
      geom_segment(
        data = vec_sig,
        aes(
          x = 0,
          y = 0,
          xend = NMDS1,
          yend = NMDS2
        ),
        inherit.aes = FALSE,
        arrow = arrow(length = unit(0.22, "cm")),
        linewidth = 0.8
      ) +
      geom_text_repel(
        data = vec_sig,
        aes(
          x = NMDS1,
          y = NMDS2,
          label = Variable
        ),
        inherit.aes = FALSE,
        size = 4,
        max.overlaps = Inf
      )
  }

  return(
    list(
      ordination = ord,
      envfit = fit,
      plot = p,
      sites = sites,
      vectors_all = vec,
      vectors_sig = vec_sig
    )
  )
}

#### 10. RUN NMDS 16S ####

res_nmds_16S_3d <- make_panel_3d(
  dist_obj = dist_soil_16S,
  meta = meta_soil_16S,
  env = env_16S_scaled,
  extra = guild_wide,
  title = "NMDS - 16S soil bacterial communities"
)

res_nmds_16S_3d$plot
