library(phyloseq)
library(vegan)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

#### 1. IMPORT ####

path <- "Data/"

dist_ITS <- readRDS(paste0(path, "Aitchison_distance_matrix_ITS.rds"))
ps_ITS <- readRDS(paste0(path, "phyloseq_no_negatives_ITS.rds"))

DATABASE <- read_excel(paste0(path, "Fungal_trait.xlsx"))
DATABASE <- as.data.frame(DATABASE)

data_env <- read.csv(
  paste0(path, "Data_env.csv"),
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

#### 2. METADATA ####

meta_ITS <- data.frame(sample_data(ps_ITS))
meta_ITS$ID <- rownames(meta_ITS)

meta_ITS$Tree_ID <- sub("^[0-9]+_", "", meta_ITS$ID)
meta_ITS$Tree_ID <- sub("_(Soil|Root|Leaf)$", "", meta_ITS$Tree_ID)

#### 3. SOIL ONLY ####

idx_ITS <- meta_ITS$Type == "Soil"

meta_ITS <- meta_ITS[idx_ITS, ]
dist_ITS <- as.dist(as.matrix(dist_ITS)[idx_ITS, idx_ITS])

#### 4. MANAGEMENT ####

mgmt <- c(
  Forest = 1,
  UW = 2,
  Park = 3,
  Lawn = 4,
  Street = 5
)

meta_ITS$Management <- mgmt[meta_ITS$Landuse]

#### 5. ENV ####

vars <- c(
  "Temperature", "Moisture", "pH", "SOM",
  "NH4_ug.g", "NO3_ug.g", "OlsenP_ug.g"
)

env_ITS <- data_env[match(meta_ITS$Tree_ID, data_env$ID), vars]

#### 6. ITS GUILDS ####

tax <- as.data.frame(tax_table(ps_ITS))
tax$ASV <- rownames(tax)
tax$Genus <- gsub("^g__", "", as.character(tax$Genus))

DATABASE$Genus <- gsub(
  "^g__",
  "",
  DATABASE[[grep("GENUS|Genus", colnames(DATABASE), value = TRUE)[1]]]
)

merge_ft <- merge(tax, DATABASE, by = "Genus", all.x = TRUE)

otu <- as(otu_table(ps_ITS), "matrix")
if(taxa_are_rows(ps_ITS)) otu <- t(otu)

otu <- otu[idx_ITS, , drop = FALSE]

guild_list <- list(
  ECM = merge_ft$ASV[
    merge_ft$primary_lifestyle == "ectomycorrhizal" |
      merge_ft$Secondary_lifestyle == "ectomycorrhizal"
  ],
  Plant_pathogen = merge_ft$ASV[
    merge_ft$primary_lifestyle == "plant_pathogen" |
      merge_ft$Secondary_lifestyle == "plant_pathogen"
  ],
  Animal_parasite = merge_ft$ASV[
    merge_ft$primary_lifestyle == "animal_parasite" |
      merge_ft$Secondary_lifestyle == "animal_parasite"
  ],
  Saprotroph = merge_ft$ASV[
    grepl("saprotroph", merge_ft$primary_lifestyle, ignore.case = TRUE) |
      grepl("saprotroph", merge_ft$Secondary_lifestyle, ignore.case = TRUE)
  ]
)

guild_df <- data.frame(row.names = rownames(otu))

for(g in names(guild_list)){
  keep_asv <- intersect(colnames(otu), guild_list[[g]])

  if(length(keep_asv) == 0){
    guild_df[[g]] <- rep(0, nrow(otu))
  } else {
    guild_df[[g]] <- rowSums(otu[, keep_asv, drop = FALSE])
  }
}

guild_df <- guild_df + 1

clr_fun <- function(x){
  log(x / exp(mean(log(x))))
}

guild_df <- as.data.frame(t(apply(guild_df, 1, clr_fun)))

#### 7. NMDS 3D + ENVFIT ####

make_panel_3d <- function(dist_obj, meta, env, extra = NULL, title){

  fit_mat <- env

  if(!is.null(extra)){
    fit_mat <- cbind(fit_mat, extra)
  }

  keep <- complete.cases(fit_mat)

  fit_mat <- fit_mat[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  dist_obj <- as.dist(as.matrix(dist_obj)[keep, keep])

  ord <- metaMDS(
    dist_obj,
    k = 3,
    trymax = 200,
    trace = FALSE,
    autotransform = FALSE
  )

  print(paste(title, "- Stress =", round(ord$stress, 3)))

  sites <- as.data.frame(scores(ord, display = "sites"))
  sites$Management <- meta$Management
  sites$Landuse <- meta$Landuse

  fit <- envfit(ord, fit_mat, permutations = 999)

  vec <- as.data.frame(scores(fit, display = "vectors"))
  vec$Variable <- rownames(vec)
  vec$p <- fit$vectors$pvals

  vec <- vec[vec$p <= 0.05, , drop = FALSE]

  if(nrow(vec) > 0){
    max_site <- max(abs(c(sites$NMDS1, sites$NMDS2)))
    max_vec  <- max(abs(c(vec$NMDS1, vec$NMDS2)))

    mult <- 0.85 * max_site / max_vec

    vec$NMDS1 <- vec$NMDS1 * mult
    vec$NMDS2 <- vec$NMDS2 * mult
  }

  p <- ggplot(sites, aes(NMDS1, NMDS2)) +
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

  if(nrow(vec) > 0){
    p <- p +
      geom_segment(
        data = vec,
        aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
        inherit.aes = FALSE,
        arrow = arrow(length = unit(0.22, "cm")),
        linewidth = 0.8
      ) +
      geom_text_repel(
        data = vec,
        aes(NMDS1, NMDS2, label = Variable),
        inherit.aes = FALSE,
        size = 4
      )
  }

  return(list(
    ordination = ord,
    envfit = fit,
    plot = p,
    sites = sites,
    vectors = vec
  ))
}

#### 8. RUN ####

res_nmds_ITS_3d <- make_panel_3d(
  dist_obj = dist_ITS,
  meta = meta_ITS,
  env = env_ITS,
  extra = guild_df,
  title = "NMDS - ITS soil fungal communities"
)

res_nmds_ITS_3d$plot
