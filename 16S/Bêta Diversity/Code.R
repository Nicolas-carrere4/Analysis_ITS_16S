library(phyloseq)
library(vegan)
library(permute)
library(ggplot2)
library(multcompView)

#### 1. IMPORT ####

path <- "Data/"

ps_16S <- readRDS(
  paste0(path, "phyloseq_CLR_normalized_16S.rds")
)

dist_16S <- readRDS(
  paste0(path, "Aitchison_distance_matrix_16S.rds")
)

#### 2. METADATA ####

meta_16S <- data.frame(sample_data(ps_16S))
meta_16S$ID <- rownames(meta_16S)

meta_16S$Tree_ID <- sub("_[^_]+$", "", meta_16S$ID)
meta_16S$Tree_ID <- as.factor(meta_16S$Tree_ID)

meta_16S$Type <- factor(
  meta_16S$Type,
  levels = c("Leaf", "Root", "Soil")
)

meta_16S$Landuse <- factor(
  meta_16S$Landuse,
  levels = c("Forest", "UW", "Park", "Lawn", "Street")
)

#### 3. COLORS ####

cols_type <- c(
  "Leaf" = "#fb8072",
  "Root" = "#00ba38",
  "Soil" = "#619cff"
)

cols_landuse <- c(
  "Forest" = "#1b5e20",
  "UW"     = "#66bb6a",
  "Park"   = "#a5d6a7",
  "Lawn"   = "#dce775",
  "Street" = "#616161"
)

#### 4. GLOBAL PERMANOVA ####

perm_16S <- how(
  nperm = 999,
  blocks = meta_16S$Tree_ID
)

adonis2(
  dist_16S ~ Type,
  data = meta_16S,
  permutations = perm_16S
)

adonis2(
  dist_16S ~ Type * Landuse,
  data = meta_16S,
  permutations = perm_16S,
  by = "margin"
)

#### 5. PERMANOVA BY COMPARTMENT ####

run_permanova_by_type <- function(type_name){

  idx <- meta_16S$Type == type_name

  meta_sub <- meta_16S[idx, ]

  dist_sub <- as.dist(
    as.matrix(dist_16S)[idx, idx]
  )

  adonis2(
    dist_sub ~ Landuse,
    data = meta_sub,
    permutations = 999
  )
}

run_permanova_by_type("Leaf")
run_permanova_by_type("Root")
run_permanova_by_type("Soil")

#### 6. DISPERSION ####

run_dispersion <- function(group_var){

  group <- meta_16S[[group_var]]

  bd <- betadisper(dist_16S, group)

  print(anova(bd))
  print(permutest(bd, permutations = 999))
  print(TukeyHSD(bd))

  return(bd)
}

bd_type_16S <- run_dispersion("Type")
bd_landuse_16S <- run_dispersion("Landuse")

#### 7. GLOBAL PCoA ####

plot_pcoa_type <- function(){

  ord <- ordinate(
    ps_16S,
    method = "PCoA",
    distance = dist_16S
  )

  plot_ordination(ps_16S, ord, color = "Type") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = cols_type) +
    theme_minimal(base_size = 14) +
    labs(
      title = "PCoA - 16S - Type",
      color = "Compartment"
    )
}

plot_pcoa_type_landuse <- function(){

  ord <- ordinate(
    ps_16S,
    method = "PCoA",
    distance = dist_16S
  )

  plot_ordination(ps_16S, ord, color = "Type", shape = "Landuse") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = cols_type) +
    theme_minimal(base_size = 14) +
    labs(
      title = "PCoA - 16S - Type and Landuse",
      color = "Compartment",
      shape = "Landuse"
    )
}

plot_pcoa_type()
plot_pcoa_type_landuse()

#### 8. PCoA BY COMPARTMENT ####

plot_pcoa_by_type <- function(type_name){

  samples <- rownames(meta_16S)[meta_16S$Type == type_name]

  ps_sub <- prune_samples(samples, ps_16S)

  dist_sub <- as.dist(
    as.matrix(dist_16S)[samples, samples]
  )

  ord <- ordinate(
    ps_sub,
    method = "PCoA",
    distance = dist_sub
  )

  plot_ordination(ps_sub, ord, color = "Landuse") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = cols_landuse) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("PCoA - 16S -", type_name),
      color = "Landuse"
    )
}

plot_pcoa_by_type("Leaf")
plot_pcoa_by_type("Root")
plot_pcoa_by_type("Soil")

#### 9. DISPERSION PLOTS WITH LETTERS ####

plot_dispersion <- function(bd, group_var, title_name){

  tuk <- TukeyHSD(bd)$group

  pvals <- tuk[, "p adj"]
  names(pvals) <- rownames(tuk)

  letters_vec <- multcompView::multcompLetters(pvals)$Letters

  letters_df <- data.frame(
    Group = names(letters_vec),
    .group = letters_vec
  )

  df <- data.frame(
    Distance = bd$distances,
    Group = meta_16S[[group_var]]
  )

  df$Group <- factor(df$Group, levels = levels(meta_16S[[group_var]]))
  letters_df$Group <- factor(letters_df$Group, levels = levels(df$Group))

  ypos <- aggregate(Distance ~ Group, df, max)

  letters_df$y <- ypos$Distance[
    match(letters_df$Group, ypos$Group)
  ] + 0.1 * diff(range(df$Distance))

  p <- ggplot(df, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    geom_text(
      data = letters_df,
      aes(x = Group, y = y, label = .group),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = title_name,
      x = group_var,
      y = "Distance to centroid"
    )

  if(group_var == "Type"){
    p <- p + scale_fill_manual(values = cols_type)
  }

  if(group_var == "Landuse"){
    p <- p + scale_fill_manual(values = cols_landuse)
  }

  return(p)
}

plot_dispersion(
  bd_type_16S,
  "Type",
  "Dispersion - 16S - Type"
)

plot_dispersion(
  bd_landuse_16S,
  "Landuse",
  "Dispersion - 16S - Landuse"
)
