library(phyloseq)
library(dplyr)
library(ggplot2)
library(FSA)
library(rcompanion)

#### 1. IMPORT ####

path <- "Data/"

ps_16S_raw <- readRDS(paste0(path, "phyloseq_no_negatives_16S.rds"))

set.seed(123)

#### 2. SETTINGS ####

landuse_levels <- c("Forest", "UW", "Park", "Lawn", "Street")

cols <- c(
  "Forest" = "#1b5e20",
  "UW"     = "#66bb6a",
  "Park"   = "#a5d6a7",
  "Lawn"   = "#dce775",
  "Street" = "#616161"
)

#### 3. RAREFACTION BY TYPE ####

get_ps_by_type <- function(ps, type_name){

  ps_sub <- prune_samples(sample_data(ps)$Type == type_name, ps)

  ps_sub <- rarefy_even_depth(
    ps_sub,
    sample.size = min(sample_sums(ps_sub)),
    replace = FALSE,
    verbose = FALSE
  )

  return(ps_sub)
}

#### 4. SHANNON DATA ####

get_shannon_16S <- function(type_name){

  ps_sub <- get_ps_by_type(ps_16S_raw, type_name)

  alpha <- estimate_richness(ps_sub, measures = "Shannon")
  meta <- data.frame(sample_data(ps_sub))

  d <- cbind(alpha, meta)
  d$ID <- rownames(d)
  d$Type <- type_name

  d$Landuse <- factor(
    d$Landuse,
    levels = landuse_levels
  )

  return(d)
}

#### 5. KRUSKAL + DUNN + PLOT ####

plot_shannon_kruskal_16S <- function(type_name){

  d <- get_shannon_16S(type_name)

  d <- d %>%
    filter(!is.na(Shannon), !is.na(Landuse))

  print(type_name)

  kw <- kruskal.test(
    Shannon ~ Landuse,
    data = d
  )

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
    summarise(
      y = max(Shannon, na.rm = TRUE),
      .groups = "drop"
    )

  letters_df <- merge(
    letters_df,
    ypos,
    by = "Landuse"
  )

  letters_df$y <- letters_df$y +
    0.1 * diff(range(d$Shannon, na.rm = TRUE))

  letters_df$Landuse <- factor(
    letters_df$Landuse,
    levels = landuse_levels
  )

  p <- ggplot(
    d,
    aes(
      x = Landuse,
      y = Shannon,
      fill = Landuse
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.85,
      width = 0.65
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.5,
      size = 1.3,
      color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 3,
      color = "black"
    ) +
    geom_text(
      data = letters_df,
      aes(
        x = Landuse,
        y = y,
        label = .group
      ),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = cols) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = paste("16S Shannon diversity -", type_name),
      x = NULL,
      y = "Shannon diversity"
    )

  print(p)

  return(
    list(
      kruskal = kw,
      dunn = dunn,
      letters = letters_df,
      plot = p,
      data = d
    )
  )
}

#### 6. RUN ####

res_leaf_16S_shannon <- plot_shannon_kruskal_16S("Leaf")
res_root_16S_shannon <- plot_shannon_kruskal_16S("Root")
res_soil_16S_shannon <- plot_shannon_kruskal_16S("Soil")

res_leaf_16S_shannon$plot
res_root_16S_shannon$plot
res_soil_16S_shannon$plot
