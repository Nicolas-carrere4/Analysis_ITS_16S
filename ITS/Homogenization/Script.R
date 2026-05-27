library(phyloseq)
library(vegan)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ggplot2)
library(DHARMa)
library(performance)

#### 1. IMPORT ####

path <- "Data/"

ps_ITS <- readRDS(paste0(path, "phyloseq_CLR_normalized_ITS.rds"))
dist_ITS <- readRDS(paste0(path, "Aitchison_distance_matrix_ITS.rds"))

#### 2. SETTINGS ####

landuse_levels <- c("Forest", "UW", "Park", "Lawn", "Street")

cols <- c(
  "Forest" = "#1b5e20",
  "UW"     = "#66bb6a",
  "Park"   = "#a5d6a7",
  "Lawn"   = "#dce775",
  "Street" = "#616161"
)

#### 3. FUNCTION ####

run_homogenization_lmm_ITS <- function(type_name){

  meta <- data.frame(sample_data(ps_ITS))
  meta$ID <- rownames(meta)

  meta$Type <- factor(meta$Type, levels = c("Leaf", "Root", "Soil"))
  meta$Landuse <- factor(meta$Landuse, levels = landuse_levels)
  meta$Neighborhood <- factor(meta$Neighborhood)

  samples <- rownames(meta)[meta$Type == type_name]
  samples <- intersect(samples, labels(dist_ITS))

  meta_sub <- meta[match(samples, rownames(meta)), ]

  meta_sub$Landuse <- droplevels(meta_sub$Landuse)
  meta_sub$Neighborhood <- droplevels(meta_sub$Neighborhood)

  dist_sub <- as.dist(as.matrix(dist_ITS)[samples, samples])

  bd <- betadisper(
    dist_sub,
    meta_sub$Landuse
  )

  d <- data.frame(
    ID = samples,
    Distance = bd$distances
  )

  d <- merge(d, meta_sub, by = "ID")

  d$Landuse <- factor(d$Landuse, levels = landuse_levels)
  d$Landuse <- droplevels(d$Landuse)
  d$Neighborhood <- droplevels(factor(d$Neighborhood))

  mod <- lmer(
    Distance ~ Landuse + (1 | Neighborhood),
    data = d
  )

  print(type_name)

  print(check_singularity(mod))
  print(anova(mod))

  sim_res <- DHARMa::simulateResiduals(
    fittedModel = mod,
    n = 1000
  )

  plot(sim_res)

  DHARMa::plotQQunif(
    sim_res,
    main = paste(type_name, "- QQ plot")
  )

  DHARMa::plotResiduals(
    sim_res,
    form = d$Landuse,
    main = paste(type_name, "- Residuals vs Landuse")
  )

  print(DHARMa::testUniformity(sim_res))
  print(DHARMa::testDispersion(sim_res))
  print(DHARMa::testOutliers(sim_res))

  print(performance::check_model(mod))

  emm <- emmeans(
    mod,
    ~ Landuse
  )

  print(
    pairs(
      emm,
      adjust = "tukey"
    )
  )

  cld <- multcomp::cld(
    emm,
    Letters = letters,
    adjust = "tukey"
  )

  cld <- as.data.frame(cld)
  cld$.group <- gsub(" ", "", cld$.group)

  ypos <- aggregate(
    Distance ~ Landuse,
    d,
    max
  )

  cld <- merge(cld, ypos, by = "Landuse")

  cld$y <- cld$Distance +
    0.12 * diff(range(d$Distance))

  p <- ggplot(
    d,
    aes(
      x = Landuse,
      y = Distance,
      fill = Landuse
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.75
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.45,
      size = 1.2,
      color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 3,
      color = "black"
    ) +
    geom_text(
      data = cld,
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
      title = paste("Homogenization LMM - ITS -", type_name),
      x = "Land use",
      y = "Distance to centroid"
    )

  print(p)

  return(
    list(
      betadisper = bd,
      model = mod,
      residuals = sim_res,
      emmeans = emm,
      letters = cld,
      data = d,
      plot = p
    )
  )
}

#### 4. RUN ####

hom_leaf_ITS <- run_homogenization_lmm_ITS("Leaf")
hom_root_ITS <- run_homogenization_lmm_ITS("Root")
hom_soil_ITS <- run_homogenization_lmm_ITS("Soil")
