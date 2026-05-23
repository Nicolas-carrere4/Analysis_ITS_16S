library(phyloseq)
library(dplyr)
library(tidyr)
library(readxl)
library(ALDEx2)
library(lme4)
library(emmeans)
library(multcomp)
library(ggplot2)
library(DHARMa)
library(performance)

#### 1. IMPORT ####

path <- "Data/"

ps <- readRDS(paste0(path, "phyloseq_no_negatives_ITS.rds"))

DATABASE <- read_excel(paste0(path, "Fungal_trait.xlsx"))
DATABASE <- as.data.frame(DATABASE)

#### 2. TAXONOMY ####

tax <- as.data.frame(tax_table(ps))

tax$ASV <- rownames(tax)

tax$Genus <- gsub(
  "^g__",
  "",
  as.character(tax$Genus)
)

DATABASE$Genus <- gsub(
  "^g__",
  "",
  DATABASE[[grep(
    "GENUS|Genus",
    colnames(DATABASE),
    value = TRUE
  )[1]]]
)

merge_ft <- merge(
  tax,
  DATABASE,
  by = "Genus",
  all.x = TRUE
)

#### 3. ANIMAL PARASITE ASVs ####

primary_col <- "primary_lifestyle"
secondary_col <- "Secondary_lifestyle"

animal_rows <- (
  merge_ft[[primary_col]] == "animal_parasite" |
  merge_ft[[secondary_col]] == "animal_parasite"
)

animal_rows[is.na(animal_rows)] <- FALSE

animal_asv <- merge_ft$ASV[animal_rows]

#### 4. OTU TABLE ####

ASV <- as(otu_table(ps), "matrix")

if(!taxa_are_rows(ps)){
  ASV <- t(ASV)
}

ASV <- matrix(
  as.numeric(ASV),
  nrow = nrow(ASV)
)

rownames(ASV) <- taxa_names(ps)
colnames(ASV) <- sample_names(ps)

#### 5. METADATA ####

meta <- data.frame(sample_data(ps))

meta$ID <- rownames(meta)

meta$Neighborhood <- as.factor(
  meta$Neighborhood
)

landuse_levels <- c(
  "Forest",
  "UW",
  "Park",
  "Lawn",
  "Street"
)

cols <- c(
  "Forest" = "#1b5e20",
  "UW"     = "#66bb6a",
  "Park"   = "#a5d6a7",
  "Lawn"   = "#dce775",
  "Street" = "#616161"
)

#### 6. ANIMAL PATHOGENS VS OTHER ####

animal_idx <- rownames(ASV) %in% animal_asv
other_idx  <- !animal_idx

guild_mat <- rbind(
  Animal_pathogens = colSums(
    ASV[animal_idx, , drop = FALSE]
  ),
  Other = colSums(
    ASV[other_idx, , drop = FALSE]
  )
)

storage.mode(guild_mat) <- "numeric"

#### 7. FUNCTION ####

run_type_animal_clr <- function(type_name){

  keep <- meta$ID[
    meta$Type == type_name
  ]

  samples <- intersect(
    colnames(guild_mat),
    keep
  )

  if(length(samples) < 3){
    return(NULL)
  }

  mat <- guild_mat[
    ,
    samples,
    drop = FALSE
  ]

  mat <- mat[
    rowSums(mat) > 0,
    ,
    drop = FALSE
  ]

  if(nrow(mat) < 2){
    return(NULL)
  }

  meta_sub <- meta[
    match(samples, meta$ID),
  ]

  cond <- as.character(
    meta_sub$Landuse
  )

  #### CLR ####

  ald <- aldex.clr(
    mat,
    cond,
    mc.samples = 128
  )

  mc <- ALDEx2::getMonteCarloInstances(
    ald
  )

  clr <- do.call(
    cbind,
    lapply(
      mc,
      function(x){

        apply(
          as.matrix(x),
          1,
          median
        )
      }
    )
  )

  colnames(clr) <- names(mc)

  #### DATAFRAME ####

  df <- as.data.frame(clr)

  df$Guild <- rownames(df)

  df <- pivot_longer(
    df,
    -Guild,
    names_to = "ID",
    values_to = "CLR"
  )

  df <- merge(
    df,
    meta_sub,
    by = "ID"
  )

  #### KEEP ANIMAL PATHOGENS ####

  d <- df[
    df$Guild == "Animal_pathogens",
  ]

  if(nrow(d) < 3){
    return(NULL)
  }

  d$Landuse <- factor(
    d$Landuse,
    levels = landuse_levels
  )

  #### MODEL ####

  mod <- lmer(
    CLR ~ Landuse + (1 | Neighborhood),
    data = d
  )

  #### DIAGNOSTICS ####

  check_singularity(mod)

  sim_res <- simulateResiduals(mod)

  plotQQunif(
    sim_res,
    main = paste(type_name, "- Animal pathogens")
  )

  testDispersion(sim_res)
  testOutliers(sim_res)

  #### EMMEANS ####

  emm <- emmeans(
    mod,
    ~ Landuse
  )

  print(type_name)

  print(
    pairs(
      emm,
      adjust = "tukey"
    )
  )

  #### LETTERS ####

  cld <- multcomp::cld(
    emm,
    Letters = letters,
    adjust = "tukey"
  )

  cld <- as.data.frame(cld)

  cld$.group <- gsub(
    " ",
    "",
    cld$.group
  )

  #### POSITION ####

  ypos <- aggregate(
    CLR ~ Landuse,
    d,
    max
  )

  cld <- merge(
    cld,
    ypos,
    by = "Landuse"
  )

  cld$y <- cld$CLR + 0.3

  #### PLOT ####

  p <- ggplot(
    d,
    aes(
      x = Landuse,
      y = CLR,
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
    scale_fill_manual(
      values = cols
    ) +
    theme_minimal(
      base_size = 14
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = paste(
        type_name,
        "Animal pathogens abundance"
      ),
      x = NULL,
      y = "CLR abundance"
    )

  return(
    list(
      model = mod,
      emmeans = emm,
      letters = cld,
      plot = p,
      data = d
    )
  )
}

#### 8. RUN ####

res_leaf_animal <- run_type_animal_clr("Leaf")
res_root_animal <- run_type_animal_clr("Root")
res_soil_animal <- run_type_animal_clr("Soil")

res_leaf_animal$plot
res_root_animal$plot
res_soil_animal$plot
