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

ps <- readRDS(
  paste0(path, "phyloseq_no_negatives_ITS.rds")
)

DATABASE <- read_excel(
  paste0(path, "Fungal_trait.xlsx")
)

DATABASE <- as.data.frame(DATABASE)

#### 2. TAXONOMY ####

tax <- as.data.frame(tax_table(ps))

tax$Genus <- gsub(
  "^g__",
  "",
  as.character(tax$Genus)
)

tax$ASV <- rownames(tax)

DATABASE$Genus <- gsub(
  "^g__",
  "",
  DATABASE[[grep(
    "GENUS|Genus",
    colnames(DATABASE),
    value = TRUE
  )[1]]]
)

lifestyle_col <- grep(
  "lifestyle",
  colnames(DATABASE),
  value = TRUE,
  ignore.case = TRUE
)[1]

#### 3. OTU TABLE ####

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

#### 4. METADATA ####

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

#### 5. ECM VS OTHER ####

merge_ft <- merge(
  tax,
  DATABASE,
  by = "Genus",
  all.x = TRUE
)

merge_ft <- merge_ft[
  match(rownames(ASV), merge_ft$ASV),
]

guild <- ifelse(
  grepl(
    "ecto",
    merge_ft[[lifestyle_col]],
    ignore.case = TRUE
  ),
  "ECM",
  "Other"
)

guild_levels <- c(
  "ECM",
  "Other"
)

guild_mat <- matrix(
  0,
  nrow = length(guild_levels),
  ncol = ncol(ASV)
)

rownames(guild_mat) <- guild_levels
colnames(guild_mat) <- colnames(ASV)

for(g in guild_levels){

  idx <- which(guild == g)

  if(length(idx) > 0){

    guild_mat[g, ] <- colSums(
      ASV[idx, , drop = FALSE]
    )
  }
}

storage.mode(guild_mat) <- "numeric"

#### 6. FUNCTION ####

run_type_ecm <- function(type_name){

  keep <- meta$ID[
    meta$Type == type_name
  ]

  samples <- intersect(
    colnames(guild_mat),
    keep
  )

  if(length(samples) < 2){
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

  #### KEEP ECM ONLY ####

  d <- df[
    df$Guild == "ECM",
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
    main = paste(type_name, "- ECM")
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
      alpha = 0.8,
      width = 0.6
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.5,
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
      size = 5,
      inherit.aes = FALSE,
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
      title = paste(type_name, "ECM"),
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

#### 7. RUN ####

res_leaf_ecm <- run_type_ecm("Leaf")
res_root_ecm <- run_type_ecm("Root")
res_soil_ecm <- run_type_ecm("Soil")

res_leaf_ecm$plot
res_root_ecm$plot
res_soil_ecm$plot
