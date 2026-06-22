library(phyloseq)
library(dplyr)
library(tidyr)
library(ALDEx2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcompView)
library(ggplot2)
library(DHARMa)
library(performance)

#### 1. IMPORT ####

path <- "Data/"

ps_16S <- readRDS(paste0(path, "phyloseq_no_negatives_16S.rds"))
func_db <- readRDS(paste0(path, "bacteria_tax_to_function.rds"))

#### 2. SETTINGS ####

landuse_levels <- c("Forest", "UW", "Park", "Lawn", "Street")
type_levels <- c("Leaf", "Root", "Soil")

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

cols <- c(
  "Forest" = "#1b5e20",
  "UW"     = "#66bb6a",
  "Park"   = "#a5d6a7",
  "Lawn"   = "#dce775",
  "Street" = "#616161"
)

mc_n <- 64

#### 3. TAXONOMY + FUNCTIONAL DATABASE ####

tax <- as.data.frame(tax_table(ps_16S))
tax$ASV <- rownames(tax)
tax$Genus <- trimws(gsub("^g__", "", as.character(tax$Genus)))

func_db$Taxon <- trimws(gsub("^g__", "", as.character(func_db$Taxon)))

func_db <- func_db[, c("Taxon", groups), drop = FALSE]
func_db[is.na(func_db)] <- 0

for(g in groups){
  func_db[[g]] <- as.numeric(func_db[[g]])
}

tax_func <- merge(
  tax,
  func_db,
  by.x = "Genus",
  by.y = "Taxon",
  all.x = TRUE
)

tax_func <- tax_func[match(tax$ASV, tax_func$ASV), ]
rownames(tax_func) <- tax_func$ASV

for(g in groups){
  tax_func[[g]][is.na(tax_func[[g]])] <- 0
  tax_func[[g]] <- as.numeric(tax_func[[g]])
}

#### 4. OTU TABLE + METADATA ####

ASV <- as(otu_table(ps_16S), "matrix")

if(!taxa_are_rows(ps_16S)){
  ASV <- t(ASV)
}

storage.mode(ASV) <- "numeric"

rownames(ASV) <- taxa_names(ps_16S)
colnames(ASV) <- sample_names(ps_16S)

meta <- data.frame(sample_data(ps_16S))
meta$ID <- rownames(meta)

meta$Landuse <- factor(meta$Landuse, levels = landuse_levels)
meta$Type <- factor(meta$Type, levels = type_levels)
meta$Neighborhood <- factor(meta$Neighborhood)

#### 5. CHECK ASV COUNTS ####

asv_function_counts <- data.frame(
  Function = groups,
  ASV_count = sapply(groups, function(g){
    sum(tax_func[[g]] == 1, na.rm = TRUE)
  })
)

print(asv_function_counts)

#### 6. FUNCTION : ONE FUNCTIONAL GUILD ####

run_function_abundance <- function(function_name, type_name){

  cat("\n====================================\n")
  cat(type_name, "-", function_name, "\n")
  cat("====================================\n")

  guild_idx <- rownames(ASV)[tax_func[rownames(ASV), function_name] == 1]
  other_idx <- setdiff(rownames(ASV), guild_idx)

  if(length(guild_idx) < 2){
    message(paste("Skipped:", function_name, type_name, "- too few ASVs"))
    return(NULL)
  }

  guild_counts <- colSums(ASV[guild_idx, , drop = FALSE])
  other_counts <- colSums(ASV[other_idx, , drop = FALSE])

  guild_mat <- rbind(
    guild_counts,
    other_counts
  )

  rownames(guild_mat) <- c(function_name, "Other")
  colnames(guild_mat) <- colnames(ASV)
  storage.mode(guild_mat) <- "numeric"

  keep <- meta$ID[meta$Type == type_name]
  samples <- intersect(colnames(guild_mat), keep)

  mat <- guild_mat[, samples, drop = FALSE]
  meta_sub <- meta[match(samples, meta$ID), ]

  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  if(!function_name %in% rownames(mat)){
    message(paste("Skipped:", function_name, type_name, "- absent after filtering"))
    return(NULL)
  }

  if(sum(mat[function_name, ] > 0) < 5){
    message(paste("Skipped:", function_name, type_name, "- too few positive samples"))
    return(NULL)
  }

  if(nrow(mat) < 2 | ncol(mat) < 10){
    message(paste("Skipped:", function_name, type_name, "- matrix too small"))
    return(NULL)
  }

  meta_sub <- meta_sub[match(colnames(mat), meta_sub$ID), ]

  #### IMPORTANT : ALDEx2 compatible call ####

  ald <- aldex.clr(
    mat,
    as.character(meta_sub$Landuse),
    mc.samples = mc_n
  )

  mc <- ALDEx2::getMonteCarloInstances(ald)

  clr <- do.call(
    cbind,
    lapply(mc, function(x){
      apply(as.matrix(x), 1, median)
    })
  )

  colnames(clr) <- names(mc)

  df <- as.data.frame(clr)
  df$Function <- rownames(df)

  df <- tidyr::pivot_longer(
    df,
    cols = -Function,
    names_to = "ID",
    values_to = "CLR"
  )

  df <- merge(df, meta_sub, by = "ID")

  d <- df[df$Function == function_name, ]

  d$Landuse <- factor(d$Landuse, levels = landuse_levels)
  d$Neighborhood <- factor(d$Neighborhood)
  d$Type <- type_name
  d$Function <- function_name

  mod <- lmer(
    CLR ~ Landuse + (1 | Neighborhood),
    data = d
  )

  singular <- performance::check_singularity(mod)
  print(singular)

  aov_tab <- anova(mod)
  print(aov_tab)

  #### DHARMA ####

  sim_res <- DHARMa::simulateResiduals(
    fittedModel = mod,
    n = 1000
  )

  dharma_file <- paste0(
    path,
    "DHARMa_16S_",
    type_name,
    "_",
    function_name,
    ".png"
  )

  png(
    filename = dharma_file,
    width = 1800,
    height = 1400,
    res = 200
  )

  plot(sim_res)

  dev.off()

  uniformity_test <- DHARMa::testUniformity(sim_res)
  dispersion_test <- DHARMa::testDispersion(sim_res)
  outlier_test <- DHARMa::testOutliers(sim_res)

  print(uniformity_test)
  print(dispersion_test)
  print(outlier_test)

  diagnostics <- data.frame(
    Function = function_name,
    Type = type_name,
    Singular = singular,
    Uniformity_p = uniformity_test$p.value,
    Dispersion_p = dispersion_test$p.value,
    Outliers_p = outlier_test$p.value,
    DHARMa_plot = dharma_file
  )

  #### EMMEANS + LETTERS ####

  emm <- emmeans(mod, ~ Landuse)

  pairwise <- pairs(
    emm,
    adjust = "tukey"
  )

  print(pairwise)

  pairwise_df <- as.data.frame(pairwise)

  pvals <- pairwise_df$p.value
  names(pvals) <- gsub(" ", "", pairwise_df$contrast)

  letters_vec <- tryCatch(
    multcompView::multcompLetters(pvals)$Letters,
    error = function(e){
      setNames(rep("a", length(landuse_levels)), landuse_levels)
    }
  )

  letters_df <- data.frame(
    Landuse = landuse_levels,
    .group = letters_vec[landuse_levels],
    stringsAsFactors = FALSE
  )

  letters_df$.group[is.na(letters_df$.group)] <- "a"

  y_pos <- d %>%
    dplyr::group_by(Landuse) %>%
    dplyr::summarise(
      y = max(CLR, na.rm = TRUE),
      .groups = "drop"
    )

  letters_df$y <- y_pos$y[
    match(letters_df$Landuse, as.character(y_pos$Landuse))
  ]

  letters_df$y <- letters_df$y +
    0.12 * diff(range(d$CLR, na.rm = TRUE))

  letters_df$Function <- function_name
  letters_df$Type <- type_name
  letters_df$Landuse <- factor(letters_df$Landuse, levels = landuse_levels)

  result <- list(
    model = mod,
    anova = aov_tab,
    diagnostics = diagnostics,
    emmeans = emm,
    pairwise = pairwise,
    letters = letters_df,
    data = d
  )

  rm(
    ald, mc, clr, df, mat, guild_mat,
    guild_counts, other_counts, sim_res
  )
  gc()

  return(result)
}

#### 7. RUN ALL MODELS ####

results_all <- list()
data_list <- list()
letters_list <- list()
diagnostics_list <- list()

for(g in groups){

  for(t in type_levels){

    res <- tryCatch(
      run_function_abundance(g, t),
      error = function(e){
        message(paste("ERROR:", g, t, "-", e$message))
        return(NULL)
      }
    )

    if(!is.null(res)){

      key <- paste(g, t, sep = "_")

      results_all[[key]] <- res
      data_list[[key]] <- res$data
      letters_list[[key]] <- res$letters
      diagnostics_list[[key]] <- res$diagnostics
    }

    rm(res)
    gc()
  }
}

if(length(results_all) == 0){
  stop("No model was successfully fitted. Check ALDEx2 call, input matrix, or functional database.")
}

data_all <- dplyr::bind_rows(data_list)
letters_all <- dplyr::bind_rows(letters_list)
diagnostics_all <- dplyr::bind_rows(diagnostics_list)

data_all$Function <- factor(data_all$Function, levels = groups)
letters_all$Function <- factor(letters_all$Function, levels = groups)

data_all$Type <- factor(data_all$Type, levels = type_levels)
letters_all$Type <- factor(letters_all$Type, levels = type_levels)

#### 8. SUMMARY TABLE OF MODEL P-VALUES ####

model_summary <- data.frame()

for(name in names(results_all)){

  a <- results_all[[name]]$anova

  model_summary <- rbind(
    model_summary,
    data.frame(
      Model = name,
      Function = as.character(results_all[[name]]$data$Function[1]),
      Type = as.character(results_all[[name]]$data$Type[1]),
      F_value = a$`F value`[1],
      p_value = a$`Pr(>F)`[1],
      Singular = results_all[[name]]$diagnostics$Singular[1],
      Uniformity_p = results_all[[name]]$diagnostics$Uniformity_p[1],
      Dispersion_p = results_all[[name]]$diagnostics$Dispersion_p[1],
      Outliers_p = results_all[[name]]$diagnostics$Outliers_p[1],
      DHARMa_plot = results_all[[name]]$diagnostics$DHARMa_plot[1]
    )
  )
}

model_summary$Signif <- "ns"
model_summary$Signif[model_summary$p_value < 0.05] <- "*"
model_summary$Signif[model_summary$p_value < 0.01] <- "**"
model_summary$Signif[model_summary$p_value < 0.001] <- "***"

model_summary <- model_summary %>%
  dplyr::arrange(Type, Function)

print(model_summary)

#### 9. GLOBAL PLOT ####

p_all_functions_16S <- ggplot(
  data_all,
  aes(x = Landuse, y = CLR, fill = Landuse)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.85,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.35,
    size = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 1.8,
    color = "black"
  ) +
  geom_text(
    data = letters_all,
    aes(x = Landuse, y = y, label = .group),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  facet_grid(Type ~ Function, scales = "free_y") +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
    strip.text.y = element_text(face = "bold")
  ) +
  labs(
    title = "Putative bacterial functional guilds across the urbanization gradient",
    x = NULL,
    y = "CLR abundance"
  )

p_all_functions_16S

#### 10. SEPARATE PLOTS BY COMPARTMENT ####

plot_compartment <- function(type_name){

  ggplot(
    subset(data_all, Type == type_name),
    aes(x = Landuse, y = CLR, fill = Landuse)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.85,
      width = 0.65
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.35,
      size = 0.7,
      color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2,
      color = "black"
    ) +
    geom_text(
      data = subset(letters_all, Type == type_name),
      aes(x = Landuse, y = y, label = .group),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    ) +
    facet_wrap(
      ~ Function,
      scales = "free_y",
      ncol = 4
    ) +
    scale_fill_manual(values = cols) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste(type_name, "bacterial functional guilds"),
      x = NULL,
      y = "CLR abundance"
    )
}

p_leaf_functions_16S <- plot_compartment("Leaf")
p_root_functions_16S <- plot_compartment("Root")
p_soil_functions_16S <- plot_compartment("Soil")

p_leaf_functions_16S
p_root_functions_16S
p_soil_functions_16S
