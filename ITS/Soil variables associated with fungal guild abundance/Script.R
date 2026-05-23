library(phyloseq)
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)

#### 1. IMPORT ####

path <- "Data/"

ps <- readRDS(paste0(path, "phyloseq_no_negatives_ITS.rds"))
DATABASE <- read_excel(paste0(path, "Fungal_trait.xlsx"))
DATABASE <- as.data.frame(DATABASE)

data_env <- read.csv(
  paste0(path, "Data_env.csv"),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

#### 2. TAXONOMY ####

tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)
tax$Genus <- gsub("^g__", "", as.character(tax$Genus))

DATABASE$Genus <- gsub("^g__", "",
  DATABASE[[grep("GENUS|Genus", colnames(DATABASE), value=TRUE)[1]]]
)

merge_ft <- merge(tax, DATABASE, by="Genus", all.x=TRUE)

#### 3. OTU ####

otu <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)) otu <- t(otu)

meta <- data.frame(sample_data(ps))
meta$ID <- rownames(meta)

meta$Tree_ID <- sub("^[0-9]+_", "", meta$ID)
meta$Tree_ID <- sub("_(Soil|Root|Leaf)$", "", meta$Tree_ID)

#### 4. SOIL ONLY ####

keep <- meta$Type == "Soil"

otu <- otu[keep, , drop=FALSE]
meta <- meta[keep, , drop=FALSE]

#### 5. GUILDS ####

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
    grepl("saprotroph", merge_ft$primary_lifestyle, ignore.case=TRUE) |
    grepl("saprotroph", merge_ft$Secondary_lifestyle, ignore.case=TRUE)
  ]
)

#### 6. GUILD COUNTS ####

df_guild <- data.frame(ID = rownames(otu))

for(g in names(guild_list)){

  asv_keep <- intersect(colnames(otu), guild_list[[g]])

  if(length(asv_keep) == 0){
    df_guild[[g]] <- 0
  } else {
    df_guild[[g]] <- rowSums(otu[, asv_keep, drop=FALSE])
  }
}

#### 7. CLR ####

guild_mat <- as.matrix(df_guild[, -1])
guild_mat <- guild_mat + 1

clr_fun <- function(x){
  log(x / exp(mean(log(x))))
}

guild_clr <- t(apply(guild_mat, 1, clr_fun))
guild_clr <- as.data.frame(guild_clr)
colnames(guild_clr) <- colnames(guild_mat)

#### 8. ENV ####

env <- data_env[match(meta$Tree_ID, data_env$ID), ]

vars <- c(
  "Temperature","Moisture","pH","SOM",
  "NH4_ug.g","NO3_ug.g","OlsenP_ug.g"
)

env <- env[, vars, drop=FALSE]

#### 9. FINAL TABLE ####

df <- cbind(guild_clr, env)
df <- na.omit(df)

#### 10. CORRELATIONS ####

guilds <- c("ECM","Plant_pathogen","Animal_parasite","Saprotroph")

res <- data.frame()

for(g in guilds){
  for(v in vars){

    test <- cor.test(
      df[[g]],
      df[[v]],
      method = "spearman",
      exact = FALSE
    )

    res <- rbind(
      res,
      data.frame(
        Guild = g,
        Variable = v,
        Rho = unname(test$estimate),
        p.value = test$p.value
      )
    )
  }
}

print(res)

#### 11. LABELS ####

res$Signif <- ""
res$Signif[res$p.value < 0.05]  <- "*"
res$Signif[res$p.value < 0.01]  <- "**"
res$Signif[res$p.value < 0.001] <- "***"

#### 12. HEATMAP ####

ggplot(res, aes(x=Variable, y=Guild, fill=Rho)) +
  geom_tile(color="white") +
  geom_text(aes(label=Signif), size=6) +
  scale_fill_gradient2(
    low="steelblue",
    mid="white",
    high="firebrick",
    midpoint=0
  ) +
  theme_minimal(base_size=14) +
  labs(
    title="Spearman correlations between CLR guild abundance and soil variables",
    x=NULL,
    y=NULL,
    fill="Rho"
  ) +
  theme(
    axis.text.x = element_text(angle=35, hjust=1)
  )
