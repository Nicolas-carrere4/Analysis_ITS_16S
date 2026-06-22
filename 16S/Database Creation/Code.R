library(data.table)
library(tidyr)
library(stringr)
library(plyr)

#### 1. PATHS ####

path_db <- "Data/Database/"
path_out <- "Data/"

#### 2. IMPORT ####

fg <- read.csv(
  paste0(path_db, "bacteria_func_groups.csv"),
  stringsAsFactors = FALSE
)

N_cyclers_raw <- read.csv(
  paste0(path_db, "Npathways_Albright2018.csv"),
  stringsAsFactors = FALSE
)

cellulolytic_raw <- read.csv(
  paste0(path_db, "cellulolytic_Berlemont.csv"),
  stringsAsFactors = FALSE
)

#### 3. FORMAT N-CYCLE DATASET ####

N_cyclers <- N_cyclers_raw[
  ,
  !colnames(N_cyclers_raw) %in% c(
    "Samplename",
    "Genome",
    "StudyName",
    "Ecosystem",
    "Ecosystem.Category",
    "Ecosystem.Subtype",
    "Ecosystem.Type",
    "Environment",
    "Genome_size_assembled",
    "Gene_Count_assembled"
  )
]

setnames(
  N_cyclers,
  old = c(
    "Nitrogen.Fixation",
    "Assimilatory.Nitrite.to.ammonia",
    "Dissimilatory.Nitrite.to.Ammonia",
    "Assimilatory.Nitrate.to.Nitrite",
    "Dissimilatory.Nitrate.to.Nitrite"
  ),
  new = c(
    "N_fixation",
    "Assim_nitrite_reduction",
    "Dissim_nitrite_reduction",
    "Assim_nitrate_reduction",
    "Dissim_nitrate_reduction"
  ),
  skip_absent = TRUE
)

N_cyclers[N_cyclers == "complete"] <- 1
N_cyclers[N_cyclers == "incomplete"] <- 0
N_cyclers[N_cyclers == "None"] <- 0

func_cols_N <- setdiff(colnames(N_cyclers), "Genus")

for(col in func_cols_N){
  N_cyclers[[col]] <- as.numeric(N_cyclers[[col]])
  N_cyclers[[col]][is.na(N_cyclers[[col]])] <- 0
}

if("Partial_Nitrification" %in% colnames(N_cyclers)){
  N_cyclers$Nitrification[
    N_cyclers$Partial_Nitrification == 1
  ] <- 1
}

partial_denit <- c("Partial_NO", "Partial_N2O", "Partial_N2")
partial_denit <- partial_denit[partial_denit %in% colnames(N_cyclers)]

if(length(partial_denit) > 0){
  N_cyclers$Denitrification[
    rowSums(N_cyclers[, partial_denit, drop = FALSE]) > 0
  ] <- 1
}

N_cyclers <- N_cyclers[
  ,
  !colnames(N_cyclers) %in% c(
    "Partial_Nitrification",
    "Partial_NO",
    "Partial_N2O",
    "Partial_N2"
  )
]

N_cyclers$Taxonomic.level <- "Genus"
N_cyclers$Taxon <- N_cyclers$Genus
N_cyclers$Genus <- NULL

N_cyclers <- unique(N_cyclers)

#### 4. FORMAT CELLULOLYTIC DATASET ####

cellulolytic <- cellulolytic_raw[
  ,
  colnames(cellulolytic_raw) %in% c(
    "Strain",
    "GH5",
    "GH6",
    "GH8",
    "GH9",
    "GH12",
    "GH44",
    "GH45",
    "GH48"
  )
]

cellulolytic$genus <- word(cellulolytic$Strain, 1)

cellulolytic$genus[
  cellulolytic$genus == "Candidatus"
] <- word(
  cellulolytic$Strain[cellulolytic$genus == "Candidatus"],
  1,
  2
)

gh_cols <- setdiff(
  colnames(cellulolytic),
  c("Strain", "genus")
)

for(col in gh_cols){
  cellulolytic[[col]] <- as.numeric(cellulolytic[[col]])
  cellulolytic[[col]][is.na(cellulolytic[[col]])] <- 0
}

cellulolytic$Cellulolytic <- as.numeric(
  rowSums(cellulolytic[, gh_cols, drop = FALSE]) > 0
)

cellulolytic <- cellulolytic[
  cellulolytic$Cellulolytic == 1,
  c("genus", "Cellulolytic")
]

colnames(cellulolytic)[1] <- "Taxon"

cellulolytic$Taxonomic.level <- "Genus"

cellulolytic <- unique(
  cellulolytic[, c("Taxonomic.level", "Taxon", "Cellulolytic")]
)

#### 5. FORMAT LITERATURE FUNCTIONAL GROUPS ####

fg_lit <- fg[
  ,
  !colnames(fg) %in% c(
    "Classification.system",
    "Source",
    "Notes"
  )
]

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
  "Copiotroph",
  "Oligotroph"
)

for(group in groups){
  fg_lit[[group]] <- as.numeric(fg_lit$Classification == group)
}

fg_lit$Classification <- NULL

#### 6. COMBINE DATABASES ####

fg_out <- plyr::rbind.fill(
  cellulolytic,
  N_cyclers,
  fg_lit
)

fg_out[is.na(fg_out)] <- 0

fg_out$Taxon <- gsub(
  "^g__",
  "",
  as.character(fg_out$Taxon)
)

fg_out$Taxon <- trimws(fg_out$Taxon)

fg_out <- unique(fg_out)

#### 7. COLLAPSE DUPLICATE TAXA ####

func_cols <- setdiff(
  colnames(fg_out),
  c("Taxonomic.level", "Taxon")
)

fg_out <- aggregate(
  fg_out[, func_cols],
  by = list(
    Taxonomic.level = fg_out$Taxonomic.level,
    Taxon = fg_out$Taxon
  ),
  FUN = max
)

#### 8. CHECK OUTPUT ####

print(dim(fg_out))
print(colSums(fg_out[, func_cols, drop = FALSE]))
print(head(fg_out))

#### 9. SAVE ####

saveRDS(
  fg_out,
  paste0(path_out, "bacteria_tax_to_function.rds")
)

write.csv(
  fg_out,
  paste0(path_out, "bacteria_tax_to_function.csv"),
  row.names = FALSE
)

colSums(fg_out[, func_cols])
