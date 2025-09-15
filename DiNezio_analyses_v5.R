########################################
##
## Di Nezio et al.
## Microbial diversity of groundwater related environments in Lanzarote
## Ecological analyses
## Script written by Alejandro Martinez and Francesco Di Nezio
##
## v1. Verbania - 01/07/2024
## v2. Verbania - 12/07/2024
## v3. Verbania - 17/07/2024
## v4. Verbania - 30/07/2024
## v5. Verbania - 10/09/2025 (community composiytion plots, abundance, shannon)
## v6. Verbania - 13/09/2025 (LOOKING LIKE THE LAST... )
##
########################################


library(tidyr)
library(dplyr)
library(ggplot2)


###### PART 1 - ARRANGING DATA FILE -----------------------------------------------------


setwd("/Users/amartinez/Dropbox/_Papers/_DiNezio_16S_Lanzarote")


### cleaning community matrix

# comm <- read.csv2("Francesco_rawdata/sequence_tab.csv")
# taxonomy <- read.csv2("Francesco_rawdata/taxonomy.csv")
# names <- readxl::read_xlsx("names.xlsx")
# names$Sample <- gsub("sum_","",names$Sample, fixed = T)
# 
# comm.names <- as.data.frame(t(comm))
# comm.names$taxonomy <- rownames(comm.names)
# 
# colnames(comm.names) <- comm.names[1,]
# comm.names <- comm.names[-1,]
# 
# comm.names <- merge(comm.names, taxonomy[c("X","ID")], by = "X")
# row.names(comm.names) <- comm.names$ID
# comm.names <- comm.names[,c(-1,-32)]
# 
# comm.names[] <- lapply(comm.names, function(x) {
#   if (is.character(x)) gsub(" ", "", x) else x
# })
# 
# comm.names[] <- lapply(comm.names, function(x) as.numeric(as.character(x)))
# 
# names$Sample %in% colnames(comm.names)
# 
# name_vector <- setNames(names$Name, names$Sample)
# colnames(comm.names) <- ifelse(colnames(comm.names) %in% names(name_vector),
#                          name_vector[colnames(comm.names)],
#                          colnames(comm.names))
# 
# write.csv2(comm.names, "Francesco_rawdata/sequence_tab_asvs.csv")
# 
# 
# 
# stations <- read.csv2("Metadat.csv", sep = ";", dec = ".")

  
###### 1.1 - Arrange taxonomy-free community matrix ----------------------------------------

setwd("/Users/amartinez/Dropbox/_Papers/_DiNezio_16S_Lanzarote")

comm <- read.csv2("Francesco_rawdata/sequence_tab_asvs.csv", header = 1)
taxonomy <- read.csv2("Francesco_rawdata/taxonomy.csv")
stations <- read.csv2("Francesco_rawdata/Metadat.csv")
traits <- read.csv2("Francesco_rawdata/pathogenicity_review_v2.csv") ### This is the file you need to update



#### 1.1.1 - Clean community matrix -----------------------

colnames(comm) <- comm[1,]
comm <- comm[-1,]

row.names(comm) <- comm[,1]
comm <- comm[,-1]

comm[] <- lapply(comm, function(x) {
  if (is.character(x)) gsub(" ", "", x) else x
})

comm[] <- lapply(comm, function(x) as.numeric(as.character(x)))
str(comm)


#### 1.1.2 - Filtering out ASVs corresponding to wrong groups ----------

taxonomy[is.na(taxonomy)] <- "Unclassified"

archaea <- taxonomy[ which (taxonomy$Kingdom == "Archaea"),"ID"]
mitochondria <- taxonomy[ which (taxonomy$Family == "Mitochondria"),"ID"]
chloroplast <- taxonomy[ which (taxonomy$Order == "Chloroplast"),"ID"]
wrong.seq <- taxonomy[ which (taxonomy$Phylum == "Unclassified"),"ID"]

drop.seqs <- c(archaea,mitochondria,chloroplast,wrong.seq)

"%ni%" <- Negate("%in%")
comm <- comm[ , which (colnames(comm) %ni% drop.seqs)]

rm(archaea,mitochondria,chloroplast,wrong.seq, drop.seqs)




#### 1.1.3 - Deal with rare community matrix ----------

### In the taxonomy-free matrix, we filter out sequences present in less than 2 samples
### and totaling less  than 3 reads.

comm <- comm[,
  colSums(comm > 0) >= 2 & colSums(comm) >= 3
]

# comm.rare <- comm[,
#              colSums(comm > 0) < 2 & colSums(comm) < 3
# ]



###### 1.2 - Complete the stations matrix ----------------------------------------


#### 1.2.1 - Calculate richness per station ----------------------------


comm.rich <- comm
comm.rich[comm.rich > 0] <- 1

richness <- as.data.frame(rowSums(comm.rich))
colnames(richness) <- "richness"
richness$Sample <- rownames(richness)

stations <- merge(stations, richness, by = "Sample")
rownames(stations) <- stations$Sample

rm(comm.rich,richness)



#### 1.2.2 -  Calculate shannon diversity per station ----------

d.shannonH <- as.data.frame(vegan::diversity(comm, index = "shannon"))
colnames(d.shannonH) <- "shannon"
d.shannonH$Sample <- rownames(d.shannonH)

stations <- merge(stations, d.shannonH, by = "Sample")
rownames(stations) <- stations$Sample

rm(d.shannonH)



#### 1.2.3 - Add number of reads ----------

reads <- read.csv("number of reads.csv")
names <- as.data.frame(readxl::read_xlsx("names.xlsx"))


names$Sample <- gsub("sum_","", names$Sample, fixed = T)
reads <- merge(reads, names, by = "Sample")

stations <- merge(stations, reads[c("Name","Reads_PostChimera")], by.x = "Sample", by.y = "Name")
colnames(stations)[colnames(stations) == "Reads_PostChimera"] <- "reads"

rm(reads, names)


#### 1.2.4 -  Count buildings and measure m of road ----------


source("functions/GIS_covariates.R")

library(sf)
library(osmdata)
library(dplyr)
library(units)

# stations$buildings <- count_buildings_osm(latitudes = as.numeric(stations$Latitude),
#                                           longitudes = as.numeric(stations$Longitude),
#                                           radius = 1000)


stations$buildings <- c(69,69,15,475,0,0,68,4,
                        141,136,136,3,0,29,0,2,
                        0,2,2,2,0,2,0,0,344,0,2,
                        2,2,2)

# stations$roads <- calculate_road_lengths(latitudes = as.numeric(stations$Latitude),
#                                          longitudes = as.numeric(stations$Longitude),
#                                          radius = 1000)

stations$roads <-  c(8579.325,8579.325,5925.543,75528.494,5717.383,
                     5182.547,7890.656,6479.438,11556.677,10361.297,
                     10361.297,10195.618,10256.365,16542.026,6820.342,
                     6413.168,0,6857.950,0,9609.712,9609.712,9609.712,
                     4451.792,4400.386,24303.961,0,8980.558,7864.168,
                     6461.287,7449.642)



##### 1.3 - Complete taxonomy matrix ----------------------------------------

#### 1.3.1 -  Add taxon name -------------

taxonomy[is.na(taxonomy)] <- "Unclassified"

taxonomy <- taxonomy %>%
  mutate(
    Taxon = case_when(
      Genus != "Unclassified" & Species != "Unclassified" ~ paste(Genus, Species),
      Genus != "Unclassified" & Species == "Unclassified" ~ Genus,
      Genus == "Unclassified" & Species == "Unclassified" ~ coalesce(
        ifelse(Family != "Unclassified", Family, NA_character_),
        ifelse(Order != "Unclassified", Order, NA_character_),
        ifelse(Class != "Unclassified", Class, NA_character_),
        ifelse(Phylum != "Unclassified", Phylum, NA_character_)
      )
    )
  )




#### 1.3.2 -  Save clean taxonomy -------------

taxonomy.checked <- taxonomy[ which (taxonomy$ID %in% colnames(comm)),]

comm.names <- as.data.frame(t(comm))
comm.names$ID <- row.names(comm.names)
taxonomy.checked <- merge(taxonomy.checked,comm.names, by ="ID", all = F)

colnames(taxonomy.checked)[2] <- "Sequence"

write.table(taxonomy.checked, "taxonomy.checked.csv", dec = ".", sep = ";", row.names = FALSE)



#### 1.3.3 -  Save ASVs sequences as fasta -------------


source("functions/dataframe2fas.R")

dataframe2fas(taxonomy.checked[c("ID","Sequence")], "taxonomy.check.fasta")


rm(comm.names)



###### 1.3.4 - Delete samples with little reads ----------------------------------------

# Pozo de los Cocoteros (CSP) and Montaña de Arena (TDA_MA)

drop <- c("CSP","TDA_MA")

stations <- stations[ which (stations$Sample %ni% drop),]
comm <- comm[ which (rownames(comm) %ni% drop),]

rm(drop)





##### 1.4 - Definition of functional groups ----------------------------------------


pathogens <- traits[ which (traits$Pathogenicity == "2" |
                              traits$Pathogenicity == "3"),]
nrow(pathogens)


pathogens.asv <- taxonomy[ which(taxonomy$Taxon %in% pathogens$Taxon ), ]
comm.pathogens <- comm[ , which ( colnames(comm) %in%  pathogens.asv$ID)]
rich.pathogens <- as.data.frame(rowSums(as.data.frame(apply(comm.pathogens, 2, function(x) as.integer(x > 0)))))
rich.pathogens$Sample <- rownames(comm.pathogens)
colnames(rich.pathogens)[1] <- "richness.pathogens"

shannonH.patho <- as.data.frame(vegan::diversity(comm.pathogens, index = "shannon"))
shannonH.patho$Sample <- rownames(shannonH.patho)
colnames(shannonH.patho)[1] <- "shannon.pathogens"

abund.pathogens <- as.data.frame(rowSums(comm.pathogens))
abund.pathogens$Sample <- rownames(abund.pathogens)
colnames(abund.pathogens)[1] <- "abundance.pathogens"



waste <- traits[ which (traits$Origin == "W"),]

waste.asv <- taxonomy[ which(taxonomy$Taxon %in% waste$Taxon ), ]
comm.waste <- comm[ , which ( colnames(comm) %in%  waste.asv$ID)]
rich.waste <- as.data.frame(rowSums(as.data.frame(apply(comm.waste, 2, function(x) as.integer(x > 0)))))
rich.waste$Sample <- rownames(comm.waste)
colnames(rich.waste)[1] <- "richness.waste"

shannonH.waste <- as.data.frame(vegan::diversity(comm.waste, index = "shannon"))
shannonH.waste$Sample <- rownames(shannonH.waste)
colnames(shannonH.waste)[1] <- "shannon.waste"

abund.waste <- as.data.frame(rowSums(comm.waste))
abund.waste$Sample <- rownames(abund.waste)
colnames(abund.waste)[1] <- "abund.waste"

comm1 <- as.data.frame(t(comm))

environmental <- traits[ which (traits$Origin == "E"),]
nrow(environmental)

environmental.asv <- taxonomy[ which(taxonomy$Taxon %in% environmental$Taxon ), ]
comm.environmental <- comm[ , which ( colnames(comm) %in%  environmental.asv$ID)]
rich.environment <- as.data.frame(rowSums(as.data.frame(apply(comm.environmental, 2, function(x) as.integer(x > 0)))))
rich.environment$Sample <- rownames(comm.environmental)
colnames(rich.environment)[1] <- "richness.environment"

shannonH.environment <- as.data.frame(vegan::diversity(comm.environmental, index = "shannon"))
shannonH.environment$Sample <- rownames(shannonH.environment)
colnames(shannonH.environment)[1] <- "shannon.environment"

abund.environment <- as.data.frame(rowSums(comm.environmental))
abund.environment$Sample <- rownames(abund.environment)
colnames(abund.environment)[1] <- "abund.environment"


stations <- merge(stations, abund.pathogens, by = "Sample")
stations <- merge(stations, rich.pathogens, by = "Sample")
stations <- merge(stations, shannonH.patho, by = "Sample")
stations <- merge(stations, abund.waste, by = "Sample")
stations <- merge(stations, rich.waste, by = "Sample")
stations <- merge(stations, shannonH.waste, by = "Sample")
stations <- merge(stations, abund.environment, by = "Sample")
stations <- merge(stations, rich.environment, by = "Sample")
stations <- merge(stations, shannonH.environment, by = "Sample")


rm(pathogens, rich.pathogens, pathogens.asv,
   waste, waste.asv, rich.waste,
   rich.environment, environmental, environmental.asv,
   shannonH.environment,shannonH.patho,shannonH.waste,
   abund.environment,abund.pathogens,abund.waste)




######################################################################################################
###### GOAL 1: DESCRIPTIVE PART --------------------------------------- 
######################################################################################################

########## G 1.1. Functional groups overall -----------------------

# --- palette for overlap-aware categories (keep stable across figures) ---


cols_combo <- c(
  "Environmental only"          = "#2c7bb6",
  "Environmental + Pathogenic"  = "#5e3c99",
  "Anthropogenic only"          = "#fdae61",
  "Anthropogenic + Pathogenic"  = "#d95f02",
  "Pathogenic (origin unknown)" = "#d7191c",
  "Unassigned"                  = "grey80"
)

# --- map traits to ASVs; keep Pathogenic orthogonal to origin


traits_by_taxon <- traits %>%
  group_by(Taxon) %>%
  summarise(
    is_patho = any(Pathogenicity %in% c("2","3"), na.rm = TRUE),
    env_flag = any(Origin == "E", na.rm = TRUE),
    ant_flag = any(Origin == "W", na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    # keep Environmental vs Anthropogenic mutually exclusive;
    # ambiguous (both TRUE or both FALSE) -> origin = NA
    origin = dplyr::case_when(
      env_flag & !ant_flag ~ "Environmental",
      ant_flag & !env_flag ~ "Anthropogenic",
      TRUE                 ~ NA_character_
    )
  ) %>%
  select(Taxon, is_patho, origin, env_flag, ant_flag)

ambig <- traits_by_taxon %>% filter(env_flag & ant_flag)
if (nrow(ambig) > 0) message("Ambiguous origin for ", nrow(ambig), " taxa; origin set to NA.")


cat_map <- taxonomy.checked %>%
  select(ID, Taxon) %>%
  left_join(traits_by_taxon, by = "Taxon", relationship = "many-to-one")

dup_ids <- cat_map %>% count(ID) %>% filter(n > 1)
if (nrow(dup_ids) > 0) {
  stop("Join duplicated ", nrow(dup_ids), " ASV IDs. Investigate taxonomy or trait duplicates.")
}

cat_map <- cat_map %>%
  mutate(
    combo = case_when(
      origin == "Environmental"  &  is_patho ~ "Environmental + Pathogenic",
      origin == "Environmental"  & !is_patho ~ "Environmental only",
      origin == "Anthropogenic"  &  is_patho ~ "Anthropogenic + Pathogenic",
      origin == "Anthropogenic"  & !is_patho ~ "Anthropogenic only",
      is_patho & is.na(origin)                 ~ "Pathogenic (origin unknown)",
      TRUE                                     ~ "Unassigned"
    )
  )

source("functions/functional-plot.R")

comm_long <- as.data.frame(comm) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ID", values_to = "reads") %>%
  left_join(cat_map, by = "ID") %>%
  mutate(
    origin2 = fct_explicit_na(origin, na_level = "Unassigned"),
    combo   = fct_relevel(combo, names(cols_combo))
  )


region_names <- c("Environmental","Anthropogenic","Pathogenic")


asv_cnt <- counts_for_eulerr(asv_flags)
fit_asv  <- eulerr::euler(asv_cnt$counts)

p_euler_asv <- plot(
  fit_asv,
  fills      = list(fill = c("#2c7bb6","#fdae61","#d7191c"), alpha = 0.6),
  labels     = list(col  = "black"),
  quantities = list(type = "percent", format = list(digits = 1))
)
grid::grid.text(
  sprintf("Unassigned: %d ASVs (%.1f%% of all)",
          asv_cnt$unassigned,
          100 * asv_cnt$unassigned / (asv_cnt$unassigned + sum(asv_cnt$counts))),
  x = grid::unit(0.02, "npc"), y = grid::unit(0.98, "npc"),
  just = c("left","top"), gp = grid::gpar(cex = 0.8)
)

## Read-weighted
reads_by_id <- comm_long %>%
  group_by(ID) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(asv_flags, by = "ID")

read_cnt <- counts_for_eulerr(reads_by_id[, region_names], weights = reads_by_id$reads)
fit_reads <- eulerr::euler(read_cnt$counts)

p_euler_reads <- plot(
  fit_reads,
  fills      = list(fill = c("#2c7bb6","#fdae61","#d7191c"), alpha = 0.6),
  labels     = list(col  = "black"),
  quantities = list(type = "percent", format = list(digits = 1))
)
grid::grid.text(
  sprintf("Unassigned: %s reads (%.1f%% of all)",
          scales::comma(read_cnt$unassigned),
          100 * read_cnt$unassigned / (read_cnt$unassigned + sum(read_cnt$counts))),
  x = grid::unit(0.02, "npc"), y = grid::unit(0.98, "npc"),
  just = c("left","top"), gp = grid::gpar(cex = 0.8)
)




########## G 1.2. Functional groups per sample -----------------------


hab_order <- c("cave","pool","salt","sea","well","pond")

## 4a) Abundance (reads)

reads_sample_combo <- comm_long %>%
  group_by(Sample, combo) %>%
  summarise(n_reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(p = n_reads / sum(n_reads)) %>%
  ungroup() %>%
  left_join(stations[, c("Sample","Type")], by = "Sample") %>%
  mutate(
    Type  = factor(Type, levels = hab_order),
    combo = fct_relevel(combo, names(cols_combo))
  ) %>%
  group_by(Type) %>%
  arrange(Sample, .by_group = TRUE) %>%
  mutate(Sample_ord = factor(Sample, levels = unique(Sample))) %>%
  ungroup()

p_hab_combo <- ggplot(reads_sample_combo, aes(x = Sample_ord, y = p, fill = combo)) +
  geom_col(color = "white", width = 0.95) +
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cols_combo, name = "Origin × Pathogenicity") +
  labs(x = "Sample", y = "Relative abundance",
       title = "Per-sample composition by origin × pathogenicity (reads)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        panel.grid.major.x = element_blank())

##  Richness (presence, share of categories present per sample)

reads_sample_combo_presence <- comm_long %>%
  mutate(present = as.integer(reads > 0)) %>%
  group_by(Sample, combo) %>%
  summarise(n = sum(present), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  left_join(stations[, c("Sample","Type")], by = "Sample") %>%
  mutate(
    Type       = factor(Type, levels = hab_order),
    combo      = fct_relevel(combo, names(cols_combo)),
    Sample_ord = factor(Sample, levels = levels(reads_sample_combo$Sample_ord))  # reuse same order
  )

p_hab_rich <- ggplot(reads_sample_combo_presence, aes(x = Sample_ord, y = p, fill = combo)) +
  geom_col(color = "white", width = 0.95) +
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cols_combo, name = "Origin × Pathogenicity") +
  labs(x = "Sample", y = "Relative richness",
       title = "Per-sample composition by origin × pathogenicity (richness)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        panel.grid.major.x = element_blank())


gridExtra::grid.arrange(
  p_hab_combo,
  p_hab_rich,
  nrow = 1, ncol = 2
)




###### G1.3. Taxonomic Community composition graphs -----------------------------------------


source("functions/taxonomy-plot.R")

library(pals)  
library(grid)

comm_long <- prep_comm_long(comm, taxonomy.checked)
comm_long.patho <- prep_comm_long(comm.pathogens, taxonomy.checked)
comm_long.anthro <- prep_comm_long(comm.waste, taxonomy.checked)
comm_long.envi <- prep_comm_long(comm.environmental, taxonomy.checked)


sample_order <- c("CUL1","CUL2","JDA_pool1","JDA_pool2","JDA_pool3",
                  "TDA_entrada","TDA_sima","TDA_LE",
                  "CHL","CHP","CHS","CHZ","MBB","MBS",
                  "FUC",
                  "CBC1","CBC2","CBL","CHG","JDA_beach1","JDA_beach2",
                  "CSS1","CSS2",
                  "FWH","FWM","GAF","SHW","TAB")

groups = c("cave","cave","cave","cave","cave",
           "cave","cave","cave",
           "pool","pool","pool","pool","pool","pool",
           "pond",
           "sea","sea","sea","sea","sea","sea",
           "salt","salt",
           "well","well","well","well","well")

stations$Sample

#### G1.3. Family-level graphs --------------------------

pf_all <- plot_comp_faceted_compact(
  dat_long = comm_long,
  tax_rank = "Family",
  top_n = 20,
  sample_order = sample_order,
  groups = groups,
  set_title = "A. All bacteria — Family",
  use_rel_abund = TRUE,
  legend_rows = 5,  
  legend_cols = NULL  
)

# pf_patho <- plot_comp_faceted_compact(
#   comm_long.patho, tax_rank = "Family",
#   top_n = 20,
#   set_title = 'C. Pathogenic bacteria by Family',
#   sample_order = sample_order,
#   groups = groups,
#   use_rel_abund = TRUE,
#   base_size = 8, x_text_size = 5, strip_size = 8,
#   legend_cols = 10, legend_key_mm = 3, legend_text = 6, legend_title = 7,
#   panel_spacing_mm = 2, bar_width = 0.85,
#   shorten_labels = TRUE
# )
# 
# 
# pf_anthr <- plot_comp_faceted_compact(
#   comm_long.anthro, tax_rank = "Family",
#   top_n = 20,
#   set_title = 'E. Anthropogenic bacteria by Family',
#   sample_order = sample_order,
#   groups = groups,
#   use_rel_abund = TRUE,
#   base_size = 8, x_text_size = 5, strip_size = 8,
#   legend_cols = 10, legend_key_mm = 3, legend_text = 6, legend_title = 7,
#   panel_spacing_mm = 2, bar_width = 0.85,
#   shorten_labels = TRUE
# )
# 
# 
# pf_envi  <- plot_comp_faceted_compact(
#   comm_long.envi, tax_rank = "Family",
#   top_n = 20,
#   set_title = 'G. Environmental bacteria by Family',
#   sample_order = sample_order,
#   groups = groups,
#   use_rel_abund = TRUE,
#   base_size = 8, x_text_size = 5, strip_size = 8,
#   legend_cols = 10, legend_key_mm = 3, legend_text = 6, legend_title = 7,
#   panel_spacing_mm = 2, bar_width = 0.85,
#   shorten_labels = TRUE
# )
# 
# gridExtra::grid.arrange(
#   pf_all$plot,
#   pf_patho$plot,
#   pf_anthr$plot,
#   pf_envi$plot,
#   nrow = 2, ncol = 2
# )




#### G1.4. Genus -level graphs --------------------------

pg_all <-  plot_comp_faceted_compact(
  dat_long = comm_long,
  tax_rank = "Genus",
  top_n = 20,
  sample_order = sample_order,
  groups = groups,
  set_title = "B. All bacteria — Genus",
  use_rel_abund = TRUE,
  legend_rows = 5,        # << ensure 4 rows
  legend_cols = NULL      # (or set legend_cols instead, but use one of them)
)
  

pg_patho <- plot_comp_faceted_compact(
  comm_long.patho, tax_rank = "Genus",
  top_n = 20,
  set_title = 'A. Pathogenic, by Genus',
  sample_order = sample_order,
  groups = groups,
  use_rel_abund = TRUE,
  base_size = 8, x_text_size = 5, strip_size = 8,
  legend_rows = 5, legend_key_mm = 3, legend_text = 6, legend_title = 7,
  panel_spacing_mm = 2, bar_width = 0.85,
  shorten_labels = TRUE
)


pg_anthr <- plot_comp_faceted_compact(
  comm_long.anthro, tax_rank = "Genus",
  top_n = 20,
  set_title = 'B. Anthropogenic, Genus',
  sample_order = sample_order,
  groups = groups,
  use_rel_abund = TRUE,
  base_size = 8, x_text_size = 5, strip_size = 8,
  legend_rows = 5, legend_key_mm = 3, legend_text = 6, legend_title = 7,
  panel_spacing_mm = 2, bar_width = 0.85,
  shorten_labels = TRUE
)


pg_envi  <- plot_comp_faceted_compact(
  comm_long.envi, tax_rank = "Genus",
  top_n = 20,
  set_title = 'C. Environmental, by Genus',
  sample_order = sample_order,
  groups = groups,
  use_rel_abund = TRUE,
  base_size = 8, x_text_size = 5, strip_size = 8,
  legend_rows = 5, legend_key_mm = 3, legend_text = 6, legend_title = 7,
  panel_spacing_mm = 2, bar_width = 0.85,
  shorten_labels = TRUE
)


#### Figure 2A: all community


gridExtra::grid.arrange(
  pf_all$plot,
  pg_all$plot,
  nrow = 1, ncol = 2
)


gridExtra::grid.arrange(
  pg_patho$plot,
  pg_anthr$plot,
  pg_envi$plot,
  nrow = 1, ncol = 3
)



rm(pf_all,pf_anthr,pf_envi,pf_patho,
   pg_all,pg_anthr,pg_envi,pg_patho,
   comm_long,comm_long.anthro,comm_long.envi,comm_long.patho)





#####################################################################################################
########## GOAL 2 - INFERENTIAL PART  -----------------------------------------------------
######################################################################################################


###### G 2.1 - Total bacterial communities in Lanzarote ----------------------------------------



###### G 2.1.1 - ASV richness global dataset ----------------------------------------


psych::pairs.panels(stations[c(4:12,14:17)]) ### Buildings and roads (0.73); Altitude and distance sea (0.93) correlated



mod.rich <- MASS::glm.nb(richness ~ Type + reads + Distance.sea + log(buildings + 1),
                       data = stations)


performance::check_model(mod.rich) 
performance::check_overdispersion(mod.rich) 
car::Anova(mod.rich)
summary(mod.rich)

emm_type <- emmeans::emmeans(mod.rich, ~ Type)
tukey_type <- pairs(emm_type, adjust = "tukey")
summary(tukey_type)





###### G 2.1.2 - ASV richness visualization ----------------------------------------


# --- Custom palette ---
my_palette <- c(
  cave = "#2D2D2D",
  sea = "#1B4F72",
  pond = "#1E8378",
  well = "#D9CBA3",
  pool = "#6E9C6D",
  salt = "#C6367B" 
)

# --- Plot ---
(p_rich.all <- ggplot(stations, aes(x = Type, y = richness, color = Type)) +
                geom_point(size = 2, position = position_jitter(width = 0, height = 0)) + # vertically aligned
                stat_summary(fun = mean, geom = "crossbar", linetype = "dashed",
                             width = 0.6, fatten = 2, color = "red") +
                scale_color_manual(values = my_palette) +
                labs(x = "", y = "Richness",
                     title = "") +
                theme_minimal() +
                theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                      axis.text.y  = element_text(size = 6),
                      axis.title.y = element_text(size = 6),
                      legend.position = "none"))





###### G 2.1.2 - ASV community composition analyses global dataset ----------------------------------------
  

beta.tax <- BAT::beta(comm, abund = T)

mean(beta.tax$Btotal)
mean(beta.tax$Brepl / beta.tax$Btotal)
mean(beta.tax$Brich / beta.tax$Btotal)

source("functions/beta_plot2.R")



(pb.all <- beta_density_plot(comm, stations$Type, beta.tax,
                              fill_values = c("Within habitat" = "#1b9e77",
                                              "Across habitats" = "#d95f02"),
                              include_total = TRUE,
                              include_nestedness = FALSE,
                              facet_cols = 2,
                              alpha = 0.6))





## Total beta diversity
model.Btotal <- beta.tax$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads

(perm.taxT <- vegan::adonis2(model.Btotal, stations, permutations = 9999, by = "terms"))


## Beta nestness
model.Bnest <- beta.tax$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi <- vegan::adonis2(model.Bnest, stations, permutations = 9999, by = "terms"))



## Beta turnover
model.Bturn <- beta.tax$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe <- vegan::adonis2(model.Bturn, stations, permutations = 9999, by = "terms"))

rm(model.Btotal,model.Bnest,model.Bturn,perm.taxRi,perm.taxRe,perm.taxT)




###### G 2.1.3 - ASV nMDS visualization ----------------------------------------

source("functions/plot_nmds.R")

library(vegan)
library(ggplot2)


my_palette <- c(
  cave = "#2D2D2D",
  sea = "#1B4F72",
  pond = "#1E8378",
  well = "#D9CBA3",
  pool = "#6E9C6D",
  salt = "#C6367B" 
)


# Prepare stations once
rownames(stations) <- stations$Sample
stations$Type <- factor(stations$Type, levels = names(my_palette))


# Plot all
(nmds1 <- plot_nmds(comm, 
                    stations, 
                    my_palette,
                    "NMDS",
                    point_size = 2.5,
                    label_size = 0,
                    label_vjust = -1))




######### G 2.1.4 - OVERALL FIGURE -------------------------------


library(patchwork)
library(ggplot2)

# Layout: row1 = A | B ; row2 = C spans both columns
design <- "
AB
CC
"

fig_small <-
  ((p_rich.all + guides(fill = "none", color = "none")) +  # hide Δ-richness legend
     nmds1 +
     pb.all) +
  plot_layout(design = design, guides = "collect", heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.title   = element_text(size = 10, face = "bold", hjust = 0,
                                margin = ggplot2::margin(b = 4)),
    strip.text   = element_text(size = 8),
    axis.title   = element_text(size = 6),
    axis.text    = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5),
    plot.subtitle= element_text(size = 6),
    plot.caption = element_text(size = 6),
    plot.tag     = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1)
  )
fig_small





######################################################################################################
###### G 2.2 - Functional bacterial groups in Lanzarote ----------------------------------------


library(glmmTMB)

###### G 2.2.1a - Pathogenes abundance and richness  -------------------------------------------------------

stations$Type <- relevel(factor(stations$Type), ref = "sea")

mod.patho.abund.prob <- glmmTMB(cbind(abundance.pathogens, reads - abundance.pathogens) ~ 
                              Type + Distance.sea + log(buildings + 1),
                      family=betabinomial(link = "logit"), data=stations)


performance::check_model(mod.patho.abund.prob) 
performance::check_overdispersion(mod.patho.abund.prob) 
car::Anova(mod.patho.abund.prob)
summary(mod.patho.abund.prob)

emm_type.patho.abund.prob <- emmeans::emmeans(mod.patho.abund.prob, ~ Type)
tukey_type.patho.abund.prob <- pairs(emm_type.patho.abund.prob, adjust = "tukey")
summary(tukey_type.patho.abund.prob)


rm(emm_type.patho.abund.prob, tukey_type.patho.abund.prob)



mod.patho.prob <- glm(cbind(richness.pathogens, richness - richness.pathogens) ~ 
                        Type + reads + log(buildings+1),
                      family=binomial(logit), data=stations)



performance::check_model(mod.patho.prob) 
performance::check_overdispersion(mod.patho.prob) 
car::Anova(mod.patho.prob)
summary(mod.patho.prob)

emm_type.patho.prob <- emmeans::emmeans(mod.patho.prob, ~ Type)
tukey_type.patho.prob <- pairs(emm_type.patho.prob, adjust = "tukey")
summary(tukey_type.patho.prob)


rm(emm_type.patho.prob, tukey_type.patho.prob)



###### G 2.2.1b - Anthropogenic abundance and richness  ---------------------------------------------------


mod.waste.abund.prob <- glmmTMB(cbind(abund.waste, reads - abund.waste) ~ 
                                  Type + Distance.sea + log(buildings + 1),
                                family=betabinomial(link = "logit"), data=stations)


performance::check_model(mod.waste.abund.prob) 
performance::check_overdispersion(mod.waste.abund.prob) 
car::Anova(mod.waste.abund.prob)
summary(mod.waste.abund.prob)

emm_type.waste.abund.prob <- emmeans::emmeans(mod.waste.abund.prob, ~ Type)
tukey_type.waste.abund.prob <- pairs(emm_type.waste.abund.prob, adjust = "tukey")
summary(tukey_type.waste.abund.prob)


rm(emm_type.waste.abund.prob, tukey_type.waste.abund.prob)




mod.waste.prob <- glm(cbind(richness.waste, richness - richness.waste) ~ 
                        Type + reads + Distance.sea + log(buildings + 1),
                      family=binomial(logit), data=stations)

summary(mod.waste.prob)
car::Anova(mod.waste.prob)

emm_type.waste.prob <- emmeans::emmeans(mod.waste.prob, ~ Type)
tukey_type.waste.prob <- pairs(emm_type.waste.prob, adjust = "tukey")
summary(tukey_type.waste.prob)


rm(emm_type.waste.prob, emm_type.waste.prob, tukey_type.waste.prob)


###### G 2.2.1c - Environment abundance and richness ------------------------------------


mod.envi.abund.prob <- glmmTMB(cbind(abund.environment, reads - abund.environment) ~ 
                                  Type + Distance.sea + log(buildings + 1),
                                family=betabinomial(link = "logit"), data=stations)


performance::check_model(mod.envi.abund.prob) 
performance::check_overdispersion(mod.envi.abund.prob) 
car::Anova(mod.envi.abund.prob)
summary(mod.envi.abund.prob)

emm_type.envi.abund.prob <- emmeans::emmeans(mod.envi.abund.prob, ~ Type)
tukey_type.envi.abund.prob <- pairs(emm_type.envi.abund.prob, adjust = "tukey")
summary(tukey_type.envi.abund.prob)



mod.env.prob <- glm(cbind(richness.environment, richness - richness.environment) ~ 
                      Type + reads +  Distance.sea + log(buildings+1),
                    family=binomial(logit), data=stations)
summary(mod.env.prob)
car::Anova(mod.env.prob)

emm_type.env.prob <- emmeans::emmeans(mod.env.prob, ~ Type)
tukey_type.env.prob <- pairs(emm_type.env.prob, adjust = "tukey")
summary(tukey_type.env.prob)



###### G 2.2.1d - Heatmap difference richness and abundace per habitat ----------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
source("functions/heatmap.R")

h.rich_patho <- mk_heat(mod.patho.prob,
                        title = "Pathogens richness",
                        ord   = c("cave","pond","pool","salt","sea","well"),
                        label_size = 2,
                        text_white_threshold = 0.045,
                        limit = 0.1,
                        tile_ratio = NULL)

h.rich_waste <- mk_heat(mod.waste.prob, 
                        title = "Anthropogenic bacteria richness",
                        ord = c("cave","pond","pool","salt","sea","well"),
                        label_size = 2,
                        text_white_threshold = 0.045,
                        limit = 0.1,
                        tile_ratio = NULL)

h.rich_envi   <- mk_heat(mod.env.prob, 
                         title = "Environmental bacteria richness",
                          ord = c("cave","pond","pool","salt","sea","well"),
                         label_size = 2,
                         text_white_threshold = 0.045,
                         limit = 0.1,
                         tile_ratio = NULL)


h.abund_patho <- mk_heat(mod.patho.abund.prob, 
                         title = "Pathogenic bacteria abundance",
                         ord = c("cave","pond","pool","salt","sea","well"),
                         label_size = 2,
                         text_white_threshold = 0.2,
                         limit = 0.327,
                         tile_ratio = NULL)

h.abund_waste <- mk_heat(mod.waste.abund.prob, 
                         title = "Anthropogenic bacteria abundance",
                         ord = c("cave","pond","pool","salt","sea","well"),
                         label_size = 2,
                         text_white_threshold = 0.2,
                         limit = 0.327,
                         tile_ratio = NULL)

h.abund.envi   <- mk_heat(mod.envi.abund.prob, 
                          title = "Environmental bacteria abundance",
                          ord = c("cave","pond","pool","salt","sea","well"),
                          label_size = 2,
                          text_white_threshold = 0.2,
                          limit = 0.327,
                          tile_ratio = NULL)





gridExtra::grid.arrange(h.abund_patho,h.abund_waste,h.abund.envi,
                        h.rich_patho,h.rich_waste,h.rich_envi, nrow=2,ncol=3)




rm(mod.environment,mod.env.prob, emm_type.env, tukey_type.env, mod.envi.probl.wr)

rm(pw1,pw2,pw3)




###### G 2.2.2a - Beta diversity pathogens  ---------------------------------

source("functions/beta_plot2.R")

beta.patho <- BAT::beta(comm.pathogens, abund=T) 
print(beta.patho)

mean(beta.patho$Btotal)
mean(beta.patho$Brepl / beta.patho$Btotal)
mean(beta.patho$Brich / beta.patho$Btotal)


(pb.patho <- beta_density_plot(comm.pathogens, stations$Type, beta.patho,
                             fill_values = c("Within habitat" = "#1b9e77",
                                             "Across habitats" = "#d95f02"),
                             include_total = TRUE,
                             include_nestedness = FALSE,
                             title = "Pathogenic β-diversity",
                             facet_cols = 1,
                             alpha = 0.6))



## Total beta diversity
model.Btotal.patho <- beta.patho$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads

(perm.taxT.patho <- vegan::adonis2(model.Btotal.patho, stations, permutations = 999, by = "terms"))

## Beta nestness
model.Bnest.patho <- beta.patho$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi.patho <- vegan::adonis2(model.Bnest.patho, stations, permutations = 999, by = "terms"))


## Beta turnover
model.Bturn.patho <- beta.patho$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe.patho <- vegan::adonis2(model.Bturn.patho, stations, permutations = 999, by = "terms"))

rm(model.Btotal.patho,model.Bnest.patho,model.Bturn.patho,perm.taxRi.patho,perm.taxRe.patho,perm.taxT.patho)


###### G 2.2.2b - Beta diversity Anthropogenic  ---------

beta.waste <- BAT::beta(comm.waste, abund=T) 

mean(beta.waste$Btotal)
mean(beta.waste$Brepl / beta.waste$Btotal)
mean(beta.waste$Brich / beta.waste$Btotal)


(pb.waste <- beta_density_plot(comm.waste, stations$Type, beta.waste,
                             fill_values = c("Within habitat" = "#1b9e77",
                                             "Across habitats" = "#d95f02"),
                             include_total = TRUE,
                             include_nestedness = FALSE,
                             title = "Anthropogenic β-diversity",
                             facet_cols = 1,
                             alpha = 0.6))



## Total beta diversity
model.Btotal.waste <- beta.waste$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads

(perm.taxT.waste <- vegan::adonis2(model.Btotal.waste, stations, permutations = 9999, by = "terms"))

## Beta nestness
model.Bnest.waste <- beta.waste$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi.waste <- vegan::adonis2(model.Bnest.waste, stations, permutations = 9999, by = "terms"))


## Beta turnover
model.Bturn.waste <- beta.waste$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe.waste <- vegan::adonis2(model.Bturn.waste, stations, permutations = 9999, by = "terms"))

rm(model.Btotal.waste,model.Bnest.waste,model.Bturn.waste,perm.taxRi.waste,perm.taxRe.waste,perm.taxT.waste)



###### G 2.2.2c - Beta diversity Environmental  ---------

beta.envi <- BAT::beta(comm.environmental, abund=T)

mean(beta.envi$Btotal)
mean(beta.envi$Brepl / beta.envi$Btotal)
mean(beta.envi$Brich / beta.envi$Btotal)

(pb.envi <- beta_density_plot(comm.environmental, stations$Type, beta.envi,
                               fill_values = c("Within habitat" = "#1b9e77",
                                               "Across habitats" = "#d95f02"),
                               include_total = TRUE,
                               include_nestedness = FALSE,
                              title = "Environmental β-diversity",
                               facet_cols = 1,
                               alpha = 0.6))



## Environmental
model.Btotal.envi <- beta.envi$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads

(perm.taxT.envi <- vegan::adonis2(model.Btotal.envi, stations, permutations = 9999, by = "terms"))

## Beta nestness
model.Bnest.envi <- beta.envi$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi.envi <- vegan::adonis2(model.Bnest.envi, stations, permutations = 9999, by = "terms"))


## Beta turnover
model.Bturn.envi <- beta.envi$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe.envi <- vegan::adonis2(model.Bturn.envi, stations, permutations = 9999, by = "terms"))


gridExtra::grid.arrange(pb.patho,pb.waste,pb.envi, nrow=1,ncol=3)
rm(model.Btotal.envi,model.Bnest.envi,model.Bturn.envi,perm.taxRi.envi,perm.taxRe.envi,perm.taxT.envi)



###### G 2.2.2d - ASV nMDS visualization  -----------------------------------------------------


source("functions/plot_nmds.R")

library(vegan)
library(ggplot2)


my_palette <- c(
  cave = "#2D2D2D",
  sea = "#1B4F72",
  pond = "#1E8378",
  well = "#D9CBA3",
  pool = "#6E9C6D",
  salt = "#C6367B" 
)


# Prepare stations once
rownames(stations) <- stations$Sample
stations$Type <- factor(stations$Type, levels = names(my_palette))

# Plot all
nmds.p <- plot_nmds(comm.pathogens, 
                    stations, 
                    my_palette,
                    "NMDS pathogenic bacteria",
                    point_size = 1,
                    label_size = 0)

nmds.w <- plot_nmds(comm.waste,
                    stations,
                    my_palette,
                    "NMDS anthropogenic bacteria",
                    point_size = 1,
                    label_size = 0)

nmds.e <- plot_nmds(comm.environmental,
                    stations,
                    my_palette,
                    "NMDS environmental bacteria",
                    point_size = 1,
                    label_size = 0)

gridExtra::grid.arrange(nmds.p,
                        nmds.w,
                        nmds.e,
                        ncol = 3,
                        nrow=1)



########## G 2.3 - FIGURE FUNCTIONAL  -------------------------------------

final <-
  (h.abund_patho | h.abund_waste | h.abund.envi) /
  (h.rich_patho  | h.rich_waste  | h.rich_envi)  /
  (pb.patho      | pb.waste      | pb.envi)      /
  (nmds.p        | nmds.w        | nmds.e) +
  plot_layout(guides = "collect", heights = c(1, 1, 1.5, 1.5)) &
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10, face = "bold", hjust = 0,
                              margin = ggplot2::margin(b = 4)),
    strip.text = element_text(size = 8),
    axis.title   = element_text(size = 6),
    axis.text    = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5),
    plot.subtitle= element_text(size = 6),
    plot.caption = element_text(size = 6),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1)
  )  

final



#####################################################################################################
########## G 2.3 - Classification  -----------------------------------------------------


########## G 2.3.0 - Generate a taxonomically-ranked community matrix -------------------------------------------


comm.tax <- as.data.frame(t(comm))
comm.tax$asvs <- row.names(comm.tax)

comm.tax <- merge(taxonomy, comm.tax, by.x = "ID", by.y = "asvs", all.y = T)

comm.tax <- comm.tax %>%
  mutate(
    Taxon = case_when(
      Genus != "Unclassified" & Species != "Unclassified" ~ paste(Genus, Species),
      Genus != "Unclassified" & Species == "Unclassified" ~ Genus,
      Genus == "Unclassified" & Species == "Unclassified" ~ coalesce(
        ifelse(Family != "Unclassified", Family, NA_character_),
        ifelse(Order != "Unclassified", Order, NA_character_),
        ifelse(Class != "Unclassified", Class, NA_character_),
        ifelse(Phylum != "Unclassified", Phylum, NA_character_)
      )
    )
  )






########## G 2.3.1 - RF Entire community-------------------------------------------

source("functions/summarytaxonomy.R")

comm.tax.list <- summarize_by_taxonomic_rank(comm.tax)


library(randomForest)
library(caret)
library(compositions)

source("functions/randomforest_onebyone.R")


rf.all <- run_rf_one_vs_rest_by_rank(
  comm.tax.list = comm.tax.list,
  stations      = stations,
  target_col    = "Type",
  ranks         = c("Family","Genus","Species"),
  ntree         = 10000,
  balance       = TRUE,
  plot          = TRUE
)



########## G 2.3.2 - RF Pathogens community-------------------------------------------


source("functions/summarytaxonomy.R")


comm.pathogens.tax <- as.data.frame(t(comm.pathogens))
comm.pathogens.tax$ID <- row.names(comm.pathogens.tax)
comm.pathogens.tax <- merge(taxonomy, comm.pathogens.tax, by = "ID")

comm.patho.list <- summarize_by_taxonomic_rank(comm.pathogens.tax,
                                               c("Family", "Genus", "Species"),
                                               min_reads = 20)




rf.pathogens <- run_rf_one_vs_rest_by_rank(
  comm.tax.list = comm.patho.list,
  stations      = stations,
  target_col    = "Type",
  ranks         = c("Family","Genus","Species"),
  ntree         = 10000,
  balance       = TRUE,
  plot          = TRUE
)




########## G 2.3.3 - RF Anthropogenic community -------------------------------------------


source("functions/summarytaxonomy.R")


comm.waste.tax <- as.data.frame(t(comm.waste))
comm.waste.tax$ID <- row.names(comm.waste.tax)
comm.waste.tax <- merge(taxonomy, comm.waste.tax, by = "ID")

comm.waste.list <- summarize_by_taxonomic_rank(comm.waste.tax,
                                               c("Family", "Genus", "Species"),
                                               min_reads = 20)




rf.waste <- run_rf_one_vs_rest_by_rank(
                comm.tax.list = comm.waste.list,
                stations      = stations,
                target_col    = "Type",
                ranks         = c("Family","Genus","Species"),
                ntree         = 10000,
                balance       = TRUE,
                plot          = TRUE
              )
               




########## G 2.3.4 -  RF Environmental -------------------------------------------

source("functions/summarytaxonomy.R")


comm.environmental.tax <- as.data.frame(t(comm.environmental))
comm.environmental.tax$ID <- row.names(comm.environmental.tax)
comm.environmental.tax <- merge(taxonomy, comm.environmental.tax, by = "ID")

comm.environmental.list <- summarize_by_taxonomic_rank(comm.environmental.tax,
                                                       c("Family", "Genus", "Species"),
                                                       min_reads = 20)



rf.environmental <- run_rf_one_vs_rest_by_rank(
                        comm.tax.list = comm.environmental.list,
                        stations      = stations,
                        target_col    = "Type",
                        ranks         = c("Family","Genus","Species"),
                        ntree         = 10000,
                        balance       = TRUE,
                        plot          = TRUE
                      )
                        
  
  
  
########## G 2.3.5 -  RF Visualization -------------------------------------------


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

source("functions/randomforest_onebyone.R")

col_low  <- "#d7191c" 
col_mid  <- "white"
col_high <- "#2c7bb6" 

hab_order   <- c("cave","pool","salt","sea","well")   
rank_order  <- c("Family","Genus","Species")
group_order <- c("All","Pathogens","Anthropogenic","Environmental")

# Build tidy metrics table
m_all <- bind_rows(
  bind_metrics(rf.all,           "All"),
  bind_metrics(rf.pathogens,     "Pathogens"),
  bind_metrics(rf.waste,         "Anthropogenic"),
  bind_metrics(rf.environmental, "Environmental")
) %>%
  # keep known habitats and drop pond in one go
  filter(habitat %in% c(hab_order, "pond")) %>%
  filter(habitat != "pond") %>%
  mutate(
    group    = factor(group,   levels = group_order),
    habitat  = factor(habitat, levels = hab_order),
    rank     = factor(rank,    levels = rank_order),
    Accuracy = 1 - OOB
  ) %>%
  # drop rows where rank not in rank_order (e.g., Order)
  filter(!is.na(rank)) %>%
  arrange(group, habitat, rank)


# AUC heatmap
(p_auc <- ggplot(m_all, aes(habitat, rank, fill = AUC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(AUC), "NA", sprintf("%.2f", AUC))), size = 3) +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0.85, limits = c(0.6, 1),
                       breaks = c(0, 0.85, 1), labels = c("0", "0.85", "1"),
                       oob = scales::squish, name = "AUC") +
  labs(title = "A. ROC AUC", x = NULL, y = NULL) +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank()))


# Accuracy heatmap
(p_acc <- ggplot(m_all, aes(habitat, rank, fill = Accuracy)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Accuracy)), size = 3) +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0.80, limits = c(0.6, 1),
                       breaks = c(0, 0.80, 1), labels = c("0", "0.80", "1"),
                       oob = scales::squish, name = "Accuracy") +
  labs(title = "B. Accuracy (1 − OOB)", x = NULL, y = NULL) +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank()))



final.rf <-
  (p_auc / p_acc) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 9),
    axis.text       = element_text(size = 9),
    strip.text      = element_text(size = 9)  # facet labels
  )
final.rf






 
#####################################################################################################
##### SUPPLEMENTARY ANALYSES --------------------------------------- 
#####################################################################################################
 
###### S 1 - Dendrograms  -------------------------------------------

# Set up a 2x2 plotting layout
par(mfrow = c(2, 2))

# Plot the four dendrograms
plot(hclust(beta.tax$Btotal, method = "average"), main = "All Taxa", xlab = "", sub = "")
plot(hclust(beta.patho$Btotal, method = "average"), main = "Pathogenic Taxa", xlab = "", sub = "")
plot(hclust(beta.waste$Btotal, method = "average"), main = "Wastewater Taxa", xlab = "", sub = "")
plot(hclust(beta.envi$Btotal, method = "average"), main = "Environmental Taxa", xlab = "", sub = "")

dev.off()



rm(comm.environmental,comm.environmental.tax,comm.environmental.list,
   comm.pathogens,comm.pathogens.tax,comm.patho.list,
   comm.waste,comm.waste.tax,comm.waste.list)


rm(pb.envi,pb.patho,pb.phylo,pb.waste,pb.all,pb.cave)






########## S 2.1 - Richness within the La Corona  -----------------------------------------------------
 
 
 ## Test for a gradient of antropogenic pollution in the cave
 
 stations.cave <- stations[ which (stations$Sample %in% c("TDA_sima",
                                                          "TDA_entrada",
                                                          "TDA_LE",
                                                          "JDA_pool1",
                                                          "JDA_pool2",
                                                          "JDA_pool3",
                                                          "CUL1",
                                                          "CUL2",
                                                          "CHL",
                                                          "CHZ") ),]
 
 
 # stations.cave$Distance.center <- c(5250,1800,250,50,0,0,0,0,50,250)
 
 
 
 mod.rich.cave <- MASS::glm.nb(richness ~ Light + reads +  Distance.sea,
                               data = stations.cave)
 
 
 performance::check_model(mod.rich.cave) 
 performance::check_overdispersion(mod.rich.cave) 
 car::Anova(mod.rich.cave)
 summary(mod.rich.cave)
 
 
 emm_light <- emmeans::emmeans(mod.rich.cave, ~ Light)
 tukey_light <- pairs(emm_light, adjust = "tukey")
 summary(tukey_light)
 
 boxplot(richness ~ Light, data = stations.cave)
 
 plot(richness ~ Distance.sea, data = stations.cave)
 
 
 rm(mod.rich.cave, emm_light, tukey_light)
 
 
########## S 2.2 - Composition within the La Corona  -----------------------------------------------------
 

 comm.cave <- comm[ which (rownames(comm) %in% stations.cave$Sample),]
 comm.cave <- comm.cave[ , colSums(comm.cave) > 0 ]
 
 beta.tax.cave <- BAT::beta(comm.cave, abund=T) 
 
 mean(beta.tax.cave$Btotal)
 mean(beta.tax.cave$Brepl / beta.tax.cave$Btotal)
 mean(beta.tax.cave$Brich / beta.tax.cave$Btotal)
 
 (pb.cave <- beta_density_plot(comm.cave, stations.cave$Type, beta.tax.cave,
                               fill_values = c("Within habitat" = "#1b9e77",
                                               "Across habitats" = "#d95f02"),
                               alpha = 0.6))
 
 
 ## Total beta diversity
 model.Btotal.cave <- beta.tax.cave$Btotal ~ Light + reads +  Distance.sea
 (perm.taxT.cave <- vegan::adonis2(model.Btotal.cave, stations.cave, permutations = 9999, by = "terms"))
 
 ## Beta nestness
 model.Bnest.cave <- beta.tax.cave$Brich ~  Light + reads +  Distance.sea
 (perm.taxRi.cave <- vegan::adonis2(model.Bnest.cave, stations.cave, permutations = 9999, by = "terms"))
 
 
 ## Beta turnover
 model.Bturn.cave <- beta.tax.cave$Brepl ~ Light + reads +  Distance.sea
 (perm.taxRe.cave <- vegan::adonis2(model.Bturn.cave, stations.cave, permutations = 9999, by = "terms"))
 
 rm(model.Btotal.cave,model.Bnest.cave,model.Bturn.cave,perm.taxRi.cave,perm.taxRe.cave,perm.taxT.cave)
 
 rm(comm.cave,stations.cave)
 
 
 
 
########## S 3 - Phylogenetic diversity  -----------------------------------------------------
 

 ########## S 3.1 - Calculate phylogenetic tree  -----------------------------------------------------

  
 # ## Read alignment
 # alignment.all <- ape::read.FASTA("phylogenetic_diversity/taxonomy.check_Mafft.fasta")
 # 
 # ## Calculate distance matrix
 # dist <- ape::dist.dna(alignment.all, pairwise.deletion = TRUE) 
 # 
 # ## Compute tree
 # tree <- ape::njs(dist)
 # tree <- phangorn::midpoint(tree)
 # ape::write.tree(tree,"phylogenetic_diversity/Lanzarote_16s.tre")
 # 
 # ### Annotate tree
 # taxonomy$label.long <- paste0(taxonomy$Phylum,"|",
 #                               taxonomy$Class,"|",
 #                               taxonomy$Order, "|",
 #                               taxonomy$Family,"|",
 #                               taxonomy$ID)
 # 
 # tree.names <- suppressWarnings(phylotools::sub.taxa.label(tree, taxonomy[c("ID", "label.long")]))
 # ape::write.tree(tree.names,"phylogenetic_diversity/Lanzarote_16s_names.tre")
 
 
########## S 3.2 - Calculate phylogenetic diversity  -----------------------------------------------------
 

 tree <- ape::read.tree("Lanzarote_16s.tre")
 
 comm_pa <- (comm > 0) * 1

 pd_sample <- BAT::alpha(comm_pa, tree = tree) 
 pd_sample <- as.data.frame(pd_sample)
 colnames(pd_sample) <- "phylo.diver"
 pd_sample$Sample <- rownames(pd_sample)
 
 stations <- merge(stations, pd_sample, by = "Sample")

 rm(pd_sample)
 
 
 ########## S 3.3 - Phylogenetic richness  -----------------------------------------------------
 
 plot(stations$phylo.diver ~ stations$richness)
 
 mod.phylo <- glm(phylo.diver ~ Type + reads + Distance.sea + log(buildings + 1),
                          data = stations,
                          family = Gamma(link="log"))
 
 performance::check_model(mod.phylo) 
 
 
 #performance::check_overdispersion(mod.phylo) 
 car::Anova(mod.phylo)
 summary(mod.phylo)
 
 emm_type <- emmeans::emmeans(mod.phylo, ~ Type)
 tukey_type <- pairs(emm_type, adjust = "tukey")
 summary(tukey_type)
 
 
 
 
 ########## S 3.4 - Phylogenetic composition  -----------------------------------------------------
 
 beta.phylo <- BAT::beta(comm_pa, tree = tree)


 mean(beta.phylo$Btotal)
 mean(beta.phylo$Brepl / beta.phylo$Btotal)
 mean(beta.phylo$Brich / beta.phylo$Btotal)
 
 (pb.phylo <- beta_density_plot(comm_pa, stations$Type, beta.phylo,
                                fill_values = c("Within habitat" = "#1b9e77",
                                                "Across habitats" = "#d95f02"),
                                alpha = 0.6))
 
 
 ## Total beta diversity
 model.Bphtotal <- beta.phylo$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads + richness
 
 (perm.phyT <- vegan::adonis2(model.Bphtotal, stations, permutations = 9999, by = "terms"))
 
 
 ## Beta nestness
 model.Bphnest <- beta.phylo$Brich ~  Type + Distance.sea + log(buildings + 1) + reads + richness
 (perm.phyRi <- vegan::adonis2(model.Bphnest, stations, permutations = 9999, by = "terms"))
 
 
 
 ## Beta turnover
 model.Bphturn <- beta.phylo$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads + richness
 (perm.phyRe <- vegan::adonis2(model.Bphturn, stations, permutations = 9999, by = "terms"))
 
 rm(model.Bphturn,model.Bphtotal,model.Bphturn,perm.phyRi,perm.phyRe,perm.phyT)
 
 
  
 
 
 
 
 
