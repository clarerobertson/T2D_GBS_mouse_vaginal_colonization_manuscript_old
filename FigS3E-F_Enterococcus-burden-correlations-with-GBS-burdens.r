library(readxl)
library(tidyr)
library(dplyr)


##GBS DATA##
# read in & clean GBS tissue data
tissue <- read_excel('MUR-ALL-GBS-tissue-CFU-combined.xlsx',
                    .name_repair='universal')

tissue$log.gbs.cfu.recovered <- as.numeric(tissue$log.gbs.cfu.recovered)

tissue_filt <- tissue %>% 
    dplyr::filter(gbs.strain == "nctc" ) %>% 
    dplyr::filter(intervention == 'mock' ) %>%
    dplyr::filter(!is.na(log.gbs.cfu.recovered))   

tissue_filt$diet <- tolower(tissue_filt$diet)
tissue_filt$exp.id <- tolower(tissue_filt$exp.id)
tissue_filt$mouse.id <- tolower(tissue_filt$mouse.id)

#to mouse.id, if there are only 2 characters, add a 0 as the second character
tissue_filt$mouse.id <- ifelse(nchar(tissue_filt$mouse.id) == 2, 
                                paste0(substr(tissue_filt$mouse.id, 1, 1), "0", substr(tissue_filt$mouse.id, 2, 2)), 
                                tissue_filt$mouse.id)
#if there are 4 characters, remove the second character but keep the 3rd and 4th
tissue_filt$mouse.id <- ifelse(nchar(tissue_filt$mouse.id) == 4, 
                                paste0(substr(tissue_filt$mouse.id, 1, 1), substr(tissue_filt$mouse.id, 3, 4)), 
                                tissue_filt$mouse.id)

#create the column unique.mouse.id by combining exp.id and mouse.id
tissue_filt <- tissue_filt %>%
    dplyr::mutate(unique.mouse.id = paste(exp.id, mouse.id, sep = "_"))


# read in & clean GBS swab data
swab <- read_excel('MUR-ALL-GBS-lumen-swab-CFU-combined.xlsx',
                    .name_repair='universal')
swab$log.gbs.cfu.ml <- as.numeric(swab$log.gbs.cfu.ml)

swab_filt <- swab %>% 
    dplyr::filter(gbs.strain == "nctc" ) %>% 
    dplyr::filter(intervention == 'mock' ) %>% 
    dplyr::filter(!is.na(log.gbs.cfu.ml))

swab_filt$diet <- tolower(swab_filt$diet)
swab_filt$exp.id <- tolower(swab_filt$exp.id)
swab_filt$mouse.id <- tolower(swab_filt$mouse.id)

#to mouse.id, if there are only 2 characters, add a 0 as the second character
swab_filt$mouse.id <- ifelse(nchar(swab_filt$mouse.id) == 2, 
                                paste0(substr(swab_filt$mouse.id, 1, 1), "0", substr(swab_filt$mouse.id, 2, 2)), 
                                swab_filt$mouse.id)
#if there are 4 characters, remove the second character but keep the 3rd and 4th
swab_filt$mouse.id <- ifelse(nchar(swab_filt$mouse.id) == 4, 
                                paste0(substr(swab_filt$mouse.id, 1, 1), substr(swab_filt$mouse.id, 3, 4)), 
                                swab_filt$mouse.id)

#create the column unique.mouse.id by combining exp.id and mouse.id
swab_filt <- swab_filt %>%
    dplyr::mutate(unique.mouse.id = paste(exp.id, mouse.id, sep = "_"))

            
##ENTERO DATA##
#read in & clean enterococcus tissue data
tissue_e <- read_excel('MUR-ALL-Enterococcus-tissue-CFU-combined.xlsx',
                    .name_repair='universal')

tissue_e$log.enterococcus.cfu.ml <- as.numeric(tissue_e$log.enterococcus.cfu.ml)

tissue_e_filt <- tissue_e %>% 
    dplyr::filter(gbs.strain == "nctc" ) %>% 
    dplyr::filter(intervention == 'mock' ) %>% 
    dplyr::filter(!is.na(log.enterococcus.cfu.ml))

tissue_e_filt$diet <- tolower(tissue_e_filt$diet)
tissue_e_filt$exp.id <- tolower(tissue_e_filt$exp.id)
tissue_e_filt$mouse.id <- tolower(tissue_e_filt$mouse.id)

#to mouse.id, if there are only 2 characters, add a 0 as the second character
tissue_e_filt$mouse.id <- ifelse(nchar(tissue_e_filt$mouse.id) == 2, 
                                paste0(substr(tissue_e_filt$mouse.id, 1, 1), "0", substr(tissue_e_filt$mouse.id, 2, 2)), 
                                tissue_e_filt$mouse.id)

#create the column unique.mouse.id by combining exp.id and mouse.id
tissue_e_filt <- tissue_e_filt %>%
    dplyr::mutate(unique.mouse.id = paste(exp.id, mouse.id, sep = "_"))


# read in & clean enterococcus swab data
swab_e <- read_excel('MUR-ALL-Enterococcus-lumen-swab-CFU-combined.xlsx',
                    .name_repair='universal')
swab_e$log.enterococcus.cfu.ml <- as.numeric(swab_e$log.enterococcus.cfu.ml)

swab_e_filt <- swab_e %>% 
    dplyr::filter(gbs.strain == "nctc" ) %>% 
    dplyr::filter(intervention == 'mock' ) %>% 
    dplyr::filter(!is.na(log.enterococcus.cfu.ml))

swab_e_filt$diet <- tolower(swab_e_filt$diet)
swab_e_filt$exp.id <- tolower(swab_e_filt$exp.id)
swab_e_filt$mouse.id <- tolower(swab_e_filt$mouse.id)

#to mouse.id, if there are only 2 characters, add a 0 as the second character
swab_e_filt$mouse.id <- ifelse(nchar(swab_e_filt$mouse.id) == 2, 
                                paste0(substr(swab_e_filt$mouse.id, 1, 1), "0", substr(swab_e_filt$mouse.id, 2, 2)), 
                                swab_e_filt$mouse.id)
#if there are 4 characters, remove the second character but keep the 3rd and 4th
swab_e_filt$mouse.id <- ifelse(nchar(swab_e_filt$mouse.id) == 4, 
                                paste0(substr(swab_e_filt$mouse.id, 1, 1), substr(swab_e_filt$mouse.id, 3, 4)), 
                                swab_e_filt$mouse.id)

#create the column unique.mouse.id by combining exp.id and mouse.id
swab_e_filt <- swab_e_filt %>%
    dplyr::mutate(unique.mouse.id = paste(exp.id, mouse.id, sep = "_"))
            
##MERGE##

#merge swab data on unique.mouse.id and timepoint
merged_swab <- dplyr::inner_join(
  swab_filt, swab_e_filt, by = c("unique.mouse.id", "dpi")) |>
  dplyr::select(unique.mouse.id, 
                log.gbs.cfu.ml, 
                log.enterococcus.cfu.ml, 
                dpi, 
                diet.x, diet.y)

#merge tissue data on unique.mouse.id and tissue

merged_tissue <- dplyr::inner_join(
    tissue_filt, tissue_e_filt, by = c("unique.mouse.id", "tissue")) |>
    dplyr::select(unique.mouse.id, 
    log.gbs.cfu.recovered, 
    log.enterococcus.cfu.ml, 
    tissue, 
    diet.x, diet.y)


##CORRELATIONS - TISSUE##


tissue_types <- c("Vagina", "Cervix","Uterus")

#initialize dataframes
tissue_results_df <- data.frame()

tissue_input_df <- data.frame()

#
for (t in tissue_types) {
  tissue_data <- merged_tissue[merged_tissue$tissue == t, ]
    tissue_input_df <- rbind(tissue_input_df, tissue_data)
  if (nrow(tissue_data) < 2) next #skip if too few data points
  
  #spearman#
  tissue_test <- cor.test(tissue_data$log.enterococcus.cfu.ml,
                           tissue_data$log.gbs.cfu.recovered,
                           method = "spearman")

  tissue_results_df <- rbind(tissue_results_df, data.frame(
    tissue = t,
    rho = tissue_test$estimate,
    p_value = tissue_test$p.value
  ))
}

#adjust p-values for multiple testing using FDR
tissue_results_df$p_value_adj <- p.adjust(tissue_results_df$p_value, method = "fdr")

#export results
write.table(tissue_results_df, 
"spearman_tissue_gbs_enterococcus_results.tsv", 
sep = "\t", row.names = FALSE) 


##CORRELATIONS - SWAB##

timepoints <- unique(swab_filt$dpi)

#initialize dataframes
swab_results_df <- data.frame()

swab_input_df <- data.frame()

#loop
for (tp in timepoints) {
  merged_swab_subset <- merged_swab[merged_swab$dpi == tp, ]
  swab_input_df <- rbind(swab_input_df, merged_swab_subset)
  if (nrow(merged_swab) < 2) next #skip if too few data points
  
  #spearman#
  swab_test <- cor.test(merged_swab_subset$log.enterococcus.cfu.ml,
                         merged_swab_subset$log.gbs.cfu.ml,
                         method = "spearman")

  #store results
  swab_results_df <- rbind(swab_results_df, data.frame(
    timepoint = tp,
    rho = swab_test$estimate,
    p_value = swab_test$p.value
  ))
}
#adjust p-values for multiple testing using FDR
swab_results_df$p_value_adj <- p.adjust(swab_results_df$p_value, method = "fdr")

#export results
write.table(swab_results_df, 
"spearman_swab_gbs_enterococcus_results.tsv", 
sep = "\t", row.names = FALSE) 