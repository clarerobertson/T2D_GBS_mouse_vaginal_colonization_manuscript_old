library(readxl)
library(tidyr)
library(dplyr)


#load metadata containing correlates
meta <- read_excel('mur-metadata.xlsx',
                    .name_repair='universal')

#to mouse.id, if there are only 2 characters, add a 0 as the second character
meta$mouse.id <- ifelse(nchar(meta$mouse.id) == 2, 
                                paste0(substr(meta$mouse.id, 1, 1), "0", substr(meta$mouse.id, 2, 2)), 
                                meta$mouse.id)

#create column unique.mouse.id with mouse.id appended to end of expt

meta <- meta %>%
    dplyr::mutate(unique.mouse.id = paste0(expt, "_", mouse.id))

meta$unique.mouse.id <- tolower(meta$unique.mouse.id)


#v.lumen.gluc.nmol.ml.w12
#v.tissue.gluc.mg.g.d7
#c.tissue.gluc.mg.g.d7
#u.tissue.gluc.mg.g.d7

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

#count number of non-na values in v.lumen.gluc.nmol.ml.w12
vaginal_lumen_glucose_count <- sum(!is.na(meta$v.lumen.gluc.nmol.ml.w12))
vaginal_lumen_glucose_count


vaginal_lumen_glucose <- meta |> 
  dplyr::select(unique.mouse.id, v.lumen.gluc.nmol.ml.w12)

#for loop - run spearman correlations between vaginal lumen glucose and log.gbs.cfu.recovered at all timepoints
swab_results_df <- data.frame()
#swab_input_df <- data.frame

for (tp in unique(swab_filt$dpi)) {
  #subset to current timepoint
  swab_filt_subset <- swab_filt[swab_filt$dpi == tp, ]
  merged_swab_gbs_gluc <- merge(vaginal_lumen_glucose, 
  swab_filt_subset, 
  by = "unique.mouse.id") |>
  dplyr::select(unique.mouse.id, v.lumen.gluc.nmol.ml.w12, log.gbs.cfu.ml) |>
  dplyr::filter(!is.na(v.lumen.gluc.nmol.ml.w12), !is.na(log.gbs.cfu.ml))

    #append input data
    #swab_input_df <- rbind(swab_input_df, merged_swab_gbs_gluc)

    
  #SPEARMAN
  swab_test <- cor.test(
    merged_swab_gbs_gluc$v.lumen.gluc.nmol.ml.w12, 
    merged_swab_gbs_gluc$log.gbs.cfu.ml, 
    method = "spearman")

    print(paste("Spearman correlation for day", tp))
    print(swab_test)

  swab_results_df <- rbind(swab_results_df, data.frame(
    timepoint = tp,
    rho = swab_test$estimate,
    p_value = swab_test$p.value
  ))
}
#adjust p-values for multiple testing using FDR
swab_results_df$p_value_adj <- p.adjust(swab_results_df$p_value, method = "fdr")

swab_results_df

#export results
write.table(swab_results_df, 
"spearman_swab_gbs_w12_vaginal_lumen_glucose_results.tsv", 
sep = "\t", row.names = FALSE) 


###tissues###
tissue_results_df <- data.frame()

for (tis in unique(tissue_filt$tissue)) {

#generate tissue initial
  tis_initial <- substr(tis, 1, 1)
#pull tissue glucose from metadata
  tis_gluc_col <- tolower(paste0(tis_initial, ".tissue.gluc.mg.g.d7"))

#subset data based on timepoint
  tissue_filt_subset <- tissue_filt[tissue_filt$tissue == tis, ]
  merged_tissue_gbs_gluc <- merge(meta, 
  tissue_filt_subset, 
  by = "unique.mouse.id") |>
  dplyr::select(unique.mouse.id, !!sym(tis_gluc_col), log.gbs.cfu.recovered) |>
  dplyr::filter(!is.na(!!sym(tis_gluc_col)), !is.na(log.gbs.cfu.recovered))

  #spearman#
  tissue_test <- cor.test(
    merged_tissue_gbs_gluc[[tis_gluc_col]], 
    merged_tissue_gbs_gluc$log.gbs.cfu.recovered, 
    method = "spearman")

    print(paste("Spearman correlation for tissue", tis))
    print(tissue_test)

  tissue_results_df <- rbind(tissue_results_df, data.frame(
    tissue = tis,
    rho = tissue_test$estimate,
    p_value = tissue_test$p.value
  ))
}
#adjust p-values for multiple testing using FDR
tissue_results_df$p_value_adj <- p.adjust(tissue_results_df$p_value, method = "fdr")

tissue_results_df

#export results
write.table(tissue_results_df, 
"spearman_tissue_gbs_d7_tissue_glucose_results.tsv", 
sep = "\t", row.names = FALSE) 