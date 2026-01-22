library(readxl)
library(tidyr)
library(dplyr)


##GTT##
gbs_and_gtt_df <- read_excel('all data for GTT GBS correlations.xlsx')

#make list of cytokine-timepoint combinations and gbs-timepoint combinations
allcols <- colnames(gbs_and_gtt_df)

gtt_tp <- allcols[grepl("gtt*", allcols)]

gbs_tp <- allcols[grepl("gbs*",allcols)]


##SPEARMAN##
cor_results <- data.frame(
  gbs_lumen_cfu_timepoint = character(),
  gtt_timepoint = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (gbs_col in gbs_tp) {
  for (gtt_col in gtt_tp) {
    df_pair <- gbs_and_gtt_df[, c(gbs_col, gtt_col)]
    df_pair <- df_pair[complete.cases(df_pair), ]  # remove NAs
    
    if (nrow(df_pair) >= 3) {  # correlation needs at least 3 points
      test <- cor.test(df_pair[[1]], df_pair[[2]], method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        gbs = gbs_col,
        gtt = gtt_col,
        rho = test$estimate,
        p_value = test$p.value
      ))
    }
  }
}

cor_results$fdr <- p.adjust(cor_results$p_value, method = "fdr")

cor_results %>%
  arrange(p_value) %>%
  head(10)

sig_results <- cor_results %>% filter(p_value < 0.05)

top_results <- cor_results %>%
  arrange(p_value) %>%
  mutate(label = paste0(gtt, " vs ", gbs))

#export
write.table(cor_results,file='spearman-corr-results_gbs-lumen-cfu-and-every-cytokine_post-infection-timepoints.tsv', quote=FALSE, sep="\t")

write.table(top_results,file='spearman-corr-results_gbs-lumen-tissue-cfu-versus-GTT-AUC.tsv', quote=FALSE, sep="\t")



##BODY MASS##
bm <- read_excel('MUR-all-diabetes-markers-combined.xlsx',
                    .name_repair='universal')

bm$bm <- as.numeric(bm$bm) 

bm$diet <- tolower(bm$diet)
bm$exp.id <- tolower(bm$exp.id)
bm$mouse.id <- tolower(bm$mouse.id)

#to mouse.id, if there are only 2 characters, add a 0 as the second character
bm$mouse.id <- ifelse(nchar(bm$mouse.id) == 2, 
                                paste0(substr(bm$mouse.id, 1, 1), "0", substr(bm$mouse.id, 2, 2)), 
                                bm$mouse.id)
#if there are 4 characters, remove the second character but keep the 3rd and 4th
bm$mouse.id <- ifelse(nchar(bm$mouse.id) == 4, 
                                paste0(substr(bm$mouse.id, 1, 1), substr(bm$mouse.id, 3, 4)), 
                                bm$mouse.id)

#create the column unique.mouse.id by combining exp.id and mouse.id
bm <- bm %>%
    dplyr::mutate(unique.mouse.id = paste(exp.id, mouse.id, sep = "_"))


bm <- bm |> 
  dplyr::select(unique.mouse.id,weeks.on.diet,bm)

#calculate percent change from W0 to W12 timepoint for each unique.mouse.id
bm_pct_change <- bm %>%
  dplyr::filter(weeks.on.diet %in% c("W0", "W12")) %>%
  dplyr::group_by(unique.mouse.id) %>%
  dplyr::summarise(
    bm_w0 = bm[weeks.on.diet == "W0"],
    bm_w12 = bm[weeks.on.diet == "W12"]
  ) %>%
  dplyr::mutate(
    pct_change = (bm_w12 - bm_w0) / bm_w0 * 100
  )

#for loop - run spearman correlations between bm_pct_change and log.gbs.cfu.rml at all timepoints

 #initialize empty dataframes
swab_input_df <- data.frame() 
swab_results_df <- data.frame()

for (tp in unique(swab_filt$dpi)) {
  #subset
  swab_filt_subset <- swab_filt[swab_filt$dpi == tp, ]
  merged_swab_gbs_bm <- merge(bm_pct_change, 
  swab_filt_subset, 
  by = "unique.mouse.id") |>
  dplyr::select(unique.mouse.id, pct_change, log.gbs.cfu.ml) |>
  dplyr::filter(!is.na(pct_change), !is.na(log.gbs.cfu.ml))

#add timepoint column
  merged_swab_gbs_bm$timepoint <- tp

  swab_input_df <- rbind(swab_input_df, merged_swab_gbs_bm)

  #spearman
  swab_test <- cor.test(
    merged_swab_gbs_bm$pct_change, 
    merged_swab_gbs_bm$log.gbs.cfu.ml, 
    method = "spearman")

    print(paste("Spearman correlation for day", tp))
    print(swab_test)

    #store results
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
"spearman_swab_gbs_w12_body_mass_increase_results.tsv", 
sep = "\t", row.names = FALSE)



tissue_results_df <- data.frame()
tissue_input_df <- data.frame() 

for (tis in unique(tissue_filt$tissue)) {

  #subset
  tissue_filt_subset <- tissue_filt[tissue_filt$tissue == tis, ]
  merged_tissue_gbs_bm <- merge(bm_pct_change, 
  tissue_filt_subset, 
  by = "unique.mouse.id") |>
  dplyr::select(unique.mouse.id, pct_change, log.gbs.cfu.recovered) |>
  dplyr::filter(!is.na(pct_change), !is.na(log.gbs.cfu.recovered))
    
    tissue_input_df <- rbind(tissue_input_df, merged_tissue_gbs_bm)

  #spearman
  tissue_test <- cor.test(
    merged_tissue_gbs_bm$pct_change, 
    merged_tissue_gbs_bm$log.gbs.cfu.recovered, 
    method = "spearman")

    print(paste("Spearman correlation for tissue", tis))
    print(tissue_test)

    #store results
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
"spearman_tissue_gbs_d7_body_mass_increase_results.tsv", 
sep = "\t", row.names = FALSE)