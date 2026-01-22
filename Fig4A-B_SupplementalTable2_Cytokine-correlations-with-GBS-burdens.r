## Spearman correlations - GBS lumen CFU & cytokines at PAIRED timepoints

library(readxl)
library(tidyr)
library(dplyr)


gbs_lumen_df <- read_excel('MUR-ALL-GBS-lumen-swab-CFU-combined.xlsx')
gbs_lumen_df <- gbs_lumen_df %>% 
  dplyr::filter(dpi %in% c(2,7)) %>%
  mutate(`unique-mouse-id` = paste0(`exp-id`,'_',`mouse-id`),
         timepoint=paste0(as.character(dpi),"dpi"))

cyto_df <- read_excel('mur-all-cyto-combined_v2.xlsx')

#make list of cytokine names
cytokine <- colnames(cyto_df[,9:ncol(cyto_df)])

#make list of cytokine timepoints to subset
cyto_tp_keep <- c("2dpi","7dpi")

#filter based on the timepoints
cyto_df_filt <- cyto_df %>% dplyr::filter(timepoint %in% cyto_tp_keep)

#determine which mice are shared between datasets
shared_mice <- intersect(gbs_lumen_df$'unique-mouse-id',cyto_df_filt$'unique-mouse-id')


#clean up cyto data in prep for df merging
#create individual columns for each cytokine-timepoint combination
cyto_df_wide <- cyto_df_filt %>% pivot_longer(
    cols= all_of(cytokine),
    names_to='cytokine',
    values_to='concentration') %>% mutate(
    cytokine_time=paste0(cytokine,"_",timepoint)) %>% select(
    'unique-mouse-id',cytokine_time,concentration) %>% pivot_wider(
    names_from = cytokine_time,
    values_from = concentration)

cyto_df_wide <- cyto_df_filt %>% pivot_longer(
    cols= all_of(cytokine),
    names_to='cytokine',
    values_to='concentration') 

#create individual columns for each gbs-timepoint combination
gbs_df_wide <- gbs_lumen_df %>% mutate(
    dpi=as.character(dpi)) %>% mutate(
    cfu_time=paste0('gbs_cfu','_d',dpi)) %>% select(
    'unique-mouse-id',cfu_time,'log-gbs-cfu-ml') %>% pivot_wider(
    names_from=cfu_time,
    values_from='log-gbs-cfu-ml')

gbs_df_wide <- gbs_lumen_df %>% mutate(
    dpi=as.character(timepoint))

#combine gbs and cyto observations for each mouse
combined_df <- inner_join(cyto_df_wide, gbs_df_wide, by=c('unique-mouse-id','timepoint'))

#make list of cytokine-timepoint combinations and gbs-timepoint combinations
allcols <- colnames(combined_df)

cyto_pattern <- paste(cytokine,collapse='|')
cyto_tp <- allcols[grepl(cyto_pattern, allcols)]

gbs_tp <- allcols[grepl("gbs*",allcols)]



cor_results <- data.frame(
  gbs_lumen_cfu_timepoint = character(),
  cytokine_and_timepoint = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (tp in unique(combined_df$timepoint)) {
  df_tp <- combined_df %>% filter(timepoint == tp) %>% filter(!is.na(gbs), !is.na(cytokine))
  if (nrow(df_tp) >= 5) {
    test <- cor.test(df_tp$gbs, df_tp$cytokine, method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      timepoint = tp,
      rho = test$estimate,
      p_value = test$p.value
    ))
  }
}

#adjust p-values for multiple testing using FDR
cor_results$p_value_adj <- p.adjust(cor_results$p_value, method = "fdr")



##SPEARMAN##
cor_results <- data.frame(
  gbs_lumen_cfu_timepoint = character(),
  cytokine_and_timepoint = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (gbs_col in gbs_tp) {
  for (cyt_col in cyto_tp) {
    df_pair <- combined_df[, c(gbs_col, cyt_col)]
    df_pair <- df_pair[complete.cases(df_pair), ]  # remove NAs
    
    if (nrow(df_pair) >= 3) {  # correlation needs at least 3 points
      test <- cor.test(df_pair[[1]], df_pair[[2]], method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        gbs = gbs_col,
        cytokine = cyt_col,
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

#export
write.table(cor_results)