library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)


#cytokines
cyto_df <- read_excel('mur-all-cyto-combined_v2.xlsx')

#metabolic markers
mark_df <- read_excel('MUR-all-diabetes-markers-combined.xlsx')


cytokine <- colnames(cyto_df[,9:ncol(cyto_df)])
marker <- colnames(mark_df[,7:ncol(mark_df)])

#data cleaning

#renaming columns
cyto_df <- cyto_df %>%
 rename(unique_sample_id = `unique-sample-id`)

mark_df <- mark_df %>%
    rename(unique_sample_id = `unique-datapoint-id`)

#timepoints
cyto_df$timepoint <- gsub("^0dpi$", "12wk", cyto_df$timepoint)

##unique sample ids
library(stringr)

cyto_df$unique_sample_id <- str_replace_all(
  cyto_df$unique_sample_id, "0dpi", "12wk")

mark_df$unique_sample_id <- str_replace_all(
  mark_df$unique_sample_id,
  "_([A-Za-z])(\\d{1,2})_",
  function(x) {
    letter <- tolower(str_match(x, "_([A-Za-z])")[,2])
    num <- sprintf("%02d", as.numeric(str_match(x, "_[A-Za-z](\\d{1,2})_")[,2]))
    paste0("_", letter, num, "_")
  }
)

#timepoint
mark_df$unique_sample_id <- gsub("W(\\d+)$", "\\1wk", mark_df$unique_sample_id)

unique(mark_df$unique_sample_id)
unique(cyto_df$unique_sample_id)

intersect(mark_df$unique_sample_id,cyto_df$unique_sample_id)
#45 samples

#create individual columns for each cytokine-timepoint combination
cyto_df_wide <- cyto_df %>% pivot_longer(
    cols= all_of(cytokine),
    names_to='cytokine',
    values_to='concentration') %>% mutate(
    cytokine_time=paste0(cytokine,"_",timepoint)) %>% select(
    'unique-mouse-id','unique_sample_id',timepoint,cytokine_time,concentration) %>% pivot_wider(
    names_from = cytokine_time,
    values_from = concentration)

#create a timepoint column
mark_df <- mark_df %>%
  mutate(timepoint = str_extract(unique_sample_id, "(?<=_)[^_]+$"))
    
mark_df$bm <- as.numeric(mark_df$bm)


#create individual columns for each marker-timepoint combination
mark_df_wide <- mark_df %>% mutate(
    timepoint=as.character(timepoint)) %>% mutate(
    marker_time=paste0('bm','_',timepoint)) %>% select(
    'unique-mouse-id','unique_sample_id',marker_time,timepoint,bm,gtt,bm_diff) %>% 
pivot_wider(
    names_from=marker_time,
    values_from=bm)

#combine gbs and cyto observations for each mouse
combined_df <- inner_join(cyto_df_wide, mark_df_wide, 
                          by=c('unique_sample_id','timepoint'))



###SPEARMAN CORR - ALL TIMEPOINTS PAIRED##


#make list of cytokine-timepoint combinations and marker-timepoint combinations
allcols <- colnames(combined_df)

cyto_pattern <- paste(cytokine,collapse='|')
cyto_tp <- allcols[grepl(cyto_pattern, allcols)]

mark_tp <- allcols[grepl("bm|gtt",allcols)]


##SPEARMAN##
cor_results <- data.frame(
  marker_and_timepoint = character(),
  cytokine_and_timepoint = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (mark_col in mark_tp) {
  for (cyt_col in cyto_tp) {
    df_pair <- combined_df[, c(mark_col, cyt_col)]
    df_pair <- df_pair[complete.cases(df_pair), ]  # remove NAs
    
    if (nrow(df_pair) >= 3) {  # correlation needs at least 3 points
      test <- cor.test(df_pair[[1]], df_pair[[2]], method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        marker = mark_col,
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

write.table(cor_results,
            file='cyto-metab-markers-corr/spearman-corr-results_gtt-bm-and-every-cytokine_PAIRED-timepoints.tsv', 
            quote=FALSE, sep="\t")




###SPEARMAN CORR - 12 WEEK ONLY##




mark_tp <- allcols[grepl("bm_12wk|gtt|bm_diff",allcols)]

##SPEARMAN##
cor_results <- data.frame(
  marker_and_timepoint = character(),
  cytokine_and_timepoint = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (mark_col in mark_tp) {
  for (cyt_col in cyto_tp) {
    df_pair <- combined_df[, c(mark_col, cyt_col)]
    df_pair <- df_pair[complete.cases(df_pair), ]  # remove NAs
    
    if (nrow(df_pair) >= 3) {  # correlation needs at least 3 points
      test <- cor.test(df_pair[[1]], df_pair[[2]], method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        marker = mark_col,
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

write.table(cor_results,
            file='cyto-metab-markers-corr/spearman-corr-results_gtt-bm-and-every-cytokine_12wk-only.tsv', 
            quote=FALSE, sep="\t")





###SPEARMAN CORR - 4- AND 8-WEEK CYTOKINES WITH 12-WEEK METABOLIC MARKERS##



#make mouse ids match
mark_df_wide$`unique-mouse-id` <- str_replace_all(
  mark_df_wide$`unique-mouse-id`,
  "_([A-Za-z])(\\d{1,2})",
  function(x) {
    letter <- tolower(str_match(x, "_([A-Za-z])")[,2])
    num <- sprintf("%02d", as.numeric(str_match(x, "_[A-Za-z](\\d{1,2})")[,2]))
    paste0("_", letter, num)
  }
)

#combine dataframes based on mouse-id
combined_df <- left_join(cyto_df_wide, mark_df_wide, 
                          by=c('unique-mouse-id'))


combined_df <- inner_join(cyto_df_wide, mark_df_wide, 
                          by=c('unique-mouse-id'))


mark_tp <- allcols[grepl("gtt|bm_diff",allcols)]

cyto_tp <- allcols[grepl(cyto_pattern, allcols)]
cyto_tp <- cyto_tp[grep("_(4wk|8wk)$",cyto_tp)]

##SPEARMAN##
cor_results <- data.frame(
  marker_and_timepoint = character(),
  cytokine_and_timepoint = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (mark_col in mark_tp) {
  for (cyt_col in cyto_tp) {
    df_pair <- combined_df[, c(mark_col, cyt_col)]
    df_pair <- df_pair[complete.cases(df_pair), ]  # remove NAs
    
    if (nrow(df_pair) >= 3) {  # correlation needs at least 3 points
      test <- cor.test(df_pair[[1]], df_pair[[2]], method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        marker = mark_col,
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
  head(50)

sig_results <- cor_results %>% filter(p_value < 0.05)

write.table(cor_results,
            file='spearman-corr-results_gtt-bm-12wk-and-every-cytokine-4wk-and-8wk.tsv', 
            quote=FALSE, sep="\t")




###SPEARMAN CORR - 12-WEEK METABOLIC MARKERS WITH 2DPI AND 7DPI CYTOKINES###



#make mouse ids match
mark_df_wide$`unique-mouse-id` <- str_replace_all(
  mark_df_wide$`unique-mouse-id`,
  "_([A-Za-z])(\\d{1,2})",
  function(x) {
    letter <- tolower(str_match(x, "_([A-Za-z])")[,2])
    num <- sprintf("%02d", as.numeric(str_match(x, "_[A-Za-z](\\d{1,2})")[,2]))
    paste0("_", letter, num)
  }
)


#combine dataframes based on mouse-id
combined_df <- left_join(cyto_df_wide, mark_df_wide, 
                          by=c('unique-mouse-id'))



mark_tp <- allcols[grepl("bm_12wk|gtt|bm_diff",allcols)]
cyto_tp <- allcols[grep("_(2dpi|7dpi)$",allcols)]



##SPEARMAN##
cor_results <- data.frame(
  marker_and_timepoint = character(),
  cytokine_and_timepoint = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (mark_col in mark_tp) {
  for (cyt_col in cyto_tp) {
    df_pair <- combined_df[, c(mark_col, cyt_col)]
    df_pair <- df_pair[complete.cases(df_pair), ]  # remove NAs
    
    if (nrow(df_pair) >= 3) {  # correlation needs at least 3 points
      test <- cor.test(df_pair[[1]], df_pair[[2]], method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        marker = mark_col,
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
  head(50)

sig_results <- cor_results %>% filter(p_value < 0.05)

write.table(cor_results,
            file='cyto-metab-markers-corr/spearman-corr-results_12wk-markers-with-post-infection-cytokines.tsv', 
            quote=FALSE, sep="\t")