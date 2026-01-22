library(dplyr)
library(ggplot2)
library(Maaslin2)

#genus level
df_input_data_g = read.delim(file = 'mur-feature-table-genus-relative-for-maaslin.tsv', header = TRUE, sep = "\t",
                            row.names = 1,
                            stringsAsFactors = TRUE);

df_input_data_g <- t(df_input_data_g);

df_input_metadata = read.delim(file = 'metadata-mur-02-12_v3.tsv', header = TRUE, sep = "\t",
                            row.names = 1,
                            stringsAsFactors = TRUE)

#create a new column containing both diet and timepoint info
df_input_metadata <- df_input_metadata %>% mutate(diet_timepoint=(paste(diet,timepoint,sep="_")));

#set variables as factors
df_input_metadata$diet_timepoint = factor(df_input_metadata$diet_timepoint,
                                              levels = c("C_W00", "C_W04", "C_W08","C_W12", "H_W00", "H_W04", "H_W08", "H_W12"))
    
fit_data_g <- Maaslin2(
    df_input_data_g, df_input_metadata, '/maaslin/output/genus',
    fixed_effects = c('diet','timepoint','diet_timepoint'),
    random_effects = c('mouse.id'),
    standardize = FALSE,
    save_models = TRUE)
#no significant associations identified

#species level
df_input_data_s = read.delim(file = 'mur-feature-table-species-relative-for-maaslin.tsv', header = TRUE, sep = "\t",
                            row.names = 1,
                            stringsAsFactors = TRUE);

df_input_data_s <- t(df_input_data_s);

df_input_metadata = read.delim(file = 'metadata-mur-02-12_v3.tsv', header = TRUE, sep = "\t",
                            row.names = 1,
                            stringsAsFactors = TRUE)
  
#create a new column containing both diet and timepoint info
df_input_metadata <- df_input_metadata %>% mutate(diet_timepoint=(paste(diet,timepoint,sep="_")));

#set variables as factors
df_input_metadata$diet_timepoint = factor(df_input_metadata$diet_timepoint,
                                              levels = c("C_W00", "C_W04", "C_W08","C_W12", "H_W00", "H_W04", "H_W08", "H_W12"))
    
    
fit_data_s <- Maaslin2(
    df_input_data_s, df_input_metadata, '/maaslin/output/species',
    fixed_effects = c('diet','timepoint','diet_timepoint'),
    random_effects = c('mouse.id'),
    standardize = FALSE,
    save_models = TRUE)

#compare diet groups at individual timepoints

df_input_metadata_w12 <- df_input_metadata %>% filter(timepoint=="W12")

fit_data_s <- Maaslin2(
    df_input_data_s, df_input_metadata_w12, '/maaslin/output/species/week12',
    fixed_effects = c('diet'),
    random_effects = c('mouse.id'),
    standardize = FALSE,
    save_models = TRUE)

df_input_metadata_w08 <- df_input_metadata %>% filter(timepoint=="W08")

fit_data_s <- Maaslin2(
    df_input_data_s, df_input_metadata_w08, '/maaslin/output/species/week08',
    fixed_effects = c('diet'),
    random_effects = c('mouse.id'),
    standardize = FALSE,
    save_models = TRUE)

df_input_metadata_w04 <- df_input_metadata %>% filter(timepoint=="W04")

fit_data_s <- Maaslin2(
    df_input_data_s, df_input_metadata_w04, '/maaslin/output/species/week04',
    fixed_effects = c('diet'),
    random_effects = c('mouse.id'),
    standardize = FALSE,
    save_models = TRUE)

#12 week only, multiple taxonomic levels

#variables
feattable_to_use_a <- "12wkonly/feature-table-asv-relative.tsv"
feattable_to_use_s <- "12wkonly/feature-table-species-relative-7.tsv"
feattable_to_use_g <- "12wkonly/feature-table-genus-relative-5.tsv"

meta_in <- "metadata/metadata-mur-02-12_v4.tsv"

out_path_a <- "maaslin/LM_DA_asv_results.tsv"
out_path_s <- "maaslin/LM_DA_species_results.tsv"
out_path_g <- "maaslin/LM_DA_genus_results.tsv"

#run LM (substitute a,s,g)
df_input_data <- read.table(file = feattable_to_use_g, header = TRUE, sep = "\t",
                            row.names = 1,
                            stringsAsFactors = TRUE)
#from the fetaure table, transposed rows and columns, deleted unnecessary empty column that might cause an error with maaslin, and labeled the column with sample ids as "sampleid" to match the metadata file
df_input_metadata <- read.table(file = meta_in, header = TRUE, sep = "\t",
                            row.names = 1,
                            stringsAsFactors = FALSE)

#create a metadata column for unique mouse id to include the mouse as a random effect
df_input_metadata |> mutate(df_input_metadata, 'unique-mouse-id'=)
    
fit_data <- Maaslin2(
    df_input_data, df_input_metadata, out_path_g,
    fixed_effects = c('diet'),
    random_effects = c('')
    standardize = FALSE,
    save_models = TRUE)
#no significant associations for any taxonomic level