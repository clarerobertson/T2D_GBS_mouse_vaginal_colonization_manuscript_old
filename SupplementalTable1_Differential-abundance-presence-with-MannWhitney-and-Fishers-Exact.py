##--differential abundance Mann-Whitney U test--#
mann-whitney U test
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

#variables
feattable_to_use_s = "python/feature-table-species-relative-7.tsv"
feattable_to_use_g = "python/feature-table-genus-relative-5.tsv"

out_path_s = "python/mannwhitney_DA_species_results.tsv"
out_path_g = "python/mannwhitney_DA_genus_results.tsv"

#feature table
features = pd.read_csv(feattable_to_use_g, sep="\t", comment="#", index_col=0)
#transpose
features = features.T

#metadata
meta_fp = "metadata/metadata-mur-02-12_v4.tsv"
metadata = pd.read_csv(meta_fp, sep="\t")

#merge feature table with metadata table
df = metadata.merge(features, left_on="id", right_index=True)

#determine diet groups
diet_groups = df["diet"].unique()
if len(diet_groups) != 2:
    raise ValueError(f"Expected exactly 2 diet groups, found: {diet_groups}")
group1, group2 = diet_groups
print(f"Comparing diet groups: {group1} vs {group2}")

#run mann whitney
results = []
for feature in features.columns:
    vals1 = df[df["diet"] == group1][feature]
    vals2 = df[df["diet"] == group2][feature]
    
    stat, p = mannwhitneyu(vals1, vals2, alternative="two-sided")
    results.append((feature, stat, p))

#convert to df
results_df = pd.DataFrame(results, columns=["Feature", "U_statistic", "p_value"])

#fdr correction
results_df["p_adj"] = multipletests(results_df["p_value"], method="fdr_bh")[1]

#sort based on p value
results_df = results_df.sort_values("p_adj")
print(results_df.head())

#export
results_df.to_csv(out_path_g, sep="\t", index=False)


##--differential presence Fisher's Exact test--#


import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

#variables
feattable_to_use_s = "python/feature-table-species-pa.tsv" #presence absence 
feattable_to_use_g = "python/feature-table-genus-pa.tsv" #presence absence
feattable_to_use_a = "python/feature-table-asv-pa.tsv" #presence absence

out_path_s = "python/fishers_DP_species_results.tsv"
out_path_g = "python/fishers_DP_genus_results.tsv"
out_path_a = "python/fishers_DP_asv_results.tsv"

#feature table
features = pd.read_csv(feattable_to_use_a, sep="\t", comment="#",index_col=0)
#transpose
features = features.T

#metadata
meta_fp = "metadata/metadata-mur-02-12_v4.tsv"
metadata = pd.read_csv(meta_fp, sep="\t")

#merge feature table with metadata table
df = metadata.merge(features, left_on="id", right_index=True)

#determine diet groups
diet_groups = df["diet"].unique()
if len(diet_groups) != 2:
    raise ValueError(f"Expected exactly 2 diet groups, found: {diet_groups}")
group1, group2 = diet_groups
print(f"Comparing diet groups: {group1} vs {group2}")

#run test
results = []
    
for feature in features.columns:
    g1_present = (df[df["diet"] == group1][feature] > 0).sum()
    g1_absent  = (df[df["diet"] == group1][feature] == 0).sum()
    g2_present = (df[df["diet"] == group2][feature] > 0).sum()
    g2_absent  = (df[df["diet"] == group2][feature] == 0).sum()
    contingency = [[g1_present, g1_absent],
                    [g2_present, g2_absent]]
    oddsratio, p = fisher_exact(contingency, alternative="two-sided")
    results.append((feature, oddsratio, p))
    
#convert to df
results_df = pd.DataFrame(results, columns=["feature", "odds_ratio", "p_value"])

#fdr correction
results_df["p_adj"] = multipletests(results_df["p_value"], method="fdr_bh")[1]

#sort based on p value
results_df = results_df.sort_values("p_adj")
print(results_df.head())

#export
results_df.to_csv(out_path_a, sep="\t", index=False)