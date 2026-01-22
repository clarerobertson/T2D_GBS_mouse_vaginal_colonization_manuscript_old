#import sequences
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path mur-16s/02-12/manifests/manifest-mur-02-12_v4.tsv \
--output-path paired-end-demux_mur-02-12.qza \
--input-format PairedEndFastqManifestPhred33V2

#these are demultiplexed paired end sequences (a forward read + a reverse read for each sample) with phred 33

qiime demux summarize \
  --i-data paired-end-demux_mur-02-12.qza \
  --o-visualization paired-end-demux_mur-02-12.qzv
#checked summary of data to determine what cutoffs to use in denoising step

#denoise the reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux_mur-02-12.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 298 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --o-representative-sequences repseqs/mur-02-12_rep-seqs.qza \
  --o-table tables/mur-02-12_feature-table.qza \
  --o-denoising-stats output/mur-02-12_denoising-stats.qza \

#tabulate the denoising stats
qiime metadata tabulate \
  --m-input-file output/mur-02-12_denoising-stats.qza \
  --o-visualization viz/mur-02-12_denoising-stats-summ.qzv

#some samples were sequenced more than once due to low amplification during first run
#group the feature table by the sample id, summing all reads that came from the same sample 

qiime feature-table group \
--i-table tables/mur-02-12_feature-table.qza \
--p-axis sample \
--m-metadata-file /mur-16s/02-12_seq_Aug2025_analysis/metadata/metadata-mur-02-12_v2.tsv \
--m-metadata-column mouse-id-full \
--p-mode sum \
--o-grouped-table tables/mur-02-12_feattable_combined-runs-per-sample.qza

qiime tools export \
--input-path tables/mur-02-12_feattable_combined-runs-per-sample.qza \
--output-path exported/mur-02-12_feattable_combined-runs-per-sample

biom convert \
-i exported/mur-02-12_feattable_combined-runs-per-sample/feature-table.biom \
-o exported/mur-02-12_feattable_combined-runs-per-sample/feature-table.tsv --to-tsv

#map reads to greengenes2 database

qiime greengenes2 non-v4-16s \
--i-table tables/mur-02-12_feattable_combined-runs-per-sample.qza \
--i-sequences repseqs/mur-02-12_rep-seqs.qza \
--i-backbone /greengenes2/2022.10.backbone.full-length.fna.qza \
--o-mapped-table tables/mur-02-12_feattable_gg2.qza \
--o-representatives repseqs/mur-02-12_repseqs_gg2.qza

#assign taxonomy

qiime greengenes2 taxonomy-from-table \
--i-reference-taxonomy /greengenes2/2022.10.taxonomy.md5.nwk.qza \
--i-table tables/mur-02-12_feattable_gg2.qza \
--o-classification tables/mur-02-12.gg2.tabletaxonomymd5.qza

#--prevalence/abundance filtering & contaminant removal--#

#look at feature table to see what ASVs are in blanks
biom convert \
-i exported/mur-02-12_feattable_gg2/feature-table.biom \
-o exported/mur-02-12_feattable_gg2/ASV_table.txt --to-tsv --header-key taxonomy

#Create a feature table with species level taxonomy
qiime taxa collapse \
--i-table tables/mur-02-12_feattable_gg2.qza \
--i-taxonomy tables/mur-02-12.gg2.tabletaxonomymd5.qza \
--p-level 7 \
--o-collapsed-table tables/mur-02-12_feattable_species_gg2_withcontaminants.qza
#export & convert to tsv
qiime tools export \
--input-path tables/mur-02-12_feattable_species_gg2_withcontaminants.qza \
--output-path exported/mur-02-12_feattable_species_withcontaminants
biom convert \
-i exported/mur-02-12_feattable_species_withcontaminants/feature-table.biom \
-o exported/mur-02-12_feattable_species_withcontaminants/feature-table-species_withcontaminants.tsv --to-tsv

#summarize feature table
qiime feature-table summarize-plus \
  --i-table tables/mur-02-12_feattable_species_gg2_withcontaminants.qza \
  --o-feature-frequencies tables/mur-02-12_feattable_gg2_feature-frequencies.qza \
  --o-sample-frequencies tables/mur-02-12_feattable_gg2_sample-frequencies.qza \
  --o-summary viz/mur-02-12_feattable_species_gg2_visual_summary.qzv

#check which species are in blanks and whether any inappropriate species are present in samples

#keep only taxa found in at least 10 samples
qiime feature-table filter-features \
  --i-table tables/mur-02-12_feattable_gg2.qza \
  --p-min-samples 10 \
  --o-filtered-table tables/mur-02-12_feattable_gg2_min-samples10.qza ;
qiime tools export \
--input-path tables/mur-02-12_feattable_gg2_min-samples10.qza \
--output-path exported/mur-02-12_feattable_gg2_min-samples10 ;
biom convert \
-i exported/mur-02-12_feattable_gg2_min-samples10/feature-table.biom \
-o exported/mur-02-12_feattable_gg2_min-samples10/feature-table_min-samples10.tsv --to-tsv

#remove contaminants
qiime feature-table filter-features \
--i-table tables/mur-02-12_feattable_gg2_min-samples10.qza \
--m-metadata-file metadata/asvs_for_decontam.tsv \
--p-exclude-ids True \
--p-filter-empty-samples True \
--o-filtered-table tables/mur-02-12_feattable_gg2.qza ;
qiime tools export \
--input-path tables/mur-02-12_feattable_gg2.qza \
--output-path exported/mur-02-12_feattable_gg2 ;
biom convert \
-i exported/mur-02-12_feattable_gg2/feature-table.biom \
-o exported/mur-02-12_feattable_gg2/feature-table.tsv --to-tsv


#--generating counts tables--#

#collapse to species level
qiime taxa collapse \
--i-table tables/mur-02-12_feattable_gg2.qza \
--i-taxonomy tables/mur-02-12.gg2.tabletaxonomymd5.qza \
--p-level 7 \
--o-collapsed-table tables/mur-02-12_feattable_species_gg2.qza

#collapse to genus level 
qiime taxa collapse \
--i-table tables/mur-02-12_feattable_gg2.qza \
--i-taxonomy tables/mur-02-12.gg2.tabletaxonomymd5.qza \
--p-level 6 \
--o-collapsed-table tables/mur-02-12_feattable_genus_gg2.qza


#--generating relative abundance tables--#

##species
#turn feature table into "relative frequency table"
qiime feature-table relative-frequency \
--i-table tables/mur-02-12_feattable_species_gg2.qza \
--o-relative-frequency-table tables/mur-02-12_feattable_species_relative.qza
##genus
#turn feature table into "relative frequency table"
qiime feature-table relative-frequency \
--i-table tables/mur-02-12_feattable_genus_gg2.qza \
--o-relative-frequency-table tables/mur-02-12_feattable_genus_relative.qza


#--generating presence/absence tables--#
#variables
ft_s=tables/mur-02-12_feattable_species_gg2.qza
ft_g=tables/mur-02-12_feattable_genus_gg2.qza
ft_a=tables/mur-02-12_feattable_gg2.qza

pa_out_s=tables/mur-02-12_feattable_species_pres_abs_gg2.qza
pa_out_g=tables/mur-02-12_feattable_genus_pres_abs_gg2.qza
pa_out_a=tables/mur-02-12_feattable_pres_abs_gg2.qza

export_out_s=exported/mur-02-12_feattable_species_pres_abs_gg2
export_out_g=exported/mur-02-12_feattable_genus_pres_abs_gg2
export_out_a=exported/mur-02-12_feattable_pres_abs_gg2

convert_in_s=exported/mur-02-12_feattable_species_pres_abs_gg2/feature-table.biom
convert_in_g=exported/mur-02-12_feattable_genus_pres_abs_gg2/feature-table.biom
convert_in_a=exported/mur-02-12_feattable_pres_abs_gg2/feature-table.biom

convert_out_s=exported/mur-02-12_feattable_species_pres_abs_gg2/feature-table.tsv
convert_out_g=exported/mur-02-12_feattable_genus_pres_abs_gg2/feature-table.tsv
convert_out_a=exported/mur-02-12_feattable_pres_abs_gg2/feature-table.tsv

#run (substitute a,s,g)
qiime feature-table presence-absence \
--i-table $ft_s \
--o-presence-absence-table $pa_out_s ;

qiime tools export \
--input-path $pa_out_s \
--output-path $export_out_s ;

biom convert \
-i $convert_in_s \
-o $convert_out_s --to-tsv

##filter, keep 12 week only - genus, species, asv level##

#variables
filt_out_s=tables/mur-02-12_feattable_species_pres_abs_12weeksonly_gg2.qza
filt_out_g=tables/mur-02-12_feattable_genus_pres_abs_12weeksonly_gg2.qza
filt_out_a=tables/mur-02-12_feattable_pres_abs_12weeksonly_gg2.qza

export_out_s=exported/mur-02-12_feattable_species_pres_abs_12weeksonly_gg2
export_out_g=exported/mur-02-12_feattable_genus_pres_abs_12weeksonly_gg2
export_out_a=exported/mur-02-12_feattable_pres_abs_12weeksonly_gg2

convert_in_s=exported/mur-02-12_feattable_species_pres_abs_12weeksonly_gg2/feature-table.biom
convert_in_g=exported/mur-02-12_feattable_genus_pres_abs_12weeksonly_gg2/feature-table.biom
convert_in_a=exported/mur-02-12_feattable_pres_abs_12weeksonly_gg2/feature-table.biom

convert_out_s=exported/mur-02-12_feattable_species_pres_abs_12weeksonly_gg2/feature-table.tsv
convert_out_g=exported/mur-02-12_feattable_genus_pres_abs_12weeksonly_gg2/feature-table.tsv
convert_out_a=exported/mur-02-12_feattable_pres_abs_12weeksonly_gg2/feature-table.tsv

metadata_ft_filt=metadata/metadata-mur-02-12_v4.tsv

#run
qiime feature-table filter-samples \
  --i-table $pa_out_a \
  --m-metadata-file $metadata_ft_filt \
  --p-where '[timepoint]='"'"'W12'"'"'' \
  --o-filtered-table $filt_out_a

qiime tools export \
--input-path $filt_out_a \
--output-path $export_out_a ;

biom convert \
-i $convert_in_a \
-o $convert_out_a --to-tsv


#--calculate diversity metrics--#

qiime diversity core-metrics-phylogenetic \
    --i-table tables/mur-02-12_feattable_gg2.qza  \
    --i-phylogeny /2022.10.phylogeny.asv.nwk.qza  \
    --p-sampling-depth 1 \
    --m-metadata-file metadata/metadata-mur-02-12_v3.tsv  \
    --output-dir diversity/all-metrics-phylogenetic-2 \
    --verbose
qiime diversity core-metrics \
    --i-table tables/mur-02-12_feattable_gg2.qza  \
    --p-sampling-depth 1 \
    --m-metadata-file metadata/metadata-mur-02-12_v3.tsv  \
    --output-dir diversity/all-metrics-non-phylogenetic \
    --verbose

#alpha diversity for week 12 timepoint only
#filter - genus, species, asv level
qiime feature-table filter-samples \
  --i-table tables/mur-02-12_feattable_genus_gg2.qza \
  --m-metadata-file metadata/metadata-mur-02-12_v3.tsv \
  --p-where '[timepoint]='"'"'W12'"'"'' \
  --o-filtered-table tables/mur-02-12_feattable_genus_gg2_week12only.qza 

qiime feature-table filter-samples \
  --i-table tables/mur-02-12_feattable_species_gg2.qza \
  --m-metadata-file metadata/metadata-mur-02-12_v3.tsv \
  --p-where '[timepoint]='"'"'W12'"'"'' \
  --o-filtered-table tables/mur-02-12_feattable_species_gg2_week12only.qza 

qiime feature-table filter-samples \
  --i-table tables/mur-02-12_feattable_gg2.qza \
  --m-metadata-file metadata/metadata-mur-02-12_v3.tsv \
  --p-where '[timepoint]='"'"'W12'"'"'' \
  --o-filtered-table tables/mur-02-12_feattable_gg2_week12only.qza 
qiime diversity core-metrics-phylogenetic \
    --i-table tables/mur-02-12_feattable_gg2_week12only.qza  \
    --i-phylogeny /2022.10.phylogeny.asv.nwk.qza  \
    --p-sampling-depth 1 \
    --m-metadata-file metadata/metadata-mur-02-12_v4.tsv  \
    --output-dir diversity/all-metrics-phylogenetic-week12only \
    --verbose


#--differential abundance ANCOM-BC--#


#add pseudocount
qiime composition add-pseudocount \
--i-table tables/mur-02-12_feattable_gg2_decontam.qza \
--o-composition-table tables/comp_mur-02-12_feattable_genus_gg2.qza

#run ancom based on diet and timepoint

qiime composition ancombc \
--i-table tables/mur-02-12_feattable_genus_gg2.qza \
--m-metadata-file metadata/metadata-mur-02-12_v3.tsv \
--p-formula 'diet + timepoint' \
--o-differentials ancom/differentials_diet-timepoint.qza \
--verbose

qiime composition da-barplot \
--i-data ancom/differentials_diet-timepoint.qza \
--o-visualization viz/ancom-bc-genus

#run ancom based on diet BY timepoint

qiime composition ancombc \
--i-table tables/mur-02-12_feattable_genus_gg2.qza \
--m-metadata-file metadata/metadata-mur-02-12_v3.tsv \
--p-formula 'diet*timepoint' \
--o-differentials ancom/differentials_diet-by-timepoint.qza \
--verbose

qiime composition da-barplot \
--i-data ancom/differentials_diet-by-timepoint.qza \
--o-visualization viz/ancom-bc-genus_diet-by-timepoint