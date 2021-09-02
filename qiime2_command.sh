#! /bin/bash

# qiime2 analysis
conda activate qiime2-2021.4;

# import the sequences
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /home/wangbeibei_amd/ibd/16s_data \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path pairend_sequences.qza;

# visualization (4m)
qiime demux summarize \
   --i-data pairend_sequences.qza \
   --o-visualization pairend_sequences.qzv;

# denoising usng dada2 ()
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs pairend_sequences.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 140 \
  --p-n-threads 20 \
  --o-table dada2_table.qza \
  --o-representative-sequences rep_seqs.qza \
  --o-denoising-stats denoising_stats.qza;

# visualization
qiime feature-table tabulate-seqs \
  --i-data rep_seqs.qza \
  --o-visualization rep_seqs.qzv;
 
qiime metadata tabulate \
  --m-input-file denoising_stats.qza \
  --o-visualization denoising_stats.qzv;

# downloading the naive bayes classifier for silva v4 region
wget https://data.qiime2.org/2021.8/common/silva-138-99-515-806-nb-classifier.qza;

qiime feature-classifier classify-sklearn \
  --i-reads rep_seqs.qza \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --p-n-jobs 20 \
  --o-classification taxonomy.qza;

# summarize tables at different levels
qiime taxa collapse \
  --i-table dada2_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table table-l2.qza;
  
qiime taxa collapse \
  --i-table dada2_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table table-l3.qza;

qiime taxa collapse \
  --i-table dada2_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table table-l4.qza;
  
qiime taxa collapse \
  --i-table dada2_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table table-l5.qza; 
  
qiime taxa collapse \
  --i-table dada2_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-l6.qza;
  
qiime taxa collapse \
  --i-table dada2_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table table-l7.qza;







