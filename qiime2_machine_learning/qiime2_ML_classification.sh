#qiime2 sample-classifier plugin to predict archaea status

qiime sample-classifier predict-classification \
  --i-table genus_table_cohort2.qza \
  --i-sample-estimator sample_estimator_100_RF.qza \
  --o-predictions predictions__100_RF_final.qza\
  --output-dir unseen_classification__100_RF_final

qiime sample-classifier confusion-matrix \
  --i-predictions predictions__100_RF_final.qza \
  --m-truth-file meta_combined.csv \
  --m-truth-column Archaea_final \
  --o-visualization unseen_confusion_matrix__100_RF_final.qzv