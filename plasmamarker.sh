module load R/3.3.3-foss-2015b-bare

Rscript 200HIV_qtl_mapping_v0.0.2.R cknmlt -S plasmamarker.cfg
Rscript 200HIV_qtl_mapping_v0.0.2.R trfm -S plasmamarker.cfg
Rscript 200HIV_qtl_mapping_v0.0.2.R qtlmp -S plasmamarker.cfg
Rscript 200HIV_qtl_mapping_v0.0.2.R qtlrpt -S plasmamarker.cfg
