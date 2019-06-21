module load R/3.3.3-foss-2015b-bare

Rscript 200HIV_qtl_mapping_v0.0.2.R cknmlt -S c_cytokine_lg10tolg2_AG.cfg
Rscript 200HIV_qtl_mapping_v0.0.2.R trfm -S c_cytokine_lg10tolg2_AG.cfg
Rscript 200HIV_qtl_mapping_v0.0.2.R qtlmp -S c_cytokine_lg10tolg2_AG.cfg
Rscript 200HIV_qtl_mapping_v0.0.2.R qtlrpt -S c_cytokine_lg10tolg2_AG.cfg
