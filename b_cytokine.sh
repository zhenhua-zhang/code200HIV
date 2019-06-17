module load R/3.3.3-foss-2015b-bare
Rscript 200HIV_qtl_mapping.R cknmlt -S b_cytokine.cfg
Rscript 200HIV_qtl_mapping.R trfm -S b_cytokine.cfg
Rscript 200HIV_qtl_mapping.R qtlmp -S b_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltmht -S b_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltdsg -S b_cytokine.cfg
