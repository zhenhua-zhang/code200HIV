module load R/3.3.3-foss-2015b-bare
Rscript 200HIV_qtl_mapping.R cknmlt -S c_cytokine.cfg
Rscript 200HIV_qtl_mapping.R trfm -S c_cytokine.cfg
Rscript 200HIV_qtl_mapping.R qtlmp -S c_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltmht -S c_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltdsg -S c_cytokine.cfg
