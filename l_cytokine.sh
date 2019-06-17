module load R/3.3.3-foss-2015b-bare
Rscript 200HIV_qtl_mapping.R cknmlt -S l_cytokine.cfg
Rscript 200HIV_qtl_mapping.R trfm -S l_cytokine.cfg
Rscript 200HIV_qtl_mapping.R qtlmp -S l_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltmht -S l_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltdsg -S l_cytokine.cfg
