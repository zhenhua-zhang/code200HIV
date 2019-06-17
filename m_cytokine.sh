module load R/3.3.3-foss-2015b-bare
Rscript 200HIV_qtl_mapping.R cknmlt -S m_cytokine.cfg
Rscript 200HIV_qtl_mapping.R trfm -S m_cytokine.cfg
Rscript 200HIV_qtl_mapping.R qtlmp -S m_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltmht -S m_cytokine.cfg
Rscript 200HIV_qtl_mapping.R pltdsg -S m_cytokine.cfg
