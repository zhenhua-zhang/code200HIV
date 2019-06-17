module load R/3.3.3-foss-2015b-bare
Rscript 200HIV_qtl_mapping.R cknmlt -S plasma_marker.cfg
Rscript 200HIV_qtl_mapping.R trfm -S plasma_marker.cfg
Rscript 200HIV_qtl_mapping.R qtlmp -S plasma_marker.cfg
Rscript 200HIV_qtl_mapping.R pltmht -S plasma_marker.cfg
Rscript 200HIV_qtl_mapping.R pltdsg -S plasma_marker.cfg
