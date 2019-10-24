#!/bin/bash

#
## Main entry of the QTLmapping
#
module load Python/3.6.3-foss-2015b
source ../_qtlmapping_venv/bin/activate

python mapper.py map \
	-w ${HOME}/Documents/projects/200HIV/outputs/HIVReservior/HIVReservior_AG \
	-p ${HOME}/Documents/projects/200HIV/inputs/datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv \
	-c ${HOME}/Documents/projects/200HIV/inputs/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
	-g ${HOME}/Documents/projects/200HIV/inputs/dosage/200HIV_dosages
