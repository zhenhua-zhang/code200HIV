#!/bin/bash
#SBATCH --time=2:00:0
#SBATCH --ntasks=1
#SBATCH --output=%j-%u-QTLmapping.log
#SBATCH --mem=40G
#SBATCH --cpus=2


#
## Main entry of the QTLmapping
#
source /apps/modules/modules.bashrc
module list
module load Python/3.6.3-foss-2015b
module list
source _qtlmapping_venv/bin/activate

my_path=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/200HIV
python mapper.py map \
	-w ${my_path}/outputs/HIVReservior/HIVReservior_AG \
	-p ${my_path}/inputs/datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv \
	-c ${my_path}/inputs/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
	-g ${my_path}/inputs/dosage/200HIV_dosages
