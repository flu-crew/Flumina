#!/bin/bash
#SBATCH --job-name=Flumina-test
#SBATCH -p priority --qos=vpru
#SBATCH -A nadc_iav
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=carl.hutter@usda.gov
#SBATCH -t 168:00:00
#SBATCH -e "stderr.%j.%N"
#SBATCH -D /Flumina_test/Flumina

date

module load miniconda/4.12.0

source activate /carl.hutter/conda/Flumina

bash Flumina config.cfg

date

# End of file
