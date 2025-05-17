#!/bin/bash
#SBATCH --job-name=model_comparison
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/suleiman.aminu/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/suleiman.aminu/Colorectal/err/slurm-%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 02_model_comparison.py

