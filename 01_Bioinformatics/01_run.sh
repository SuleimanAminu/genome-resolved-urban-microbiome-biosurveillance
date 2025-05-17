#!/bin/bash

##########################################
# SLURM Job Configuration
##########################################

#SBATCH --nodes=1                     # Number of nodes to request
#SBATCH --ntasks=56                   # Number of tasks (max = nodes * 56)
#SBATCH --time=35:00:00               # Max walltime (HH:MM:SS)
#SBATCH --partition=compute           # Partition to use
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/out/slurm-%j.out2
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/err/slurm-%j.err2

# Load Conda
eval "$(conda shell.bash hook)"

##########################################
# Input Parameters and Path Definitions
##########################################

project=$1

SRA_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/DATA/$project
sh_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/used_scripts
collect=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/collect_results/$project

phix_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Shotgun-metagenomics/fastq_screen_contamination_genomes/PhiX/GCF_000819615.1_ViralProj14015_genomic.fna
bbmap_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/bbmap
adapters=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Shotgun-metagenomics/bbmap/resources
database=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Shotgun-metagenomics/bacteria_Ref_genomes/database
protein_nr=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Shotgun-metagenomics/protein_nr
diamond_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/shared/

##########################################
# Data Directory Setup & Dataset Download
##########################################

mkdir -p $SRA_path
cd $SRA_path

# Activate environment and download SRA project
conda activate pysradb
pysradb download -y -p $project

# Add SRA Toolkit to PATH
export PATH=$PATH:/home/user/sratoolkit.3.1.1-ubuntu64/bin

##########################################
# Sample Processing
##########################################

for ids in `less $SRA_path/SRA_ids`
do
    mkdir $SRA_path/$ids                   # Create folder for each run (only first time)
    cd $SRA_path/$ids

    # Download and check sample
    prefetch $ids -O .
    printf "\n $ids is OK \n\n"

    # Extract and clean directory names
    ls -d */ > ids2.txt
    sed -i 's/\///g' ids2.txt

    for fname in `less ids2.txt`
    do
        id=$(basename $fname ".sra")

        ##########################################
        # Step 1: Generate FASTQ files
        fasterq-dump $id
        printf "\n $ids/$id is OK \n\n"

        ##########################################
        # Step 2: Quality Control with FASTQC
        conda activate fastqc
        fastqc "$id"_1.fastq "$id"_2.fastq -f fastq -o ./
        printf "\n $ids : FASTQC is Done \n\n"

        ##########################################
        # Step 3: Read Preprocessing with BBTools

        # Reformat reads
        bash $bbmap_path/reformat.sh in1="$id"_1.fastq in2="$id"_2.fastq out1="$id"_1.processed.fq out2="$id"_2.processed.fq
        printf "\n $ids : reformat.sh is Done \n\n"

        # Merge paired-end reads
        bash $bbmap_path/bbmerge.sh in1="$id"_1.processed.fq in2="$id"_2.processed.fq out="$id"_merged.fq outu1="$id"_unmerged1.fq outu2="$id"_unmerged2.fq
        printf "\n '$id'_merged.fq : bbmerge.sh is Done \n\n"

        # Remove optical duplicates
        bash $bbmap_path/clumpify.sh in="$id"_merged.fq out="$id"_clumped.fq groups=16
        bash $bbmap_path/clumpify.sh in="$id"_unmerged1.fq out="$id"_clumped_unmerged1.fq groups=16
        bash $bbmap_path/clumpify.sh in="$id"_unmerged2.fq out="$id"_clumped_unmerged2.fq groups=16
        printf "\n '$id' clumpify.sh is Done \n\n"

        # Quality trim to Q20 and minimum length of 25
        bash $bbmap_path/bbduk.sh in="$id"_clumped.fq out="$id"_clean.fq qtrim=rl trimq=20 minlen=25
        bash $bbmap_path/bbduk.sh in="$id"_clumped_unmerged1.fq out="$id"_unmerged1_clean.fq qtrim=rl trimq=20 minlen=25
        bash $bbmap_path/bbduk.sh in="$id"_clumped_unmerged2.fq out="$id"_unmerged2_clean.fq qtrim=rl trimq=20 minlen=25
        printf "\n '$id' quality trimming is Done \n\n"

        # Adapter trimming if necessary (based on FASTQC output)
        bash $bbmap_path/bbduk.sh in="$id"_clean.fq out="$id"_clean_adapters.fq ref=$adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
        bash $bbmap_path/bbduk.sh in="$id"_unmerged1_clean.fq out="$id"_unmerged1_clean_adapters.fq ref=$adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
        bash $bbmap_path/bbduk.sh in="$id"_unmerged2_clean.fq out="$id"_unmerged2_clean_adapters.fq ref=$adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
        printf "\n $ids : Adapter trimming is Done \n\n"

        ##########################################
        # Step 4: Contaminant Screening & Mapping
        # You can choose one or both methods below depending on your setup

        sbatch $sh_path/fastq_screen.sh $ids $id $SRA_path $bbmap_path $sh_path $database $protein_nr $collect $diamond_path
        bash $sh_path/fastq_screen.sh $ids $id $SRA_path $bbmap_path $sh_path $database $protein_nr $collect $diamond_path

    done
done

