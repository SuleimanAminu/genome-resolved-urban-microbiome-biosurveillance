#!/bin/bash

################################################################################
# SLURM Job Configuration
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/out/slurm-%j.out2
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/err/slurm-%j.err2

# Activate Conda environment system
eval "$(conda shell.bash hook)"

################################################################################
# INPUTS AND DIRECTORY SETUP
################################################################################

project=$1

SRA_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/DATA/$project
collect=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/collect_AMR_Env/$project
sh_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/used_scripts

# Create result collection subfolders
mkdir -p $collect/{FQS_Before,FQS_After,quast,Bin,Bin1,Process,Process2,bracken_species,bracken_genus,bracken_family,CARD,VFDB}

fastq_screen_Before=$collect/FQS_Before
fastq_screen_After=$collect/FQS_After
quast_result=$collect/quast
Bin_result=$collect/Bin
Bin_result1=$collect/Bin1
Process_screen=$collect/Process
Process_screen2=$collect/Process2
bracken1=$collect/bracken_species
bracken2=$collect/bracken_genus
bracken3=$collect/bracken_family
CARD=$collect/CARD
VFDB=$collect/VFDB

################################################################################
# MAIN COLLECTION LOOP ACROSS SAMPLES
################################################################################

for ids in $(< $SRA_path/SRA_ids); do
    cd $SRA_path/$ids

    for fname in $(< ids2.txt); do
        id=$(basename $fname)

        echo "$id" >> $collect/SRA_ids

        ###########################################################
        # STEP 1: FASTQ Screen Results Before and After Filtering
        ###########################################################

        cp "$SRA_path/$ids/${id}_clean_adapters_screen.txt" $fastq_screen_Before
        cp "$SRA_path/$ids/${id}_clean_adapters.tagged.tagged_filter_screen.txt" $fastq_screen_After

        # ------------------------- BEFORE ------------------------
        cd $fastq_screen_Before
        ls > ids2.txt
        sed -i 's/\///g' ids2.txt

        for fname in $(< ids2.txt); do
            SRR=$(basename $fname "_clean_adapters_screen.txt")
            awk -F "\t" -v org="$SRR" '{print org"\t"$1"\t"$2"\t"$3}' "$SRR"_clean_adapters_screen.txt > "$Process_screen/${SRR}_screen.txt"
            sed -i '1d;1d;s/%//g' "$Process_screen/${SRR}_screen.txt"
        done

        cd $Process_screen
        for screen_file in *_screen.txt; do
            SRR=$(basename $screen_file "_screen.txt")
            while IFS= read -r ids; do
                printf "$ids\t" >> "${SRR}_FQS_before.txt"
            done < "$screen_file"
        done

        sed -i -e '$a\\' *_FQS_before.txt
        cat *_FQS_before.txt > $fastq_screen_Before/Fastqscreen_before.csv

        # ------------------------- AFTER -------------------------
        cd $fastq_screen_After
        ls > ids2.txt
        sed -i 's/\///g' ids2.txt

        for fname in $(< ids2.txt); do
            SRR=$(basename $fname "_clean_adapters.tagged.tagged_filter_screen.txt")
            awk -F "\t" -v org="$SRR" '{print org"\t"$1"\t"$2"\t"$3}' "$SRR"_clean_adapters.tagged.tagged_filter_screen.txt > "$Process_screen2/${SRR}_screen2.txt"
            sed -i '1d;1d;s/%//g' "$Process_screen2/${SRR}_screen2.txt"
        done

        cd $Process_screen2
        for screen_file in *_screen2.txt; do
            SRR=$(basename $screen_file "_screen2.txt")
            while IFS= read -r ids; do
                printf "$ids\t" >> "${SRR}_FQS_After.txt"
            done < "$screen_file"
        done

        sed -i -e '$a\\' *_FQS_After.txt
        cat *_FQS_After.txt > $fastq_screen_After/Fastqscreen_After.csv

        ###########################################################
        # STEP 2: Copy QUAST Results
        ###########################################################
        src_quast="$SRA_path/$ids/quast/transposed_report.tsv"
        dest_quast="$quast_result/${id}_transposed_report.tsv"

        if [[ -f "$src_quast" ]]; then
            cp "$src_quast" "$dest_quast"
            echo "Copied $src_quast to $dest_quast"
        else
            echo "$src_quast not found. Skipping."
        fi

        ###########################################################
        # STEP 3: Copy Binning Results
        ###########################################################
        src_bin="$SRA_path/$ids/Bin/${id}_binning_result/Selected_bin.csv"
        dest_bin="$Bin_result/${id}_Selected_bin.csv"

        if [[ -f "$src_bin" ]]; then
            cp "$src_bin" "$dest_bin"
            echo "Copied $src_bin to $dest_bin"
        else
            echo "$src_bin not found. Skipping."
        fi

        ###########################################################
        # STEP 4: Copy Taxonomic Classification Outputs
        ###########################################################
        cp "$SRA_path/$ids/bracken/${id}_merged_bracken_S.more1000" $bracken1
        cp "$SRA_path/$ids/bracken/${id}_merged_bracken_G.more1000" $bracken2
        cp "$SRA_path/$ids/bracken/${id}_merged_bracken_F.more1000" $bracken3

        ###########################################################
        # STEP 5: Copy Functional Annotations (CARD & VFDB)
        ###########################################################
        cp "$SRA_path/$ids/ARG/${id}_merged_final_output.tsv_ARG_updated.tsv" $CARD
        cp "$SRA_path/$ids/VFDB/${id}_merged_final_output.tsv_vfdb_updated" $VFDB

    done
done

###########################################################
# Optional: Merge all files (uncomment if needed)
###########################################################

# cat $bracken3/*.more1000 > $bracken3/merged_family
# cat $bracken2/*.more1000 > $bracken2/merged_genus
# cat $bracken1/*.more1000 > $bracken1/merged_species
# cat $quast_result/*.tsv > $quast_result/Quast_merged
# cat $Bin_result/*.csv > $Bin_result/Selected_bin_merged
# cat $CARD/*.tsv > $CARD/CARD_merged
# cat $VFDB/*_updated > $VFDB/VFDB_merged
# echo "Concatenation complete."

