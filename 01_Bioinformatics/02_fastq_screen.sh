#!/bin/bash

################################################################################
# SLURM Job Configuration
################################################################################

#SBATCH --nodes=1                     # Number of compute nodes requested
#SBATCH --ntasks=56                   # Total number of CPU tasks (cores)
#SBATCH --time=35:00:00               # Walltime (35 hours)
#SBATCH --partition=compute           # Queue partition
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/out/slurm-%j.out2
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Environmental_Metagenomes/err/slurm-%j.err2

# Initialize Conda shell
eval "$(conda shell.bash hook)"

################################################################################
# INPUT PARAMETERS (Passed from run.sh or user)
################################################################################

ids="$1"                             # Directory identifier (SRA sample)
id="$2"                              # Unique sample ID (individual file name)
SRA_path="$3"                        # Base path to downloaded SRA files
bbmap_path="$4"                      # Path to BBMap tools
sh_path="$5"                         # Path to helper scripts
database="$6"                        # Kraken/Bracken database path
protein_nr="$7"                      # Protein annotation DB
collect="$8"                         # Output collection directory
diamond_path="$9"                    # DIAMOND binary path

# Fixed reference databases
swiss_prot_bacteria=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Shotgun-metagenomics/Swiss-Prot_bacteria
UniProtKB_all_bacteria_proteines=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/UniProtKB_all_bacteria_proteines
VFDB_diamond_database=/srv/lustre01/project/mmrd-cp3fk69sfrq/shared/VFDB
ARG_diamond_database=/srv/lustre01/project/mmrd-cp3fk69sfrq/shared/ARGS_CARD

# Navigate to sample directory
cd $SRA_path/$ids

################################################################################
# STEP 1: Create All Necessary Output Directories
################################################################################

mkdir -p $SRA_path/$ids/{prokka,bamfile,Bin,Final_bin,diamond,kraken,bracken,VFDB,ARG,quast}

prokka=$SRA_path/$ids/prokka
bamfile=$SRA_path/$ids/bamfile
Bin=$SRA_path/$ids/Bin
Final_bin=$SRA_path/$ids/Final_bin
diamond=$SRA_path/$ids/diamond
kraken=$SRA_path/$ids/kraken
bracken=$SRA_path/$ids/bracken
VFDB=$SRA_path/$ids/VFDB
ARG=$SRA_path/$ids/ARG
quast=$SRA_path/$ids/quast

################################################################################
# STEP 2: Contaminant Screening using FastQ Screen
################################################################################

conda activate fastq_screen

# Run fastq_screen for each of the cleaned and adapter-trimmed read sets
for read_type in "" "_unmerged1" "_unmerged2"; do
    fastq_screen --tag "$id${read_type}_clean_adapters.fq"
    fastq_screen --nohits "$id${read_type}_clean_adapters.tagged.fastq"
    fastq_screen "$id${read_type}_clean_adapters.tagged.tagged_filter.fastq"
done

################################################################################
# STEP 3: Post-Filtering Quality Control using FastQC
################################################################################

conda activate fastqc

# Perform FastQC on the filtered reads
for read_type in "" "_unmerged1" "_unmerged2"; do
    fastqc "$id${read_type}_clean_adapters.tagged.tagged_filter.fastq" -f fastq -O .
done

################################################################################
# STEP 4: Metagenome Assembly with SPAdes (metaSPAdes)
################################################################################

conda activate spades

spades.py --meta \
    --merged "$id"_clean_adapters.tagged.tagged_filter.fastq \
    -s "$id"_unmerged1_clean_adapters.tagged.tagged_filter.fastq \
    -s "$id"_unmerged2_clean_adapters.tagged.tagged_filter.fastq \
    -t 56 \
    -o metaspades

################################################################################
# STEP 5: Header Simplification + Length Filtering + QUAST
################################################################################

conda activate prokka
perl $sh_path/simplifyFastaHeaders.pl \
    $SRA_path/$ids/metaspades/contigs.fasta \
    xx \
    $prokka/"$id"_contigs_modified_headers.fasta \
    $Bin/"$id"_contigs_header.map

# Filter out contigs less than 1,000 bp
conda activate seqtk
seqtk seq -L 1000 $prokka/"$id"_contigs_modified_headers.fasta > $prokka/"$id"_contigs_modified_headers.more1000.fasta
printf "\n Contigs <1000bp removed.\n"

# Assess contig quality
conda activate quast
python3 /home/user/anaconda3/envs/quast/bin/quast \
    $prokka/"$id"_contigs_modified_headers.more1000.fasta \
    --threads 56 \
    -o $quast
printf "\n QUAST completed.\n"

################################################################################
# STEP 6: Read Mapping and Binning Preparation (MetaBAT2)
################################################################################

# Generate BAM file using BBMap
bash $bbmap_path/bbmap.sh \
    ref=$prokka/"$id"_contigs_modified_headers.fasta \
    in="$id"_1.fastq \
    in2="$id"_2.fastq \
    out=$bamfile/alignment.bam
printf "\n BAM file created for $ids.\n"

# Sort and index BAM file
samtools sort -o $bamfile/alignment_sorted.bam $bamfile/alignment.bam
samtools index $bamfile/alignment_sorted.bam

# Clean up intermediate BAM files
rm $bamfile/alignment.bam
rm $bamfile/alignment_sorted.bam.bai


###########################################
# STEP 7: Metagenomic Binning with MetaBAT2
###########################################

conda activate metabat2

# Summarize read depth per contig for MetaBAT2
jgi_summarize_bam_contig_depths \
    --outputDepth $bamfile/depth.txt \
    $bamfile/alignment_sorted.bam

# Run MetaBAT2 for binning (MAG recovery)
metabat2 \
    -i $prokka/"$id"_contigs_modified_headers.more1000.fasta \
    -a $bamfile/depth.txt \
    -o $prokka/"$id"_binning_result/bin \
    -t 56 \
    -m 1500 --maxP 90 --minS 80 --maxEdges 500 --minCV 1

# Rename binned files by prefixing sample ID for traceability
for fa_file in "$prokka/${id}_binning_result"/*.fa; do
    base_name=$(basename "$fa_file")
    mv "$fa_file" "$prokka/${id}_binning_result/${id}_${base_name}"
done

###########################################
# STEP 8: Bin Quality Assessment with CheckM
###########################################

conda activate comebin_env

# Run CheckM to evaluate bin completeness and contamination
checkm lineage_wf -t 56 -x fa \
    $prokka/"$id"_binning_result/ \
    $prokka/"$id"_binning_result/

# Summarize CheckM output to TSV format
checkm qa \
    $prokka/"$id"_binning_result/lineage.ms \
    $prokka/"$id"_binning_result \
    > $prokka/"$id"_binning_result/checkm_quality_report.tsv

# Clean report: remove header and dashed lines
sed '1,6d; /^-----------------------------------------------------------------------------------------------------------------------------------------------------------------------/d' \
    $prokka/"$id"_binning_result/checkm_quality_report.tsv \
    > $prokka/"$id"_binning_result/cleaned_checkm_report.tsv

# Filter bins: ≥90% completeness AND ≤10% contamination
input_file="$prokka/$id"_binning_result/cleaned_checkm_report.tsv
output_file="$prokka/$id"_binning_result/Selected_bin.csv

{
    read
    while read -r line; do
        completeness=$(echo "$line" | awk '{print $(NF-2)}')
        contamination=$(echo "$line" | awk '{print $(NF-1)}')

        if (( $(echo "$completeness >= 90" | bc -l) )) && (( $(echo "$contamination <= 10" | bc -l) )); then
            echo "$line"
        fi
    done
} < "$input_file" > "$output_file"

echo "Filtered bins saved to $output_file"

# Backup bins and contig files for next steps
cp $prokka/"$id"_contigs_modified_headers.more1000.fasta $Bin
cp -R $prokka/"$id"_binning_result $Bin

###########################################
# STEP 9: Clean Intermediate Files
###########################################

rm -rf $SRA_path/$ids/prokka
rm -rf $SRA_path/$ids/metaspades
rm -rf $SRA_path/$ids/ref
rm -rf $SRA_path/$ids/bamfile

###########################################
# STEP 10: Move Selected High-Quality Bins
###########################################

selected_bin_csv="$Bin/"$id"_binning_result/Selected_bin.csv"

if [[ ! -f "$selected_bin_csv" ]]; then
    echo "Selected_bin.csv not found!"
    exit 1
fi

# Move selected .fa bins to Final_bin
while IFS=, read -r id_col _; do
    id_trimmed=$(echo "$id_col" | awk '{print $1}')
    fa_file=$Bin/"$id"_binning_result/"${id_trimmed}.fa"

    if [[ -f "$fa_file" ]]; then
        cp "$fa_file" "$Final_bin"
        echo "Moved $fa_file to $Final_bin"
    else
        echo "File $fa_file not found!"
    fi
done < "$selected_bin_csv"

echo "All selected bins have been moved to Final_bin."

###########################################
# STEP 11: Taxonomic Assignment with Kraken2 + Bracken
###########################################

conda activate kraken2

for fa_file in "$Final_bin"/*.fa; do
    id=$(basename "$fa_file" .fa)

    # Kraken2 for classification
    kraken2 --db "$database" --threads 56 "$fa_file" \
        --classified-out "$kraken/$id.classified-out.more1000" \
        --unclassified-out "$kraken/$id.unclassified-out.more1000" \
        --use-names \
        --report "$kraken/$id.report.more1000" \
        --use-mpa-style \
        --report-zero-counts > "$kraken/$id.screen.more1000"

    kraken2 --db "$database" --threads 56 "$fa_file" \
        --use-names \
        --report "$kraken/$id.report2.more1000" > "$kraken/$id.screen2.more1000"

    awk -F "\t" '{print $2"\t"$3}' "$kraken/$id.screen.more1000" > "$kraken/$id.contigs_species.names"

    echo "Kraken2 completed for $fa_file"

    ###########################################
    # Bracken for quantitative classification
    ###########################################

    conda activate bracken

    for level in S F G; do
        bracken -d "$database" \
            -i "$kraken/$id.report2.more1000" \
            -o "$bracken/$id.bracken_${level}.more1000" \
            -l $level

        # Prepend sample ID to results
        if [[ -f "$bracken/$id.bracken_${level}.more1000" ]]; then
            awk -v id="$id" '{print id "\t" $0}' "$bracken/$id.bracken_${level}.more1000" > "${bracken}/${id}.bracken_${level}.temp"
            mv "${bracken}/${id}.bracken_${level}.temp" "$bracken/$id.bracken_${level}.more1000"
        else
            echo "Missing Bracken output for level $level"
        fi
    done

    echo "Bracken processing done for $id"
done

###########################################
# STEP 12: Merge Bracken Outputs
###########################################

mapfile -t ids < ids2.txt

for id in "${ids[@]}"; do
    id_name=$(basename "$id" ".sra")
    for level in S F G; do
        merged_output_file="$bracken/${id_name}_merged_bracken_${level}.more1000"
        awk 'NR==1 || FNR>1 {print "'"$id_name"'" "\t" $0}' "$bracken"/*.bracken_${level}.more1000 > "$merged_output_file"
        echo "Merged Bracken files at level $level for $id_name"
    done
done


###############################################
# STEP 13: VFDB and ARG Functional Annotation
###############################################

# Loop over each .fa file in the Final_bin folder (binned genomes)
for fa_file in "$Final_bin"/*.fa; do
    id=$(basename "$fa_file" .fa)  # Extract sample ID

    #####################
    # VFDB Annotation
    #####################
    # Run DIAMOND BLASTX against VFDB for virulence gene detection
    $diamond_path blastx --threads 56 --sensitive --header --max-target-seqs 1 \
        --evalue 1e-05 --outfmt 6 -d "$VFDB_diamond_database/VFDB_diamond_database" \
        -q "$fa_file" -o "$VFDB/${id}_VFDB_output.csv"

    sed -i 's/|/\t/g' "$VFDB/${id}_VFDB_output.csv"  # Optional: make output tab-delimited
    grep ">" "$VFDB_diamond_database/VFDB_setB_pro.fas" > "$VFDB/${id}_VFDB_output_file2"
    sed 's/>//g' "$VFDB/${id}_VFDB_output_file2" > "$VFDB/${id}_ids_vfdb"
    echo "VFDB annotation completed for $id"

    #####################
    # ARG Annotation
    #####################
     $diamond_path blastx --threads 56 --sensitive --header --max-target-seqs 1 \
        --evalue 1e-05 --outfmt 6 -d "$ARG_diamond_database/ARG_diamond_database" \
        -q "$fa_file" -o "$ARG/${id}_ARG_output.csv"

    sed -i 's/|/\t/g' "$ARG/${id}_ARG_output.csv"
    awk -F "\t" '{print $3}' "$ARG/${id}_ARG_output.csv" > "$ARG/${id}_ARG_diamond_swiss_prot_bacteria.protein.csv"
    grep -f "$ARG/${id}_ARG_diamond_swiss_prot_bacteria.protein.csv" "$ARG_diamond_database/aro_index.tsv" > "$ARG/${id}_ARG_diamond.grep"
    echo "ARG annotation completed for $id"
done

###############################################
# STEP 14: Processing VFDB OUTPUT FILES
###############################################

# Loop over each .csv file in the VFDB folder
for vfdb_file in "$VFDB"/*_VFDB_output.csv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$vfdb_file" _VFDB_output.csv)

    # Check if the VFDB file exists
    if [[ ! -f "$vfdb_file" ]]; then
        echo "Error: VFDB file $vfdb_file does not exist."
        continue
    fi

    # Remove the first two lines (the header and description) and save cleaned output
    sed '1,2d' "$vfdb_file" > "$VFDB/${id}_VFDB_output_cleaned.csv"

    # Extract the accession numbers (e.g., VFG013387) from cleaned output
    #awk -F "\t" '{split($2,a,"("); print a[1]}' "$VFDB/${id}_VFDB_output_cleaned.csv" > "$VFDB/${id}_accession_list_vfdb"
    awk -F "\t" '{split($2, a, "\\("); print a[1]}' "$VFDB/${id}_VFDB_output_cleaned.csv" > "$VFDB/${id}_accession_list_vfdb"

    echo "Processed accession list saved for $id."

    # Check if the accession list exists
    if [[ ! -f "$VFDB/${id}_accession_list_vfdb" ]]; then
        echo "Error: Accession list ${id}_accession_list_vfdb does not exist for $id."
        continue
    fi

    # Check if ids_vfdb exists
    if [[ ! -f "$VFDB/${id}_ids_vfdb" ]]; then
        echo "Error: Required file ids_vfdb does not exist."
        continue
    fi

    # Create or empty the matched_vfdb file
    matched_vfdb_file="$VFDB/${id}_matched_vfdb"
    > "$matched_vfdb_file"

    # Read each ID from the accession list and match with ids_vfdb
    while IFS= read -r accession_id; do
        accession_id=$(echo "$accession_id" | xargs)  # Trim any leading/trailing whitespace
        if [[ -z "$accession_id" ]]; then
            continue
        fi
        # Append all lines from ids_vfdb that contain the ID
       grep -F "$accession_id" "$VFDB/${id}_ids_vfdb" >> "$matched_vfdb_file"
    done < "$VFDB/${id}_accession_list_vfdb"

    echo "Matching lines have been written to $matched_vfdb_file."
done

   # Loop over each matched vfdb file in the VFDB folder
for matched_vfdb_file in "$VFDB"/*_matched_vfdb; do
    # Extract the base name to determine the ID
    id=$(basename "$matched_vfdb_file" _matched_vfdb)

    # Define the output file path for processed results
    output_file="$VFDB/${id}_processed_matched_vfdb.tsv"  # Output file path
    > "$output_file"  # Create or empty the output file

    # Check if matched_vfdb_file is empty
    if [[ ! -s "$matched_vfdb_file" ]]; then
        echo "Warning: $matched_vfdb_file is empty. Skipping processing."
        continue
    fi

    # Process the matched vfdb file
    while IFS= read -r line; do
        # Split by parentheses
        IFS='()' read -r vfg accession description <<< "$line"

        # Trim whitespace
        vfg=$(echo "$vfg" | xargs)
        accession=$(echo "$accession" | xargs)
        description=$(echo "$description" | xargs)

        # Write to the output file, separating by tabs
        echo -e "$vfg\t$accession\t$description" >> "$output_file"
    done < "$matched_vfdb_file"

    echo "Processed matched VFDB file saved to $output_file."
done

# Loop over each .csv file in the VFDB folder
for vfdb_file in "$VFDB"/*_VFDB_output.csv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$vfdb_file" _VFDB_output.csv)

    # Define the input and output files
    processed_matched_file="$VFDB/${id}_processed_matched_vfdb.tsv"
    cleaned_vfdb_file="$VFDB/${id}_VFDB_output_cleaned.csv"
    merged_output_file="$VFDB/${id}_merged_output.tsv"

    # Check if necessary files exist
    if [[ ! -f "$processed_matched_file" ]]; then
        echo "Error: Processed matched file $processed_matched_file does not exist."
        continue
    fi

    if [[ ! -f "$cleaned_vfdb_file" ]]; then
        echo "Error: Cleaned VFDB file $cleaned_vfdb_file does not exist."
        continue
    fi

    # Create or empty the merged output file
    > "$merged_output_file"

    # Read the processed matched file and create an associative array for VFG IDs
    declare -A vfg_map
    while IFS=$'\t' read -r vfg_id accession_id description; do
        # Store each line in the associative array using VFG ID as the key
        vfg_map["$vfg_id"]="$accession_id\t$description"
    done < "$processed_matched_file"

    # Read the cleaned VFDB file and perform the merge
    while IFS= read -r line; do
        # Extract the VFG ID from the line (second column)
        vfg_id=$(echo "$line" | awk '{print $2}' | sed 's/(.*//;s/gb//')  # Clean to get just the VFG ID

        # Check if the VFG ID exists in the associative array
        if [[ -n "${vfg_map[$vfg_id]}" ]]; then
            # Merge the lines
            echo -e "$line\t${vfg_map[$vfg_id]}" >> "$merged_output_file"
        fi
    done < "$cleaned_vfdb_file"

    echo "Merged output saved to $merged_output_file."
done

# Loop over each .csv file in the VFDB folder
for vfdb_file in "$VFDB"/*_VFDB_output.csv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$vfdb_file" _VFDB_output.csv)

    # Define the input files and the output file for the final merge
    merged_output_file="$VFDB/${id}_merged_output.tsv"
    kraken_file="$kraken/$id.contigs_species.names"
    final_output_file="$VFDB/${id}_final_output_with_species.tsv"

    # Check if necessary files exist
    if [[ ! -f "$merged_output_file" ]]; then
        echo "Error: Merged output file $merged_output_file does not exist."
        continue
    fi

    if [[ ! -f "$kraken_file" ]]; then
        echo "Error: Kraken species names file $kraken_file does not exist."
        continue
    fi

    # Create or empty the final output file
    > "$final_output_file"

    # Read the Kraken species file and create an associative array for contig IDs
    declare -A contig_species_map
    while IFS=$'\t' read -r contig_id species_info; do
        # Store each contig ID and species info in the associative array
        contig_species_map["$contig_id"]="$species_info"
    done < "$kraken_file"

    # Read the merged VFDB output file and perform the final merge
    while IFS= read -r line; do
        # Extract the contig ID (first column)
        contig_id=$(echo "$line" | awk '{print $1}')

        # Check if the contig ID exists in the associative array
        species="${contig_species_map[$contig_id]}"

        # Append the species information as a separate column
        if [[ -n "$species" ]]; then
            echo -e "$line\t$species" >> "$final_output_file"
        else
            echo -e "$line\tN/A" >> "$final_output_file"
        fi
    done < "$merged_output_file"

    echo "Final output with species saved to $final_output_file."
done


# Loop over each final output file in the VFDB folder
for final_output_file in "$VFDB"/*_final_output_with_species.tsv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$final_output_file" _final_output_with_species.tsv)

    # Check if the file exists
    if [[ ! -f "$final_output_file" ]]; then
        echo "Error: Final output file $final_output_file does not exist."
        continue
    fi

    # Create a new file with id prepended
    final_output_with_id="$VFDB/${id}_final_output_with_species_and_id.tsv"
    > "$final_output_with_id"  # Create or empty the new output file

    # Prepend the id to each line and save to the new file
    while IFS= read -r line; do
        echo -e "$id\t$line" >> "$final_output_with_id"
    done < "$final_output_file"

    echo "Final output with ID prepended saved to $final_output_with_id."
done

# Loop through each final output file
for final_output_with_id in "$VFDB"/*_final_output_with_species_and_id.tsv; do
    # Extract the base name (ID)
    id=$(basename "$final_output_with_id" _final_output_with_species_and_id.tsv)

    # Define the output file for the modified content
    output_file="$VFDB/${id}_final_output_with_species_and_id_separated.tsv"

    # Use awk to separate the 16th column
    awk -F"\t" '{
        # Capture the original 16th column
        original_column = $16;
        
        # Use regex to match the pattern
        if (match(original_column, /\(([^)]+)\)/, m)) {
             m[1] contains the text inside parentheses
            print $0 "\t" m[1];  # Print the original line and the new column
        } else {
            print $0 "\t";  # Print the original line without adding a new column
        }
    }' OFS="\t" "$final_output_with_id" > "$output_file"
done

# Read IDs from ids2.txt into an array
mapfile -t ids < ids2.txt

# Merge all separated files into a single temporary file
cat "$VFDB"/*_final_output_with_species_and_id_separated.tsv > "$VFDB/temporary_merged_output.tsv"

# Loop through the IDs and rename the merged output
for id in "${ids[@]}"; do
    # Remove the .sra extension to get the basename
    id_name=$(basename "$id" ".sra")
    
    # Rename the merged output file with the current ID
    mv "$VFDB/temporary_merged_output.tsv" "$VFDB/${id_name}_merged_final_output.tsv"
done
# Loop over each final output file in the VFDB folder
for final_output_with_id in "$VFDB"/*_merged_final_output.tsv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$final_output_with_id" _merged_final_output.tsv)

    # Check if the file exists
    if [[ ! -f "$final_output_with_id" ]]; then
        echo "Error: Final output file $final_output_with_id_separated does not exist."
        continue
    fi

    # Create a new file for storing the VF hit counts
    vf_hit_counts_file="$VFDB/${id}_vf_hit_counts.tsv"
    > "$vf_hit_counts_file"  # Create or empty the new output file

    # Count the number of hits for each VF based on the bins
    awk -F"\t" '{ print $18 }' "$final_output_with_id" | sort | uniq -c | while read count vf; do
        echo -e "$vf\t$count" >> "$vf_hit_counts_file"
    done

    echo "VF hit counts saved to $vf_hit_counts_file."
done

# Loop over each final output file in the VFDB folder
for final_output_with_id in "$VFDB"/*_merged_final_output.tsv; do
     #Extract the base name of the file (prefix ID)
    id=$(basename "$final_output_with_id" _merged_final_output.tsv)

    # Define the vf_hit_counts file
    vf_hit_counts_file="$VFDB/${id}_vf_hit_counts.tsv"

    # Check if the vf_hit_counts file exists
    if [[ ! -f "$vf_hit_counts_file" ]]; then
        echo "Error: VF hit counts file $vf_hit_counts_file does not exist."
       continue
    fi

    # Create the output file with the updated data
    updated_output_file="${final_output_with_id}_vfdb_updated"
    > "$updated_output_file"  # Create or empty the new output file

    # Append the second column from vf_hit_counts_file based on matching last column of final_output_with_id
    while IFS= read -r line; do
        # Get the last column from the final_output_with_id
        last_column=$(echo "$line" | awk '{print $NF}')
        
        # Get the corresponding count from vf_hit_counts_file
        count=$(awk -v vf="$last_column" '$1 == vf {print $2}' "$vf_hit_counts_file")

        # If a count is found, append it to the line along with the ID in the first column
        if [[ -n "$count" ]]; then
            echo -e "${id}\t${line}\t${count}" >> "$updated_output_file"
        else
           echo -e "${id}\t${line}\t0" >> "$updated_output_file"  # Append 0 if no match is found
        fi
    done < "$final_output_with_id"

    echo "Appended second column from $vf_hit_counts_file based on last column to $final_output_with_id, including the ID as the first column."
done

###############################################
# STEP 15: Processing ARG OUTPUT FILES
###############################################

# Loop over each .csv file in the ARG folder
for ARG_file in "$ARG"/*_ARG_output.csv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$ARG_file" _ARG_output.csv)

    # Check if the VFDB file exists
    if [[ ! -f "$ARG_file" ]]; then
        echo "Error: ARG file $ARG_file does not exist."
        continue
    fi

    # Remove the first two lines (the header and description) and save cleaned output
    sed '1,2d' "$ARG_file" > "$ARG/${id}_ARG_output_cleaned.csv"

    # Extract the accession numbers from cleaned output
     awk -F "\t" '{print $4}' "$ARG/${id}_ARG_output_cleaned.csv" > "$ARG/${id}_accession_list_ARG"


    echo "Processed accession list saved for $id."

    # Check if the accession list exists
    if [[ ! -f "$ARG/${id}_accession_list_ARG" ]]; then
        echo "Error: Accession list ${id}_accession_list_ARG does not exist for $id."
        continue
    fi

    # Check if {id}_ARG_diamond.grep exists
    if [[ ! -f "$ARG/${id}_ARG_diamond.grep" ]]; then
        echo "Error: Required file {id}_ARG_diamond.grep does not exist."
        continue
    fi

    # Create or empty the matched_ARG file
    matched_ARG_file="$ARG/${id}_matched_ARG"
    > "$matched_ARG_file"

    # Read each ID from the accession list and match with {id}_ARG_diamond.grep
    while IFS= read -r accession_id; do
        accession_id=$(echo "$accession_id" | xargs)  # Trim any leading/trailing whitespace
        if [[ -z "$accession_id" ]]; then
            continue
        fi
        # Append all lines from {id}_ARG_diamond.grep that contain the ID
       grep -F "$accession_id" "$ARG/${id}_ARG_diamond.grep" >> "$matched_ARG_file"
    done < "$ARG/${id}_accession_list_ARG"

    echo "Matching lines have been written to $matched_ARG_file."
done


# Loop over each .csv file in the VFDB folder
for ARG_file in "$ARG"/*_ARG_output.csv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$ARG_file" _ARG_output.csv)

    # Define the input and output files
    processed_matched_file="$ARG/${id}_matched_ARG"
    cleaned_ARG_file="$ARG/${id}_ARG_output_cleaned.csv"
    merged_output_file="$ARG/${id}_merged_output.tsv"

    # Check if necessary files exist
    if [[ ! -f "$processed_matched_file" ]]; then
        echo "Error: Processed matched file $processed_matched_file does not exist."
        continue
    fi

    if [[ ! -f "$cleaned_ARG_file" ]]; then
        echo "Error: Cleaned ARG file $cleaned_ARG_file does not exist."
        continue
    fi

    # Create or empty the merged output file
    > "$merged_output_file"

    # Read the processed matched file and create an associative array for ARO IDs
    declare -A ARO_map
    while IFS=$'\t' read -r ARO_id accession_id description; do
        # Store each line in the associative array using ARO ID as the key
        ARO_map["$ARO_id"]="$accession_id\t$description"
    done < "$processed_matched_file"

    # Read the cleaned ARG file and perform the merge
    while IFS= read -r line; do
        # Extract the ARO ID from the line (FOURTH column)
        ARO_id=$(echo "$line" | awk '{print $4}')   # Clean to get just the ARO ID



        # Check if the ARO ID exists in the associative array
        if [[ -n "${ARO_map[$ARO_id]}" ]]; then
            # Merge the lines
            echo -e "$line\t${ARO_map[$ARO_id]}" >> "$merged_output_file"
        fi
    done < "$cleaned_ARG_file"

    echo "Merged output saved to $merged_output_file."
done


# Loop over each .csv file in the ARG folder
for ARG_file in "$ARG"/*_ARG_output.csv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$ARG_file" _ARG_output.csv)

    # Define the input files and the output file for the final merge
    merged_output_file="$ARG/${id}_merged_output.tsv"
    kraken_file="$kraken/$id.contigs_species.names"
    final_output_file="$ARG/${id}_final_output_with_species.tsv"

    # Check if necessary files exist
    if [[ ! -f "$merged_output_file" ]]; then
        echo "Error: Merged output file $merged_output_file does not exist."
        continue
    fi

    if [[ ! -f "$kraken_file" ]]; then
        echo "Error: Kraken species names file $kraken_file does not exist."
        continue
    fi

    # Create or empty the final output file
    > "$final_output_file"

    # Read the Kraken species file and create an associative array for contig IDs
    declare -A contig_species_map
    while IFS=$'\t' read -r contig_id species_info; do
        # Store each contig ID and species info in the associative array
        contig_species_map["$contig_id"]="$species_info"
    done < "$kraken_file"

    # Read the merged ARG output file and perform the final merge
    while IFS= read -r line; do
        # Extract the contig ID (first column)
        contig_id=$(echo "$line" | awk '{print $1}')

        # Check if the contig ID exists in the associative array
        species="${contig_species_map[$contig_id]}"

        # Append the species information as a separate column
        if [[ -n "$species" ]]; then
            echo -e "$line\t$species" >> "$final_output_file"
        else
            echo -e "$line\tN/A" >> "$final_output_file"
        fi
    done < "$merged_output_file"

    echo "Final output with species saved to $final_output_file."
done

# Loop over each final output file in the ARG folder
for final_output_file in "$ARG"/*_final_output_with_species.tsv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$final_output_file" _final_output_with_species.tsv)

    # Check if the file exists
    if [[ ! -f "$final_output_file" ]]; then
        echo "Error: Final output file $final_output_file does not exist."
        continue
    fi

    # Create a new file with id prepended
    final_output_with_id="$ARG/${id}_final_output_with_species_and_id.tsv"
    > "$final_output_with_id"  # Create or empty the new output file

    # Prepend the id to each line and save to the new file
    while IFS= read -r line; do
        echo -e "$id\t$line" >> "$final_output_with_id"
    done < "$final_output_file"

    echo "Final output with ID prepended saved to $final_output_with_id."
done

# Read IDs from ids2.txt into an array
mapfile -t ids < ids2.txt

# Merge all separated files into a single temporary file
cat "$ARG"/*_final_output_with_species_and_id.tsv > "$ARG/temporary_merged_output.tsv"

# Loop through the IDs and rename the merged output
for id in "${ids[@]}"; do
    # Remove the .sra extension to get the basename
    id_name=$(basename "$id" ".sra")
    
    # Rename the merged output file with the current ID
    mv "$ARG/temporary_merged_output.tsv" "$ARG/${id_name}_merged_final_output.tsv"
done

# Loop over each final output file in the ARG folder
for final_output_with_id in "$ARG"/*_merged_final_output.tsv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$final_output_with_id" _merged_final_output.tsv)

    # Check if the file exists
    if [[ ! -f "$final_output_with_id" ]]; then
        echo "Error: Final output file $final_output_with_id_separated does not exist."
        continue
    fi

    # Create a new file for storing the ARG hit counts
    ARG_hit_counts_file="$ARG/${id}_ARG_hit_counts.tsv"
    > "$ARG_hit_counts_file"  # Create or empty the new output file

    # Count the number of hits for each ARG based on the bins
    awk -F"\t" '{ print $6 }' "$final_output_with_id" | sort | uniq -c | while read count ARG; do
        echo -e "$ARG\t$count" >> "$ARG_hit_counts_file"
    done

    echo "ARG hit counts saved to $ARG_hit_counts_file."
done


# Loop over each final output file in the ARG folder
for final_output_with_id in "$ARG"/*_merged_final_output.tsv; do
    # Extract the base name of the file (prefix ID)
    id=$(basename "$final_output_with_id" _merged_final_output.tsv)

    # Define the ARG_hit_counts file
    ARG_hit_counts_file="$ARG/${id}_ARG_hit_counts.tsv"

    # Check if the ARG_hit_counts file exists
    if [[ ! -f "$ARG_hit_counts_file" ]]; then
        echo "Error: ARG hit counts file $ARG_hit_counts_file does not exist."
        continue
    fi

    # Create the output file with the updated data
    updated_output_file="${final_output_with_id}_ARG_updated.tsv"  # Use .tsv extension for clarity
    > "$updated_output_file"  # Create or empty the new output file

    # Append the second column from ARG_hit_counts_file based on matching sixth column of final_output_with_id
    while IFS= read -r line; do
        # Get the sixth column from the final_output_with_id
        sixth_column=$(echo "$line" | awk '{print $6}')  # Change to get the sixth column

        # Get the corresponding count from ARG_hit_counts_file
        count=$(awk -v ARG="$sixth_column" '$1 == ARG {print $2}' "$ARG_hit_counts_file")

        # If a count is found, append it to the line along with the ID in the first column
        if [[ -n "$count" ]]; then
            echo -e "${id}\t${line}\t${count}" >> "$updated_output_file"
        else
            echo -e "${id}\t${line}\t0" >> "$updated_output_file"  # Append 0 if no match is found
        fi
    done < "$final_output_with_id"

    echo "Appended second column from $ARG_hit_counts_file based on sixth column to $final_output_with_id, including the ID as the first column."
done

#################################End of Processing Step for ARG #####################################################################################################################################

###############################################
#done
#done









