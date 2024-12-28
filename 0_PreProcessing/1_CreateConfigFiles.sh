# Nurefsan and Peter, 241227
# Create sample sheets automatically

cd /broad/vangalenlab/sariipek/cellranger/

# Input sample sheet format. One for each NovaSeq run (five)
SAMPLE_SHEET="1_SampleSheet_221214.txt"; SEQDATE=221214 # NOT DONE YET
SAMPLE_SHEET="1_SampleSheet_230302.txt"; SEQDATE=230302 # NOT DONE YET
SAMPLE_SHEET="1_SampleSheet_230317.txt"; SEQDATE=230317 # STARTED BY PETER
SAMPLE_SHEET="1_SampleSheet_241015.txt"; SEQDATE=241015 # STARTED BY PETER
SAMPLE_SHEET="1_SampleSheet_241217.txt"; SEQDATE=241217 # STARTED BY NUREFSAN

# Make output folder
OUTDIR=multi_configs_${SEQDATE}
mkdir $OUTDIR

# Parse the sample sheet and group entries by sample ID (first make sure it's cleared)
unset SAMPLE_CONFIGS
declare -A SAMPLE_CONFIGS

# Read the sample sheet
while IFS=',' read -r sample_id fastq_id fastqs feature_types; do
  # Append this line to the correct sample's config data
  if [[ -n "${SAMPLE_CONFIGS[$sample_id]}" ]]; then
    SAMPLE_CONFIGS[$sample_id]+="\n$fastq_id,$fastqs,$feature_types,"
  else
    SAMPLE_CONFIGS[$sample_id]="$fastq_id,$fastqs,$feature_types,"
  fi
done < $SAMPLE_SHEET

# Generate config csv files for each sample
for sample_id in "${!SAMPLE_CONFIGS[@]}"; do
  SAMPLE_CONFIG=/broad/vangalenlab/sariipek/cellranger/$OUTDIR/${sample_id}_config.csv

# Write the config csv file
cat << EOF > "$SAMPLE_CONFIG"
[gene-expression]
reference,/broad/vangalenlab/vangalen/Genomes/GRCh38.221223/GRCh38
create-bam,true
[vdj]
reference,/broad/vangalenlab/sariipek/my_vdj_ref
[libraries]
fastq_id,fastqs,feature_types
$(echo -e "${SAMPLE_CONFIGS[$sample_id]}" | sed '/^$/d')
EOF

done

