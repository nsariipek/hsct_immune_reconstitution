# Nurefsan and Peter, 241227
# Run Cellranger multi for TP53 Immune Escape project

cd /broad/vangalenlab/sariipek/cellranger

# Now run the jobs
use UGER

# Test with one sample:
#sample_id="P1764_MIX"
#qsub -N $sample_id -o ${sample_id}.log cellranger.sh $sample_id

# File containing sample IDs (first column). Run each separately
SAMPLE_SHEET="1_SampleSheet_221214.txt"; SEQDATE=221214 # NOT DONE YET
SAMPLE_SHEET="1_SampleSheet_230302.txt"; SEQDATE=230302 # NOT DONE YET
SAMPLE_SHEET="1_SampleSheet_230317.txt"; SEQDATE=230317 # STARTED BY PETER
SAMPLE_SHEET="1_SampleSheet_241015.txt"; SEQDATE=241015 # STARTED BY PETER
SAMPLE_SHEET="1_SampleSheet_241217.txt"; SEQDATE=241217 # STARTED BY NUREFSAN

# Loop through each line in the file
while IFS=',' read -r sample_id _; do
  # Skip empty lines
  if [[ -z "$sample_id" ]]; then
    continue
  fi

  # Submit the job using qsub
  qsub -N $sample_id -o ${sample_id}.log -j y 3_Cellranger.sh $sample_id $SEQDATE
done < $SAMPLE_SHEET

# Monitor progress
watch "qstat; echo; tail *.log"

