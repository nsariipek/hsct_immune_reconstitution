# Peter van Galen, 250209
# Quantify T cell programs using starCAT (https://www.biorxiv.org/content/10.1101/2024.05.03.592310v1, https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vignette_R.ipynb)  
# This should be run after creating matrix.mtx, barcodes.tsv, and features.tsv as described in 3.2_PvG-TCAT.R

# Installation (only run once)
#sudo apt install python3-pip
#pip install starcatpy
# Add to ~/.bashrc
#export PATH=$HOME/.local/bin:$PATH

# Determine program usage
starcat --reference "TCAT.V1" \
        --counts "/home/unix/vangalen/TP53_ImmuneEscape/3_DGE/starcat/matrix.mtx.gz" \
        --output-dir "/home/unix/vangalen/TP53_ImmuneEscape/3_DGE/starcat/" \
        --name "results"

# Delete cache folder that was generated