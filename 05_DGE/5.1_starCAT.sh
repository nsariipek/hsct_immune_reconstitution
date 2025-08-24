# Peter van Galen, 250209
# Quantify T cell programs using starCAT (https://www.biorxiv.org/content/10.1101/2024.05.03.592310v1, https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vignette_R.ipynb)  
# This should be run after creating matrix.mtx, barcodes.tsv, and features.tsv as described in 3.1_PvG-TCAT.R
# This script was run on a Google Cloud Platform virtual machine with Ubuntu 22.04.5 LTS

# Installation (only run once)
#sudo apt install python3-pip
#pip install starcatpy
# Add the following to ~/.bashrc
#export PATH=$HOME/.local/bin:$PATH

# Determine program usage
starcat --reference "TCAT.V1" \
        --counts "/home/unix/vangalen/hsct_immune_reconstitution/AuxiliaryFiles/starCAT/matrix.mtx.gz" \
        --output-dir "/home/unix/vangalen/hsct_immune_reconstitution/05_DGE/5.1_starCAT/"

# Delete cache folder that was generated