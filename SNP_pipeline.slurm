#!/bin/bash -l
#SBATCH --account=<account>
#SBATCH --clusters=<cluster>
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --job-name=SNPCall_Alzheimers

# Exporting environment.
# <export PATH: .....>

# Creating workflow.rst for the report's landing page.
echo "" > workflow.rst
echo "**SNPCall Pipeline for Alzheimer's**" >> workflow.rst
echo "" >> workflow.rst
echo "" >> workflow.rst
echo "This is a pipeline for SNP Variant calling in Alzheimer's patients looking specifically at genes: "APP", "SOD1" and "DYRK1A". " >> workflow.rst
echo "" >> workflow.rst
echo "All of the analysis will be performed on Chromosome 21." >> workflow.rst
echo "" >> workflow.rst
echo "Functioning DAG_ graph." >> workflow.rst

# Generating correct DAG plot (Snakemake's current version always displays a vertical unclear plot).
# This generates the old version's one to be inserted in the report.
snakemake --forceall --rulegraph | dot -Tpng > dag.png

# Running snakemake
snakemake -c1