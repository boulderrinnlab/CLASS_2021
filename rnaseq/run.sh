#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=HEPG2_rna_seq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vipa5343@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=20:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "Here we go You've requested $SLURM_CPUS_ON_NODE core."

module load singularity/3.1.1

nextflow run nf-core/rnaseq -r 3.0 \
-resume \
-profile singularity \
--input design.csv \
--aligner star_salmon \
--fasta /Shares/rinn_class/data/genomes/human/gencode/v32/GRCh38.p13.genome.fa \
--gtf /Shares/rinn_class/data/genomes/human/gencode/v32/gencode.v32.annotation.gtf \
--gencode \
--email vipa5343@colorado.edu \
-c nextflow.config

date