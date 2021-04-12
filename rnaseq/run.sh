#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=HEPG2_rna_seq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nina.ripin@colorado.edu
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
--fasta /scratch/Shares/rinnclass/data/genomes/GRCh38.p13.genome.fa \
--gtf /scratch/Shares/rinnclass/data/genomes/gencode.v32.annotation.gtf \
--gencode \
--email nina.ripin@colorado.edu \
-c nextflow.config

date
