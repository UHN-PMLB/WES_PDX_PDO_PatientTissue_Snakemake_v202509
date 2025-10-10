#! /bin/bash
#SBATCH -J WES-scheduler
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p all
#SBATCH --mem=2gb

snakemake \
--jobs 20 \
--profile slurm \
--cluster-config slurm/cluster.json \
--latency-wait 60 \
--rerun-incomplete
