#!/bin/bash
# Full-scale embeddings on a Princeton Adroit GPU node (optional: the Week 6 run used the local Apple GPU).
# Usage, from the repository root on Adroit after `make setup` and copying data/processed/mammalia/corpus.parquet:
#   sbatch scripts/slurm/embed_adroit.sh
# Outputs go to data/processed/mammalia/embeddings/, the same files scripts/06_embed.py writes locally.
#SBATCH --job-name=lc-embed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/embed_%j.out

set -euo pipefail
module purge
module load anaconda3/2024.6
source .venv/bin/activate
export HF_HOME="${HF_HOME:-/scratch/network/$USER/hf}"   # model weights cache off the home quota
python scripts/06_embed.py --clade mammalia
