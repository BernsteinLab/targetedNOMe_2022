#!/bin/bash

source /broad/software/scripts/useuse
use Anaconda3

#$ -l h_vmem=8G
#$ -l h_rt=04:00:00
#$ -l os=RedHat7

eval "$(conda shell.bash hook)"

conda create -y -n whatshap python=3.6
source activate whatshap
conda install -c bioconda whatshap nomkl