#!/bin/bash

source /broad/software/scripts/useuse
use Anaconda3

#$ -l h_vmem=8G
#$ -l h_rt=04:00:00
#$ -l os=RedHat7

eval "$(conda shell.bash hook)"
conda create -y -n convert_bam python=3.6
conda activate /seq/epiprod02/kdong/environments/convert_bam

pip install --upgrade pip
pip install pysam
pip install numpy
pip install biopython
