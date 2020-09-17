#!/bin/bash


SIMG=/containers/PE-Dask/rapidsai_0.8-cuda10.0-devel-centos7-py37.sif

srun --mem=150G --gpus=11 --cpus-per-gpu=4 -p dgx2a-birthright \
  singularity exec --nv --cleanenv --pwd /home/osz/cgn-amber99-analysis \
  ${SIMG} bash -c \
'
source ~/admd-homedir.bashrc
echo $(which python)
./02-calculate-pipeline-msm.py
'
