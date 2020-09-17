#!/bin/bash

source ~/admd-homedir.bashrc

cd /home/osz/cgn-amber99-analysis

jupyter-lab > jupyter-lab.log &
JLABPROC="$!"

ssh -L 8888:localhost:8888 or-dgx-login01 -N -f 

wait "$JLABPROC"
