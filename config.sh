#!/bin/bash

CURRENT_DIR=$(pwd)
PATH=$PATH:$CURRENT_DIR/.bin
export PATH

source ~/anaconda3/bin/activate
conda activate ./env_cpa
