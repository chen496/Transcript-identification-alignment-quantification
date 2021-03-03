#!/bin/bash
source ~/.bashrc
module remove python/2.7.3
module switch gcc/4.4.7 gcc/4.9.4
conda activate salmon
