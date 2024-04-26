#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH -p public2023
#SBATCH --mem=4000
#SBATCH --time=0-01:00:00
#SBATCH --mail-user=nfiege@uni-kassel.de
#SBATCH --mail-type=ALL

bench_type=${1}
base_dir=/home/groups/fb16-digi/projects/constmult/sat_mcm/
script_base=benchmark/linux-cluster
results_dir=${base_dir}benchmark/results/cmm_rnd

cd ${base_dir};python3 ${script_base}/cleanup_intermediate.py ${results_dir} ${bench_type}
