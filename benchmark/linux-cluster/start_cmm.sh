#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH -p public2023
#SBATCH --mem=251000
#SBATCH --time=0-00:01:00
#SBATCH --mail-user=nfiege@uni-kassel.de
#SBATCH --mail-type=ALL

bench_type=${1}
base_dir=/home/groups/fb16-digi/projects/constmult/sat_mcm/
script_base=benchmark/scripts

if [[ "${bench_type}" == "cmm" ]];
then
  cd ${base_dir};python3 ${script_base}/cmm_benchmark.py 47
fi

if [[ "${bench_type}" == "pcmm" ]];
then
  cd ${base_dir};python3 ${script_base}/pcmm_benchmark.py 47
fi

if [[ "${bench_type}" == "complex" ]];
then
  cd ${base_dir};python3 ${script_base}/complex_mult_benchmark.py 47
fi

if [[ "${bench_type}" == "pcomplex" ]];
then
  cd ${base_dir};python3 ${script_base}/pcomplex_mult_benchmark.py 47
fi

if [[ "${bench_type}" == "conv" ]];
then
  cd ${base_dir};python3 ${script_base}/conv_core_benchmark.py 47
fi

if [[ "${bench_type}" == "pconv" ]];
then
  cd ${base_dir};python3 ${script_base}/pconv_core_benchmark.py 47
fi