#!/bin/sh
#
#SBATCH -t 5:00:00
#SBATCH -J run-jfock
#SBATCH --partition=dev
#SBATCH --gres gpu:8


# NOTE: These modules are specific to XStream.
module load GCC/4.9.2-binutils-2.25
module load CMake/3.5.2
module load CUDA/9.0.176
module load icc/2015.5.223-GNU-4.9.2-2.25
module load impi/5.0.3.049
module load Python/2.7.10



export LEGION_SRC=/cstor/stanford/toddmtz/users/seemah/legion_master
export REGENT=$LEGION_SRC/language/regent.py
alias regent=$REGENT

#regent top_jfock.rg -L P -i tests/h2o -v tests/h2o/output.dat -p 2
#regent top_jfock.rg -L S -i tests/h2 -v tests/h2/output.dat
#regent top_jfock.rg -L P -i tests/co2 -v tests/co2/output.dat
#regent top_jfock.rg -L D -i tests/fe -v tests/fe/output.dat -p 7 -ll:gpu 7
regent top_jfock.rg -L P -i tests/dna-pair -v tests/dna-pair/output.dat -p 7 -ll:gpu 7
#regent top_jfock.rg -L P -i tests/dna-pair -v tests/dna-pair/output.dat -p 3 -ll:cpu 3
