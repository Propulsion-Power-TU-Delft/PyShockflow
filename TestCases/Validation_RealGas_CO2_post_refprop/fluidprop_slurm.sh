#!/bin/sh
#
#SBATCH --job-name=ArFix0.2
#SBATCH --partition=compute
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=research-ae-fpt

cd /home/fneri/FluidProp
./Install_FluidProp


. /home/fneri/FluidProp/SetEnv_FluidProp.sh

cd /scratch/fneri/ASTER_Test_Gamma_0/fluidprop/Validation_RealGas_CO2_post
python main.py 
