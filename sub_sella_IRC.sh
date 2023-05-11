#!/bin/bash --login

#SBATCH --partition=shared
#SBATCH --job-name="mGs-C_IRC"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --time=07-00:00:00
#SBATCH --constraint=op
#SBATCH -o %j.o 
#SBATCH -e %j.e
#SBATCH --exclusive
##SBATCH --reservation=chemistry_new

cd $SLURM_SUBMIT_DIR

module load nwchem
module load anaconda3/2019.03
conda activate py310

export ASE_NWCHEM_COMMAND="mpirun -np $SLURM_NTASKS nwchem PREFIX.nwi > PREFIX.nwo"
python3 sella_IRC.py