#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH --threads-per-core=1
#SBATCH -c 1
#SBATCH --mem-per-cpu=50GB
#SBATCH --requeue
#SBATCH --ntasks=1

#SBATCH -A dror-account
#SBATCH -p drorq

#SBATCH --job-name=$1

#SBATCH --output=out.txt
#SBATCH --error=err.txt

##SBATCH --open-mode=append

##SBATCH --array=


## script for job

trap 'scontrol requeue ${SLURM_JOB_ID}; exit 15' 15
srun ~/.juliaup/bin/julia run_pluto.jl $1 &
wait
EOT
