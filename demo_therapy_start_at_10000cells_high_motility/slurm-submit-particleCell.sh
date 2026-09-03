#!/bin/sh
#SBATCH --partition=ncpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --job-name="D6e5"
#SBATCH --output=RMM_demo_high_motility.out
#SBATCH --mail-user=xiao.fu@crick.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH --requeue

# #SBATCH --constraint=IB
# #SBATCH --mem=20000

# ------------ variables, paths and modules -----------
echo "SLURM_JOBID = "$SLURM_JOBID
echo "SLURM_JOB_NODELIST = "$SLURM_JOB_NODELIST
echo "SLURM_NNODES = "$SLURM_NNODES
echo "SLURMTMPDIR = "$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

module list

# --------------- launch application -----------------

NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo NPROCS=$NPROCS
echo "Submitting jobs with srun ..."
export OMP_NUM_THREADS=32

srun ./virtualTumour 
wait

echo "Done!"

