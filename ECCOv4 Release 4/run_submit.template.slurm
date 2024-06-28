#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
# are set outside this script to allow easy changes and testing.


#SBATCH --cpus-per-task=1

# Replace [budget code] below with your project code (e.g. t01)

#SBATCH --account=n01-OceanBOUND

#SBATCH --export=none

#SBATCH --cpu-freq=2250000

# Setup the job environment (this module needs to be loaded before any other modules)

#. ../scripts/case_setup
module list
ftn --version

# Set the number of threads to 1

#   This prevents any threaded system libraries from automatically 
#   using threading.

export OMP_NUM_THREADS=1
export FI_OFI_RXM_SAR_LIMIT=64K
export MPICH_MPIIO_STATS=1

export SCRATCH=/mnt/lustre/a2fs-nvme/work/n01/n01/afstyles/

# Launch the parallel job
cd ${SCRATCH}/___RUNDIR___

#   srun picks up the distribution from the sbatch options

timestart="$(date +%s)"

#echo srun --distribution=block:cyclic --hint=nomultithread ./mitgcmuv_ad
echo srun --distribution=cyclic:cyclic --hint=nomultithread ./___EXECUTABLE___
#srun --distribution=block:cyclic --hint=nomultithread ./mitgcmuv_ad
srun --distribution=cyclic:cyclic --hint=nomultithread ./___EXECUTABLE___

# Below is just to capture the run details into a file tracesummary.txt

timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo Run end `date` 
echo Run-time seconds $elapsedtotal

stdo=STDOUT.0000
x=$(stat --format=%z $stdo)

tracesummary1="| ${x:0:16} | $SLURM_JOB_ID | $SLURM_JOB_NAME | $SLURM_JOB_NUM_NODES | $SLURM_JOB_QOS | $SLURM_JOB_CPUS_PER_NODE | $SLURM_TASKS_PER_NODE "
        traces2=" |  $MITGCM_FFLAGS $MITGCM_GENM_FLAGS | $elapsedtotal |"
tracesummary=${tracesummary1}${traces2}

echo $tracesummary 
echo $tracesummary > traceSummary.txt

#../scripts/ecse_ecco_check.sh > eec.txt

