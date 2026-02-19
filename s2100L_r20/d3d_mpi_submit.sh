#!/bin/bash
#SBATCH -J test            # Job name
#SBATCH -n 64                 # Number of MPI tasks (processes/cores)
#SBATCH -N 1                  # Number of nodes to be allocated
#SBATCH -t 8:00:00          # Wall time limit (hh:mm:ss)
#SBATCH -p single             # Partition/queue name
#SBATCH -A hpc_d3d2025     # Allocation/project name
#SBATCH -o delft3dfm.o        # Output file
#SBATCH -e delft3dfm.err      # Error file

date > run.begin
srun --overlap -n 1 singularity exec -B /work /home/admin/singularity/delft3dfm_r142632.sif run_dflowfm.sh autostartstop --partition:ndomains=64:icgsolver=6 FlowFM.mdu 
mpiexec -n 64 singularity exec -B /work:/work /home/admin/singularity/delft3dfm_r142632.sif  dflowfm --autostartstop FlowFM.mdu >out.txt 2>err.txt
date > run.end
exit 0

