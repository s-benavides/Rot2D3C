#!/bin/bash
#SBATCH --nodes=1           # Number of nodes (1 node = 16 cores)
#SBATCH --tasks-per-node=1     # 16 cores per node?
##SBATCH -p newnodes      # partition name
#SBATCH -p sched_mit_hill
#SBATCH --time=0-12:00:00   # 1 day and 3 hours
##SBATCH -p sched_any_quicktest
##SBATCH --time=0-00:15:00   # 1 day and 3 hours
#SBATCH -J Rot2D3C_backup  # sensible name for the job


set -eu   ## Stop on errors and on undefined variables

echo Beginning

~/.local/bin/rsync --stats --append-verify --exclude './test*' -azhR ./* /pool001/santiago_b/Rot2D3C/

echo done

