#!/bin/bash
#SBATCH --account=def-dfuller   # replace this with your own account
#SBATCH --mem-per-cpu=10240M      # memory; default unit is megabytes
#SBATCH --array=21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37
#SBATCH --time=0-11:55           # time (DD-HH:MM)
#SBATCH --ntasks=11
#SBATCH --mail-user=hiroshimamiya@usask.ca
#SBATCH --mail-type=ALL

module load r/4.1.2 
#                                     useArg(binary),sampling(int), ndpost(int), thin(int), core(int),  burnin(int), devLoc(string), ID of data (1:4), note (string)     participant (filepath), scaling of var(bin), drop correlated var (bin),  Index for monitoing (text)                          nTree (Int), bool_looCV
Rscript ~/scratch/BART_test/BART_src/BART_BatchCode.R   1            10000          2000         100         8           2000        "pock"           5                 "LOOCV_2"        $SLURM_ARRAY_TASK_ID                    1          1                     "~/scratch/BART_test/BART_src/monitorIndex.R"       50    1 












