#!/bin/bash
#SBATCH --account=b1042 ## Required: your allocation/account name, i.e. eXXXX, pXXXX or bXXXX
#SBATCH --partition=genomicsguestA ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=3:00:00 ## Required: How long will the job need to run (remember different partitions h$
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on per computer/node (default va$
#SBATCH --mem=20G ## how much RAM do you need per computer/node (this affects your FairShare score so b$
#SBATCH --job-name=run_simfiles ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=20240212_fullpopulation_simfiles_0.01_noLD.log ## standard out and standard error goes to this file
#SBATCH --mail-type=ALL ## you can receive e-mail alerts from SLURM about your job


module load singularity

nextflow run prepare_sims.nf --maf 0.00 --ld false --out 20240212_fullpopulation_simfiles_noLD_0.00