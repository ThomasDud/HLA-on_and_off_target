#!/bin/bash
#SBATCH --job-name=good_targets_others
#SBATCH --mail-type=END
#SBATCH --mail-user=username@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --output=%x_%A.log #default to <jobname>_<jobid>.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# === ARRAY JOBS ===
# Set the number of jobs in the array using --array
# and replace the output with %x_%A_%a.log
# then the job index becomes available as $SLURM_ARRAY_TASK_ID
#SBATCH --array=1-64
#SBATCH --output=logs/%x_%A_%a.log

# This ensure the script is interrupted when
# -e: any of the steps fails
# -u: a used variable was not set
# -o pipefail: any command in a piped chain fails
set -euo pipefail

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

#Load modules
module load samtools/1.18

# Use the cram list that I provide

file=$(head -n $SLURM_ARRAY_TASK_ID cram.list.txt | tail -1)
fname=$(echo $file | cut -d "/" -f 11)
outdir_idxstats="outdir_idxstats"
mkdir -p $outdir_idxstats

# samtools idxstats

samtools idxstats $file > $outdir_idxstats/${fname}_stats.tsv

# off target

genes=(HLA-A HLA-B HLA-C HLA-DQA1 HLA-DQB1 HLA-DRB1)
outdir_off="off_targets"
mkdir -p $outdir_off

for gene in ${genes[@]}
do
    result=$(grep $gene $outdir_idxstats/${fname}_stats.tsv | awk '{sum += $3} END {print sum}')
    echo -e "$fname\t$result" >> $outdir_off/${gene}_off_sum.tsv
done

# On target

regions=(chr6:29863213-29974447 chr6:31335469-31369490 chr6:31250837-31286044 chr6:32627874-32652782 chr6:32652603-32669077 chr6:32578775-32589847)
outdir_on="on_targets"
mkdir -p $outdir_on

for i in ${!genes[@]}
do
    gene=${genes[$i]}
    region=${regions[$i]}
    idxstats="/localscratch/idxstats_${fname}_${region}.tsv"
    bam="/localscratch/${fname}_${region}_subset.bam"
    samtools view -b --write-index -o $bam $file $region 
    samtools idxstats $bam > $idxstats
    chr=$(echo $region | cut -d ":" -f 1)
    number=$(grep -P "^${chr}\t" $idxstats| cut -f 3)
    echo -e "$fname\t$number" >> $outdir_on/${gene}_target.tsv
    rm $bam $idxstats 
done

#Closing statement
echo "========================"
echo "Completed: $(date)"