#### DONT SUBMIT THIS SCRIPT ITSELF TO SLURM ####

#activate conda env
module load anaconda3/2022.10
conda activate snakemake

#tmux stuff
systemd-run --scope --user tmux #start new window locally that will continue to run after logging out
tmux ls #list existing windows
tmux attach -t $(int) #attach existing window
#ctrl + b (release), d to exit
#ctrl + b (release) :kill-session to kill session
#ctrl + b (release) [ to scroll

#make sure tmux is in the right directory and has conda env activated
#use the below command to submit snakemake jobs as separate slurm jobs
#adjust the --jobs parameter as needed (especially important for the variant calling step: there are a total of 690 jobs to run for the callHaps rule alone so you want to make sure you can run lots of stuff in parallel)
#latency wait is in seconds: check the rules in the snakefile to adjust accordingly (latency wiat should be > longest expected run time for a job)
#run entire pipeline
snakemake --executor cluster-generic --cluster-generic-submit-cmd "sbatch -N {resources.nodes} -p {resources.slurm_partition} {resources.slurm_extra}" \
--cluster-generic-cancel-cmd scancel --use-conda --jobs 30 --immediate-submit --notemp --printshellcmds --latency-wait 21600 --rerun-triggers mtime

#run with target rule
snakemake -R --until "INSERT TARGET RULE HERE" --executor cluster-generic --cluster-generic-submit-cmd "sbatch -N {resources.nodes} -p {resources.slurm_partition} {resources.slurm_extra}" \
--cluster-generic-cancel-cmd scancel --use-conda --jobs 30 --immediate-submit --notemp --printshellcmds --latency-wait 21600 --rerun-triggers mtime

#dry run
snakemake --executor cluster-generic --cluster-generic-submit-cmd "sbatch -N {resources.nodes} -p {resources.slurm_partition} {resources.slurm_extra}" \
--cluster-generic-cancel-cmd scancel --use-conda --jobs 30 --immediate-submit --notemp --printshellcmds --latency-wait 21600 --rerun-triggers mtime -n




