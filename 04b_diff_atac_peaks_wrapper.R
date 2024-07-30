"%&%" <- function(a,b) paste(a,b, sep = "")
wd <- '/project2/gilad/daraujo/scRNA_scATAC'

for (c in seq(0:26)){
  command_line <- '#!/bin/sh' %&% '\n' %&% '#SBATCH --time=36:00:00' %&% '\n' %&%
    '#SBATCH --mem=400G' %&% '\n' %&% '#SBATCH --partition=bigmem2' %&% '\n' %&%
    '#SBATCH --account=pi-gilad' %&% '\n' %&% '#SBATCH --error=' %&%
    wd %&% '/DApeaks_cluster_' %&% as.character(c-1) %&% '.error' %&% '\n' %&%
    '#SBATCH --out=' %&% wd %&% '/DApeaks_cluster_' %&% as.character(c-1) %&% '.out' %&% '\n\n' %&%
    'cd /project2/gilad/daraujo/01_scripts/ \n' %&% 'module load R/4.2.0' %&% '\n' %&%
    'Rscript 04a_diff_atac_peaks.R -c ' %&% as.character(c-1)
  
  cat(command_line, file='DApeaks_cluster_'%&%as.character(c-1)%&%'.sbatch')
  system('sbatch DApeaks_cluster_'%&%as.character(c-1)%&%'.sbatch')
  print('submitted: DApeaks_cluster_'%&%as.character(c-1)%&%'.sbatch')
}