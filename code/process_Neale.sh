#!/bin/bash
#SBATCH --partition sgg
#SBATCH --workdir /data/sgg2/jenny/projects/proxyMR
#SBATCH --job-name process_Neale
#SBATCH --output analysis/process_Neale.out

Rscript ${SGG_generic}/scripts/merge_Neale.r
### reupdate Neale_SGG_dir using ${SGG_generic}/scripts/merge_Neale.r

process_IV_extractions=$(sbatch $SGG_generic/scripts/UKBB/gen_rslurm_Neale_extraction.sh | cut -f 4 -d' ')
until [[ -f ${SGG_generic}/_rslurm_Neale_extraction/submit.sh ]];
do
  sleep 2m
done
get_Neale_IV=$(sbatch --dependency=afterok:$process_IV_extractions ${SGG_generic}/_rslurm_Neale_extraction/submit.sh | cut -f 4 -d' ')

sbatch --dependency=afterok:$get_Neale_IV ${SGG_generic}/scripts/UKBB/MAKE_Neale_IVs_clump.sh
