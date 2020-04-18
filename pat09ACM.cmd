#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J pat09_ACM_endo_1_beat_deleted_files_in_savedir

#SBATCH --nodes=40
#SBATCH --ntasks-per-node=40

#SBATCH -p general
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D /dss/dsshome1/05/di39wun2/Chaste2017
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=peter.marinov@cs.ox.ac.uk
# Wall clock limit:
#SBATCH --time=05:30:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --account=pn34qa
#SBATCH --get-user-env
 
module load slurm_setup
module list 

mpiexec ../build/projects/template_chaste/test/./TestARVC09_ACM --fibrosis 0 --scale_cond 1 --lv_basal_anterior_root 1 --lv_apical_posterior_root 1 --lv_mid_posterior_root 0 --lv_septal_root 1 --rv_mid_anterolateral_root 1 --rv_basal_posterolateral_root 1 --rv_septal_root 1 --save_dir "TestHeart09ACMEndoCalibrated" --regions_file  "/hppfs/work/pn34qa/di39wun2/meshes/chaste09/regions_full.txt" --apexbase "/hppfs/work/pn34qa/di39wun2/meshes/chaste09/apexbase" --mindists "/hppfs/work/pn34qa/di39wun2/meshes/chaste09/" --edge_nodes "/hppfs/work/pn34qa/di39wun2/meshes/chaste09/" 
