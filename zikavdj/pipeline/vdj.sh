#!/bin/sh
#
#SBATCH --job-name=FZ_VDJ
#SBATCH --partition=quake,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --time=1-0:0:0
#
# Sample names
#sid=anti-CD137_7dpi-VDJ
#sname=M_GV-7_dpi_pos-anti-CD137_mAb-7_dpi_pos-anti-CD137_mAb-lib2

#sid=isotype_control-VDJ
#sname=M_GV-7_dpi_pos-isotype_control_mAb-7_dpi_pos-isotype_control_mAb-lib2

sid=M_GV-Na_ve-Na_ve-VDJ
sname=M_GV-Na_ve-Na_ve-lib2

echo "Sample ID: $sid"
echo "Sample name: $sname"

fastqpath=${sname%lib2}VDJ
echo "Fastq path: $fastqpath"

echo "Assembling VDJ from sample ID: $sid"
cwd=$(pwd)
cd /oak/stanford/groups/quake/fzanini/sequencing_data/zika_cd137/vdj
/oak/stanford/groups/quake/fzanini/programs/cellranger-3.0.2/cellranger \
  vdj \
  --id=$sid \
  --reference=../genome/vdj_IMGT_mouse \
  --fastqs=../fastq/$fastqpath \
  --sample=$sname \
  --localcores=16 \
  --localmem=128
cd $cwd
echo 'done!'
