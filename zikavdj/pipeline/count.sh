#!/bin/sh
#
#SBATCH --job-name=FZ_mRNA
#SBATCH --partition=quake,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --time=1-0:0:0
#
# Sample names
#sid=anti-CD137_7dpi-mRNA
#sname=M_GV-7_dpi_pos-anti-CD137_mAb-7_dpi_pos-anti-CD137_mAb-lib1

#sid=isotype_control-mRNA
#sname=M_GV-7_dpi_pos-isotype_control_mAb-7_dpi_pos-isotype_control_mAb-lib1

sid=M_GV-Na_ve-Na_ve-mRNA
sname=M_GV-Na_ve-Na_ve-lib1

echo "Sample ID: $sid"
echo "Sample name: $sname"

fastqpath=${sname%lib1}mRNA
echo "Fastq path: $fastqpath"

echo "Counting mRNA gene expression from sample ID: $sid"
cwd=$(pwd)
cd /oak/stanford/groups/quake/fzanini/sequencing_data/zika_cd137/mRNA
/oak/stanford/groups/quake/fzanini/programs/cellranger-3.0.2/cellranger count \
  --id=$sid \
  --transcriptome=../genome/refdata-cellranger-mm10-3.0.0 \
  --fastqs=../fastq/$fastqpath \
  --sample=$sname \
  --localcores=16 \
  --localmem=128 \
  --nosecondary
cd $cwd
echo 'done!'
