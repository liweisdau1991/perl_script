#!/bin/bash
perl -i -pe 's/\*\n/\n/' *.fa
perl -i -pe 's/\.\n/\n/' *.fa
for name in $(cat $1)
do
perl ~/scripts/bin/add_chr.pl /home/liwei/space2/hebei/yangcao/RepeatAnalysis/yangcao.bed ${name} > ${name}.tmp
done
cd ..
perl ~/scripts/bin/bsub_mafft.pl -d Single_Copy_Orthologue_Sequences|sh
perl ~/scripts/bin/parse_bsub_mafft_out.pl -d out/ -o mafft.line.out
~/SoftWare/trimal/trimAl/bin/trimal -in mafft.line.out -out mafft.line.meg -mega -strictplus
nohup perl ~/SoftWare/annotation/RNA-seq/pasa/PASApipeline-v2.5.1/Launch_PASA_pipeline.pl -c alignAssembly.config --ALIGNERS gmap,blat --MAX_INTRON_LENGTH 100000 -C -R -g QZ.final.fasta -t P02TYR19202042-1_r64030_20200218_024115_1_A02.ccs.fasta.fasta --CPU 28 
