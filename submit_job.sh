#!/bin/bash
for name in $(cat $1)
do
echo -e "#!/bin/sh
#PBS -N exonerate_$name
#PBS -q cu
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12000:00:00
#PBS -V
#PBS -S /bin/bash
source /opt/intel/compilers_and_libraries_2018.1.163/linux/bin/compilervars.sh intel64
source /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/bin/mklvars.sh intel64
source /opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpivars.sh
cat \$PBS_NODEFILE > /tmp/nodefile.\$\$
echo \"process will start at : \"
date
echo \"++++++++++++++++++++++++++++++++++++++++\"
cd \$PBS_O_WORKDIR
exonerate --model protein2genome --minintron 20 --maxintron 30000 --showtargetgff --showvulgar 0 --showalignment 0 --softmasktarget TRUE --score 100 --percent 70 --ryo '#qi %qi length=%ql alnlen=%qal\\\n#ti %ti length=%tl alnlen=%tal\\\n' --query ./chunk/All_protein.fasta.f2.$name.seq --target Macadamia_ternifolia_softmasked.fasta > ./out/jianguo_exonr_$name.out
echo \"++++++++++++++++++++++++++++++++++++++++\"
echo \"processs will sleep 30s\"
sleep 30
echo \"process end at : \"
date
rm -f /tmp/nodefile.\$\$
rm -f /tmp/nodes.\$\$" > $name.qsub
qsub $name.qsub
done
