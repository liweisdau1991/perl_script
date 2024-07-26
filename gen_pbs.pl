#/usr/bin/perl
my $id=0;
my $file = $ARGV[0];
open IN,'<',$file;
while(<IN>){
  chomp;
  open OUT,'>',$id.'.pbs';
  print OUT "#BSUB -J job_$id
#BSUB -n 1
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R span[hosts=1]
#BSUB -q compute\n";
print OUT "$_\n";
$id++;
}
