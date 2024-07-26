#/usr/bin/perl
my $cmdnum = $ARGV[0];
my $outPrefix = $ARGV[1];
my $genome = $ARGV[2];
my $protein = $ARGV[3];;
system `mkdir chunk out`;
system `perl /share/home/liwei/scripts/divide_fasta.pl -i $protein -od ./chunk -f2 -N $cmdnum`;
open COM, '>', "$outPrefix.command" or die $!;

$cnt = 1;
while($cnt <= $cmdnum){
  print COM "exonerate --model protein2genome --minintron 20 --maxintron 500000 --showtargetgff --showvulgar 0 --showalignment 0 --softmasktarget TRUE --score 100 --percent 70 --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n' --query ./chunk/$protein.f2.$cnt.seq --target $genome > ./out/$outPrefix_exonr_$cnt.out"."\n";
  $cnt++;
  }
