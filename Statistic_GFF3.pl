#!/usr/bin/perl

my $file = $ARGV[0]; # gff3 file;

my(@arr,$cnt,$len,$gcnt,$glen,$mcnt,$mlen,$ccnt,$clen,$ecnt,$elen,$icnt,$ilen,%CDS,$id);

open IN,'<',$file;
while(<IN>){
  chomp;
  @arr = split/\t/,$_;
  if($arr[2] eq 'gene'){
    $gcnt++;
    $glen+=(abs($arr[4]-$arr[3])+1);
   }

  if($arr[2] eq 'mRNA'){
    $mcnt++;
    $mlen+=(abs($arr[4]-$arr[3])+1);
   }

  if($arr[2] eq 'CDS'){
    $ccnt++;
    $clen+=abs($arr[4]-$arr[3]);
    ($id) = $arr[8] =~ /Parent=(.*)$/;
    $CDS{$id}+=(abs($arr[4]-$arr[3])+1);
   }

  if($arr[2] eq 'exon'){
    $ecnt++;
    $elen+=(abs($arr[4]-$arr[3])+1);
   }
}
close IN;

foreach (keys %CDS){
 $cnt++;
 $len+=$CDS{$_};
 }

$icnt = $ecnt - $gcnt;  
$ilen = $glen - $elen;
#print "Total number of genes: $gcnt\n";
print "Total number of genes: $gcnt\nTotal number of models: $mcnt\nAverage gene length (bp): ".sprintf("%.0f",$glen/$gcnt)."\nAverage CDS length (bp): ".sprintf("%.0f",$len/$mcnt)."\nAverage exon length (bp): ".sprintf("%.0f",$elen/$ecnt)."\nAverage exon per gene: ".sprintf("%.1f",$ecnt/$mcnt)."\nAverage intron length (bp): ".sprintf("%.0f",$ilen/$icnt)."\n";
