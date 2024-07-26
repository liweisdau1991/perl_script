#!/usr/bin/perl

my $file = $ARGV[0];  # repeatMasker.out
my $size = $ARGV[1];  # genome length

my($st,$end,$type,%CNT,$cnt,%LEN,%SUM_cnt,$sum,%SUM_len);
open IN,'<',$file;
while(<IN>){
   chomp;
   s/^\s+//g;
   ($st,$end,$type) = (split/\s+/,$_)[1,2,4];
   #print $_."\n".$type."\t".$st."\t".$end."\n";
   next if $type eq '';
   ($type) = (split/-/,$type)[0];
   $type =~ s/\?//g;
   $CNT{$type}++;
   $LEN{$type}+=(abs($end-$st));
   $sum = $type;
   ($sum) = (split/\//,$sum)[0];
    $SUM_cnt{$sum}++;
    $SUM_len{$sum}+=(abs($end-$st));
    
}
close IN;

$cnt =0;
$sum =0;

foreach (keys %SUM_cnt){
  chomp;
  $cnt+=$SUM_cnt{$_};
  $sum+=$SUM_len{$_};
  print $_."\t".$SUM_cnt{$_}."\t".$SUM_len{$_}."\t".sprintf("%.2f",100*($SUM_len{$_}/$size))."\%\n";
}


foreach (keys %CNT){
  chomp;
  print $_."\t".$CNT{$_}."\t".$LEN{$_}."\t".sprintf("%.2f",100*($LEN{$_}/$size))."\%\n";
}

print "Total\t$cnt\t$sum\t".sprintf("%.2f",100*($sum/$size))."\%\n";

