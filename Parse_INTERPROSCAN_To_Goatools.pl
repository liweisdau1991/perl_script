#!/usr/bin/perl

my $tsv = $ARGV[0]; # interproscan.tsv
my $study = $ARGV[1]; # original gene id file for GSEA; only those with functions are output

if(!$tsv || !$study){
  print "\nperl $0 interproscan.tsv genelist.txt \n\n";
  print "     Nine files (*_study.txt and *_background.txt) are output for goatools\n\n";
  exit;
}

my(@arr,%GO,%YESGO,%KEGGYES,%KEGG,%IPR,%IPRYES,%IPRDES,@tmp,$member,$cnt);

open IN,'<',$tsv;
while(<IN>){
  chomp;
  @arr = split/\t/,$_;
  if($arr[3] eq 'Pfam'){
    $IPR{$arr[0]}{$arr[4]} = $arr[4];
    $IPRDES{$arr[4]} = $arr[5];
    $IPRYES{$arr[0]} = 1;
    }

  if($arr[13] =~ /GO/){
    $GOYES{$arr[0]} = 1;
    @tmp = split/\|/,$arr[13];
    foreach $member (@tmp){
      chomp($member);
      $GO{$arr[0]}{$member} = $member;
      }
    } 

  if($arr[14] =~ /KEGG/){
    $KEGGYES{$arr[0]} = 1;
    @tmp = $arr[14] =~ /KEGG\:\s+(\d+?)\+/g;
    foreach $member (@tmp){
      chomp($member);
      $member = 'ko'.$member;
      $KEGG{$arr[0]}{$member} = $member;
      }
    }
}
close IN;

open GOOUT,'>','go_study.txt';
open IPROUT,'>','ipr_study.txt';
open KEGGOUT,'>','kegg_study.txt';

open IN,'<',$study;
while(<IN>){
  chomp;
  print GOOUT $_."\n" if $GOYES{$_} == 1;
  print IPROUT $_."\n" if $IPRYES{$_} == 1;
  print KEGGOUT $_."\n" if $KEGGYES{$_} == 1;
}
close IN;
close GOOUT;
close IPROUT;
close KEGGOUT;

open GOOUT,'>','go_background.txt';
open IPROUT,'>','ipr_background.txt';
open KEGGOUT,'>','kegg_background.txt';

open OUT1,'>','go_association.txt';
open OUT2,'>','ipr_association.txt';
open OUT3,'>','kegg_association.txt';

foreach (keys %GO){
  chomp;
  print GOOUT $_."\n";
  print OUT1 $_;
  $cnt = 0;
  foreach $member (keys %{$GO{$_}}){
   $cnt++;
   if($cnt == 1){
     print OUT1 "\t".$GO{$_}{$member};
     }
     else{
      print OUT1 ";".$GO{$_}{$member};
      }
   }
   print OUT1 "\n";
}
close OUT1;
close GOOUT;

foreach (keys %IPR){
  chomp;
  print IPROUT $_."\n";
  print OUT2 $_;
  $cnt = 0;
  foreach $member (keys %{$IPR{$_}}){
   $cnt++;
   if($cnt == 1){
     print OUT2 "\t".$IPR{$_}{$member};
     }
     else{
      print OUT2 ";".$IPR{$_}{$member};
      }
   }
   print OUT2 "\n";
}
close OUT2;
close IPROUT;

foreach (keys %KEGG){
  chomp;
  print KEGGOUT $_."\n";
  print OUT3 $_;
  $cnt = 0;
  foreach $member (keys %{$KEGG{$_}}){
   $cnt++;
   if($cnt == 1){
     print OUT3 "\t".$KEGG{$_}{$member};
     }
     else{
      print OUT3 ";".$KEGG{$_}{$member};
      }
   }
   print OUT3 "\n";
}
close OUT3;
close KEGGOUT;

open OUT,'>','IPR_DES.txt';
foreach (keys %IPRDES){
  chomp;
  print OUT $_."\t".$IPRDES{$_}."\n";
}
close OUT;
