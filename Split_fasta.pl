#!/usr/bin/perl


my $file = $ARGV[0];

my($id);
open IN,'<',$file;
while(<IN>){
  chomp;
  if(/^>(\S+)/){
  $id = $1;
  open OUT,'>',$id.'.fasta';
  print OUT $_."\n";
  }
  else{
  print OUT $_."\n";
  }
}
close IN;
