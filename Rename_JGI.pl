#!/usr/bin/perl

my $file = $ARGV[0];

my($id);
open IN,'<',$file;
while(<IN>){
  chomp;
  if(/^>/){
  ($id) = (split/\|/,$_)[-2];
  print '>'.$id."\n";
  }
  else{
  print $_."\n";
  }
}
close IN;
