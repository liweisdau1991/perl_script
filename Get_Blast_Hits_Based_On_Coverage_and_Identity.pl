#!/usr/bin/perl

my $blast = $ARGV[0]; # blast.ps
my $iden_cutoff = $ARGV[1]; # identity cutoff, such as 80
my $cov_cutoff = $ARGV[2]; # coverage cutoff, such as 80
my $both = $ARGV[3]; # wether use both coverage for query and target, such as yes or no

chomp($both);

if(!$blast || !$iden_cutoff || !$cov_cutoff){
  print "\nperl $0 blast_out.ps identity_cutoff coverage_cutoff both_cov\?\n\n";
  exit(0);
  }

my($qst,$qend,$tst,$tend,$qlen,$tlen,$qiden,$tiden,$qcov,$tcov);

open IN,'<',$blast;
while(<IN>){
  chomp;
  ($qlen,$tlen,$qst,$qend,$tst,$tend,$qiden,$alignlength) = (split/\t/,$_)[1,3,8,9,10,11,13,14];
  next if $qlen !~ /\d+/ || $qst !~ /\d+/;
  $qiden =~ s/\%//;
  $qcov = $qiden*$alignlength/$qlen;
  $tcov = 100*abs($tend - $tst + 1)/$tlen;
  if($both eq 'no'){
    next if $qiden < $iden_cutoff || $qcov < $cov_cutoff;
    print $_."\t".sprintf("%.1f",$qcov)."\%\n";
    }
    else{
     next if $qiden < $iden_cutoff || $qcov < $cov_cutoff || $tcov < $cov_cutoff;
     print $_."\t".sprintf("%.1f",$qcov)."\%\t".sprintf("%.1f",$tcov)."\%\n";
     }
}
close IN;

