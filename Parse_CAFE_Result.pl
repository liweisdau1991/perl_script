#!/usr/bin/perl

my $cafefile = $ARGV[0]; # resultfile.cafe
my $groupfile = $ARGV[1]; # groups.txt
my $targetspecies = $ARGV[2]; # the three-letter short name used in OrthoMCL for the species interest
my $cols = $ARGV[3]; # which two branchs (refer to phylogenic tree) do you want to compare, such as: 3,4

if(!$cafefile || !$groupfile || !$targetspecies || !$cols){
  print "\nperl $0 resultfile.cafe groups.txt OSI 3,4\n\n";
  exit;
}

chomp($cols);
my($sp1, $sp2) = split/\,/,$cols;

$sp1 = $sp1 - 1;
$sp2 = $sp2 - 1;
$pos = $pos - 1;

my($str,@arr,$cnt,$gfid,$gid,$type,%info,$pval,@values);

open IN,'<',$groupfile;
while(<IN>){
  chomp;
  ($gfid,@arr) = split/\s+/,$_;
  $gfid =~ s/\://;
  $info{$gfid} = $_;
}
close IN;

open OUT,'>','Expansion.txt';
open GF,'>','Expanded_GF.txt';
open IN,'<',$cafefile;
while(<IN>){
  chomp;
  $cnt++;
  next if /^\#/ || $cnt <= 11;
  ($gfid,$str) = (split/\s+/,$_)[0,1];
  ($pval) = (split/\s+/,$_)[3];
  $pval =~ s/\(|\)//g;
  @values = split/\,/,$pval;
#  next if $values[$pos] > $pval_cutoff || $values[$pos] !~ /\d+/;
  @arr = $str =~ /\_(\d+)/g;
  #print $str."\n".join("\t",@arr)."\n";
  next if $arr[$sp2] >= $arr[$sp1];
  print GF $gfid."\t".join("\t",@arr)."\t".$values[$pos]."\n";
  ($gfid,@arr) = split/\s+/,$info{$gfid};
  foreach $member (@arr){
   ($type,$gid) = split/\|/,$member;
   next if $type ne $targetspecies;
   print OUT $gid."\n";
   }
}
close IN;
close OUT;
close GF;
open OUT,'>','Contraction.txt';
open GF,'>','Contracted_GF.txt';
open IN,'<',$cafefile;
while(<IN>){
  chomp;
  $cnt++;
  next if /^\#/ || $cnt <= 11;
  ($gfid,$str) = (split/\s+/,$_)[0,1];
  ($pval) = (split/\s+/,$_)[3];
  $pval =~ s/\(|\)//g;
  @values = split/\,/,$pval;
  @arr = $str =~ /\_(\d+)/g;
  next if $arr[$sp2] <= $arr[$sp1];
  print GF $gfid."\t".join("\t",@arr)."\t".$values[$pos]."\n";
  ($gfid,@arr) = split/\s+/,$info{$gfid};
  foreach $member (@arr){
   ($type,$gid) = split/\|/,$member;
   next if $type ne $targetspecies;
   print OUT $gid."\n";
   }
}
close IN;
close OUT;
close GF;

