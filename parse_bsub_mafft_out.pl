#!/usr/bin/perl
# by xiaenhua @ KIB
# 25/8/2012
################
use Getopt::Std;
getopts("d:o:");

my $dir=$Getopt::Std::opt_d;
my $out=$Getopt::Std::opt_o;
if (!$dir || !$out) {
print "\nUsage: perl $0 -d dir -o mafft.line.out\n";
exit;
}

my(%SEQ,$id,$gid,$sp,$file,@FILE,$str);

opendir(DIR,$dir);

@FILE=grep(/\d+\.out/,readdir(DIR));

foreach $file (sort {$a <=> $b } @FILE) {
	chomp($file);
	open IN,'<',$dir.$file;
        while($str = <IN>){
        chomp($str);
        if($str=~/^>/){
         ($sp,$gid) = (split/\|/,$str)[0,1];
          }
          else{
           $SEQ{$sp}.=$str;
           }
        }
        close IN;
	print $file.'.. done!'."\n";
	}
closedir(DIR);
print "-------------Finished!\n";

open OUT,'>',$out;
foreach (keys %SEQ){
  chomp;
  print OUT "$_"."\n".$SEQ{$_}."\n";
}
close OUT;
