#!/usr/bin/perl
my $usage = "Usage:\n\tperl sortAndRename_evmPasa_gff3.pl evmPasa.gff3 genePrefix > genome.gff3\n";
if (@ARGV == 0) { die $usage }
my $genePrefix;
if ($ARGV[1]) { $genePrefix = $ARGV[1] }
else          { $genePrefix = "g" }

open GFF, '<', $ARGV[0] or die $!;
my ($geneSym, %geneAnnot);
while ( <GFF> ) {
        next if /^#/;
        next if /^\s/;
        if (/([^\t]+).*?\tgene\t(\d+)\t(\d+)\t(\S+)\t(\S+).*ID=(\S+)\;Name/) { $geneSym = "$6"; 
	$geneAnnot{$geneSym} .= $_ }
        else                         { $geneAnnot{$geneSym} .= $_ }
}
open ID, '<', $ARGV[1] or die $!;
while(<ID>){
	chomp;
	print "$geneAnnot{$_}\n";
}
