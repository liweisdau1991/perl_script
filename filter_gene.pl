#!/usr/bin/perl
my ($geneSym,%geneAnnot);
open GFF, '<', $ARGV[0] or die $!;
my ($geneSym, %geneAnnot);
while ( <GFF> ) {
        next if /^#/;
        next if /^\s/;
	$chr = (split/\t/,$_)[0];
	$type = (split/\t/,$_)[2];
	$start = (split/\t/,$_)[3];
	$end = (split/\t/,$_)[4];
	$strand = (split/\t/,$_)[6];
        if ($type !~ /CDS/ && $type !~ /exon/ && $type !~ /mRNA/)
	{$geneSym = "$chr\t$type\t$start\t$end\t$strand"; $geneAnnot{$geneSym} .= $_ }
        else                         { $geneAnnot{$geneSym} .= $_ }
}
foreach (keys %geneAnnot){
	$type = (split/\t/,$_)[1];
	print "$geneAnnot{$_}" if ($type eq "gene");
}
