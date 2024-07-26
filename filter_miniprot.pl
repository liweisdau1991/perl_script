#!/usr/bin/perl
open GFF, '<', $ARGV[0] or die $!;
my ($geneSym, %geneAnnot);
while ( <GFF> ) {
	next if /^#/;
	next if /^\s/;
	if (/([^\t]+).*?\tmRNA\t(\d+)\t(\d+)\t(\S+)\t(\S+)/) { $geneSym = "$1\t$2\t$3\t$5";
	$geneAnnot{$geneSym} .= $_ }
	else { $geneAnnot{$geneSym} .= $_ }
}
for (keys %geneAnnot) {
	my @array = split /\n/, $geneAnnot{$_};
	($identity,$coverage) = $array[0] =~ /Identity\=(\d+\.\d+)\;Positive\=(\d+\.\d+)/;
	print "$geneAnnot{$_}" if (($identity >= 0.9) && ($coverage >=0.9));
}
