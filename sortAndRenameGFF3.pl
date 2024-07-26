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
        if (/([^\t]+).*?\tgene\t(\d+)\t(\d+)\t(\S+)\t(\S+)/) { $geneSym = "$1\t$2\t$3\t$5"; 
	$geneAnnot{$geneSym} .= $_ }
        else                         { $geneAnnot{$geneSym} .= $_ }
}

my @geneSym = keys %geneAnnot;
my (%geneSymSort1, %geneSymSort2, %geneSymSort3, %geneSymSort4);
foreach ( @geneSym ) { $geneSymSort1{$_} = $1 if /(\d+)/  }
foreach ( @geneSym ) { $geneSymSort2{$_} = $1 if /\t(\d+)/ }
foreach ( @geneSym ) { $geneSymSort3{$_} = $2 if /\t(\d+)\t(\d+)/ }
foreach ( @geneSym ) { $geneSymSort4{$_} = $4 if /\t(\d+)\t(\d+)\t(\S)\t(\S)/ }
@geneSym = sort { $geneSymSort1{$a} <=> $geneSymSort1{$b} or $geneSymSort2{$a} <=> $geneSymSort2{$b} or $geneSymSort3{$a} <=> $geneSymSort3{$b} or $geneSymSort4{$a} <=> $geneSymSort4{$b} } @geneSym;

my ($geneNum,$scaffoldName,$scaffoldGeneNum);
foreach (@geneSym) {
        my $geneAnnot = $geneAnnot{$_};
        $geneNum ++;
        $geneNum = '00000'."$geneNum" if length $geneNum == 1;
        $geneNum = '0000'."$geneNum" if length $geneNum == 2;
        $geneNum = '000'."$geneNum" if length $geneNum == 3;
        $geneNum = '00'."$geneNum" if length $geneNum == 4;
        $geneNum = '0'."$geneNum" if length $geneNum == 5;
        my $scaffoldNum = $1 if /(\d+)/;
        my $scaffoldname = $1 if /(.*?)\t/;
        if ($scaffoldName eq $scaffoldname) { $scaffoldGeneNum ++ }
        else                                { $scaffoldName=$scaffoldname; $scaffoldGeneNum = 1 }

        my @geneAnnot = split /\n/, $geneAnnot;
        my ($transcriptIdNum, $fiveUTR_Num, $threeUTR_Num, $exonNum, $cdsNum);
        foreach (@geneAnnot) {
                if (/\tgene\t/) { s/ID=.*/ID=$genePrefix$geneNum;Name=$genePrefix$geneNum/; print "$_\n" }
                if (/\tmRNA\t/) { $transcriptIdNum ++; $fiveUTR_Num=0; $threeUTR_Num=0; $exonNum=0; $cdsNum=0; s/ID=.*/ID=$genePrefix$geneNum\.t$transcriptIdNum;Parent=$genePrefix$geneNum;Name=$genePrefix$geneNum\.t$transcriptIdNum/; print "$_\n" }
                if (/\tfive_prime_UTR\t/)  { $fiveUTR_Num ++; s/ID=.*/ID=$genePrefix$geneNum\.t$transcriptIdNum.utr5p$fiveUTR_Num;Parent=$genePrefix$geneNum\.t$transcriptIdNum/; print "$_\n" }
                if (/\tthree_prime_UTR\t/) { $threeUTR_Num ++; s/ID=.*/ID=$genePrefix$geneNum\.t$transcriptIdNum.utr3p$threeUTR_Num;Parent=$genePrefix$geneNum\.t$transcriptIdNum/; print "$_\n" }
                if (/\texon\t/) { $exonNum ++; s/ID=.*/ID=$genePrefix$geneNum\.t$transcriptIdNum.exon$exonNum;Parent=$genePrefix$geneNum\.t$transcriptIdNum/; print "$_\n" }
                if (/\tCDS\t/)  { $cdsNum ++; s/ID=.*/ID=$genePrefix$geneNum\.t$transcriptIdNum.cds$cdsNum;Parent=$genePrefix$geneNum\.t$transcriptIdNum/; print "$_\n" }
        }
                print "\n";
}

