#/usr/bin/perl
open GFF, '<', $ARGV[0] or die $!;
my ($geneSym, %geneAnnot,$tmp,$chr,$geneid);
#
while(<GFF>){
    next if /^#/;
      next if /^\s/;
        if (/([^\t]+).*?\tgene\t(\d+)\t(\d+)\t(\S+)\t(\S+).*ID=(\S+)\;Name/) { $geneSym = $6; 
        	$geneAnnot{$geneSym} .= $_ }
        	        else                         { $geneAnnot{$geneSym} .= $_ }
        	        }
my $geneid = $ARGV[1];
print "$geneAnnot{$geneid}";
