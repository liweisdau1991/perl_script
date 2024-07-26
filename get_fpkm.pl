#!/usr/bin/perl
my ($id,$gtf,$geneid,$fpkm,%hash,@array);
open ID, "<$ARGV[0]"  || die $!;
while(<ID>){
        chomp;
        $id = (split/\s+/,$_)[0];
	push @array, $id;
        $gtf = (split/\s+/,$_)[1];
	open(GTF,"$gtf.tab");
	while(<GTF>){
	chomp;
	next if ($_ =~ /^Gene/);
	my @array = split/\t/,$_;
	$geneid = $array[0];
	$tpm = $array[8];
	$hash{$geneid}{$id} = $tpm;
	}
	close GTF;
}
#
print "geneID\t";
for ($i=0; $i <= $#array; $i++){
	print "$array[$i]\t";
}
print "\n";
foreach my $key (keys %hash){
	print "$key\t";
	for ($i=0; $i <= $#array; $i++){
	my $value = $array[$i];
	print "$hash{$key}{$value}\t";
	}
print "\n";
}
