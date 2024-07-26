#!/usr/bin/perl
my ($chr,$annot,@array);
open BED, "<$ARGV[0]"  || die $!;
my $chrlist = $ARGV[1];
while(<BED>){
        chomp;
	$annot = $_;
	$chr = (split/\t/,$_)[0];
	@array = split/\,/,$chrlist;
	if ( grep { $_ eq $chr } @array) {
	print "$annot\n";
	}
}

