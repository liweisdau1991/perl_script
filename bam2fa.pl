#!/usr/bin/perl
open NAME, "<$ARGV[0]"  || die $!;
while(<NAME>){
	chomp;
	my $file=(split(/\//,$_))[-1];
	my $prefix=(split(/\.bam/,$file))[0];
	print "bam2fastq -u -o $prefix $_\n";
}
