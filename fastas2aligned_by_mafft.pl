#!/usr/bin/perl
my ($species,%hash);
my $usage = <<USAGE;
Usage:
	perl $0 OriginalData SingleCopyOrthologGroups 24[this is CPU threads]

USAGE
die $usage if @ARGV != 3;
my @files = `ls $ARGV[0]`;
foreach (@files) {
	chomp;
	next unless m/(\S+)\.pep$/;
	$species = $1;
	open INPUT, "$ARGV[0]/$_" or die $!;
	while(<INPUT>){
	chomp;
	if (/>(\S+)/){$header = $1;
	$hash{$header} = $species;
	#print "$header\t$hash{$header}\n";
	}
	}
	close INPUT;
}
my @files = `ls $ARGV[1]`;

open CMDS, '>', "mafft_cmds" or die $!;

foreach (@files) {
	chomp;
	next unless m/\.fa$/;
	print CMDS "mafft --auto $ARGV[1]/$_ > $ARGV[1]/$_.aln \&\& ";
	print CMDS "seqkit sort $ARGV[1]/$_.aln|seqkit seq -w 0 > $ARGV[1]/$_.format\n";
}

close CMDS;

my $parafly_result =`ParaFly -c mafft_cmds -CPU $ARGV[2] -shuffle -v`;
my @files = `ls $ARGV[1]`;

my %sequences;
my $sequenceId;
my $sequenceLength;

foreach (@files) {
	chomp;
	next unless m/\.format$/;

	open INPUT, "$ARGV[1]/$_" or die $!;

	while (<INPUT>) {
		chomp;
		if (/>(\S+)/) { 
		$header= $1;
		$sequenceId = $hash{$header}; }
		else         { $sequences{$sequenceId} .= $_; $sequenceLength += length $_ }
	}
	close INPUT;
}

my @sequences = keys %sequences;
@sequences = sort { $a cmp $b } @sequences;

open RESULT, '>', 'allSingleCopyOrthologsAlign.fasta' or die $!;

foreach (@sequences) {
	print RESULT ">$_\n$sequences{$_}\n";
}
print "\nAligned length is $sequenceLength\n";
