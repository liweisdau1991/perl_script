#!usr/bin/perl -w
use strict;

#####################################################################
### Concatenate sequences from the same species into
### the same file for 1:1:1 ortholog genes in the order of orthologs
### suitable for orthofinder analysis result
#####################################################################

unless(@ARGV >= 3) {
  print STDERR "perl $0 [MSA dir] [single-copy dir] [output file]\n";
  die "Parameters not enough!$!";
}


my $indir = $ARGV[0];      ### ortholog MSA result directory
my $singledir = $ARGV[1];  ### single-copy ortholog directory
my $outfile = $ARGV[2];    ### output super sequence file


sub readMSAlist {
	opendir INDIR, "$indir" || die "Cannot open this directory!$!";
	my @files = ();
	while(my $file = readdir INDIR) {
		next if($file =~ /^\./);
		if($file =~ /\.fa/) {
			push @files, $file;
		}

	}
	closedir INDIR;
}


### core part
my $files = &readSingleList($singledir);
printf "A total of %s single-copy orthologs\n", scalar @$files;

my %superseqs = ();
foreach my $i(0 .. $#$files) {
	my @headers = ();
	my $seq = &readFasta("$indir/$files->[$i]");
	foreach my $head(sort keys %$seq) {
		$head =~ /([A-Za-z]+)_([A-Za-z]+)_(\S+)/;
		my $species = "$1_$2";
		$superseqs{$species} .= $seq->{$head};
	}
	
}

open OUT, ">$outfile" || die "Cannot write to this file!$!";
foreach my $species(sort keys %superseqs) {
	write2fasta($superseqs{$species}, $species, *OUT);

}
close OUT;


sub readSingleList {
	my $singledir = shift;
	opendir DIR, "$singledir" || die "Cannot open this directory!$!";
	my @singlelist = ();
	while(my $file = readdir DIR) {
		next if($file =~ /^\./);
		push @singlelist, $file;
	}
	close IN;
	return \@singlelist;
}


sub readFasta {
	my $in = shift;
	my %seqs = ();
	my $header;
	open IN, "<$in" || die "Cannot open this file!$!";
	while(<IN>) {
		chomp;
		if(/^>(\S+)/) {
			$header = $1;
		} elsif(/^[^>]+/) {
			$seqs{$header} .= $_;
		}
	}
	close IN;
	return \%seqs;
}


sub write2fasta {
	my ($seq, $header, $fh) = @_;
	print $fh ">$header\n";
	for(my $i=0; $i < length($seq); $i+=60) {
		printf $fh "%s\n", substr($seq, $i, 60);
	}
}
