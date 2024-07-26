#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 species.gff species.collinearity 15 > link_for_circos.txt
		⣿⣿⣿⣿⣿⠟⠋⠄⠄⠄⠄⠄⠄⠄⢁⠈⢻⢿⣿⣿⣿⣿⣿⣿⣿
		⣿⣿⣿⣿⣿⠃⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠈⡀⠭⢿⣿⣿⣿⣿
		⣿⣿⣿⣿⡟⠄⢀⣾⣿⣿⣿⣷⣶⣿⣷⣶⣶⡆⠄⠄⠄⣿⣿⣿⣿
		⣿⣿⣿⣿⡇⢀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣧⠄⠄⢸⣿⣿⣿⣿
		⣿⣿⣿⣿⣇⣼⣿⣿⠿⠶⠙⣿⡟⠡⣴⣿⣽⣿⣧⠄⢸⣿⣿⣿⣿
		⣿⣿⣿⣿⣿⣾⣿⣿⣟⣭⣾⣿⣷⣶⣶⣴⣶⣿⣿⢄⣿⣿⣿⣿⣿
		⣿⣿⣿⣿⣿⣿⣿⣿⡟⣩⣿⣿⣿⡏⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
		⣿⣿⣿⣿⣿⣿⣹⡋⠘⠷⣦⣀⣠⡶⠁⠈⠁⠄⣿⣿⣿⣿⣿⣿⣿
		⣿⣿⣿⣿⣿⣿⣍⠃⣴⣶⡔⠒⠄⣠⢀⠄⠄⠄⡨⣿⣿⣿⣿⣿⣿
		⣿⣿⣿⣿⣿⣿⣿⣦⡘⠿⣷⣿⠿⠟⠃⠄⠄⣠⡇⠈⠻⣿⣿⣿⣿
		⣿⣿⣿⣿⡿⠟⠋⢁⣷⣠⠄⠄⠄⠄⣀⣠⣾⡟⠄⠄⠄⠄⠉⠙⠻
		⡿⠟⠋⠁⠄⠄⠄⢸⣿⣿⡯⢓⣴⣾⣿⣿⡟⠄⠄⠄⠄⠄⠄⠄⠄
		⠄⠄⠄⠄⠄⠄⠄⣿⡟⣷⠄⠹⣿⣿⣿⡿⠁⠄⠄⠄⠄⠄⠄⠄⠄
 	
USAGE
if (@ARGV==0){die $usage}

open GFF, $ARGV[0] or die $!;
my %gene_stats;
while (<GFF>) {
	chomp;
	@_ = split /\t/;
	$gene_stats{$_[1]}{"scaffold"} = $_[0];
	$gene_stats{$_[1]}{"start"} = $_[2];
	$gene_stats{$_[1]}{"end"} = $_[3];
}
close GFF;

open MCSCAN, $ARGV[1] or die $!;
my (%strand);
while (<MCSCAN>) {
	if (m/^## Alignment/){
	my $string = (split/\s+/,$_)[2];
	my $id = (split/\:/,$string)[0];
	$strand{$id} = (split/\s+/,$_)[-1];
	}
}
close MCSCAN;

open MCSCAN, $ARGV[1] or die $!;
my (%mcscan, @mcscan);
while (<MCSCAN>) {
	next if m/^#/;
	if (m/^\s*(\d+)-\s*\d+:\s*(\S+)\s*(\S+)/) {
		push @mcscan, $1 if ! exists $mcscan{$1};
		$mcscan{$1}{"size"} ++;
		$mcscan{$1}{"strand"} = $strand{$1};
		$mcscan{$1}{1}{$2} = 1;
		$mcscan{$1}{2}{$3} = 1;
	}
}
close MCSCAN;

foreach my $aligment_code (@mcscan) {
	if ($mcscan{$aligment_code}{"size"} >= $ARGV[2]) {
		my ($scaffold_name, @sites, $scaffold_name2, @sites2);
		foreach (keys %{$mcscan{$aligment_code}{1}}) {
			push @sites, $gene_stats{$_}{"start"};
			push @sites, $gene_stats{$_}{"end"};
			$scaffold_name = $gene_stats{$_}{"scaffold"};
		}
		foreach (keys %{$mcscan{$aligment_code}{2}}) {
                        push @sites2, $gene_stats{$_}{"start"};
                        push @sites2, $gene_stats{$_}{"end"};
                        $scaffold_name2 = $gene_stats{$_}{"scaffold"};
                }

		@sites = sort {$a <=> $b} @sites;
		@sites2 = sort {$a <=> $b} @sites2;
		my $start = $sites[0];
		my $end = $sites[-1];
		print "$scaffold_name\t$start\t$end\t";
		my $start2 = $sites2[0];
		my $end2 = $sites2[-1];
		if ($mcscan{$aligment_code}{"strand"} eq "plus"){
			print "$scaffold_name2\t$start2\t$end2\n";
		}
		elsif ($mcscan{$aligment_code}{"strand"} eq "minus"){
			print "$scaffold_name2\t$end2\t$start2\n";
		}
	}
}
