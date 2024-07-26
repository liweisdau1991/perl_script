#!/usr/bin/perl
my (%hash,%fa,$name,$i);
open BED, "<$ARGV[0]"  || die $!;
while(<BED>){
        chomp;
	my $chr =(split/\t/,$_)[0];
	my $id = (split/\t/,$_)[1];
	$fa{$id} = $chr;
}
open FASTA, "<$ARGV[1]"  || die $!;
while(<FASTA>){
	chomp;
	if (/>(\S+)/){
	$name=$1;}
	else{
	chomp($_);
        $hash{$name}.=$_;
        }
}
for (keys %hash) {
        if ($fa{$_}) {
                print ">$fa{$_}\|$_\n";
                print "$hash{$_}\n";
        }
	else{print ">Other\|$_\n";
	print "$hash{$_}\n";
	}
}

