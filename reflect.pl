#!/usr/bin/perl
my (%hash,%fasta,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
while(<FA>){
        chomp;
	my ($id1,$id2) = split/\t/,$_;
	$hash{$id2} = $id1;
}
open FB, "<$ARGV[1]"  || die $!;
while(<FB>){
        chomp;
        if (/>(\S+)/){
                 $name=$1;
        }else{
        chomp($_);
        $fasta{$name}.=$_;
        }
}
for (sort keys %hash) {
                print ">$hash{$_}\n";
                print "$fasta{$_}\n";
}
