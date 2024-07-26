#!/usr/bin/perl
my (%hash,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
while(<FA>){
        if (/>(\S+)/){
                 $name=$1;       
        }else{
        $hash{$name}.=$_;
        }
}
open LOCUS, "<$ARGV[1]"  || die $!;
while(<LOCUS>){
	my $id=(split/\t/, $_)[2];
	my $te=(split/\t/, $_)[6];
	print ">$id\n"."$hash{$id}" if ($te eq "N");
}
