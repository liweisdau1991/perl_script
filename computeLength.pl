#!/usr/bin/perl
my (%hash,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
while(<FA>){
        chomp;
        if (/>(\S+)/){
                 $name=$1;       
        }else{
        chomp($_);
        $hash{$name}.=$_;
        }
}
for (sort keys %hash) {
	$i=length($hash{$_});	
	print "$_\t$i\n";
}
close FA;

