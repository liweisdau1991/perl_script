#!/usr/bin/perl
my (%hash,$i);
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
my $chrlist = $ARGV[1];
@array = split/\,/,$chrlist;
for ($i = 0; $i <= $#array; $i++){
	$chr = $array[$i];
	print ">$chr\n$hash{$chr}\n";
}
