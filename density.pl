#!/usr/bin/perl
my $input =shift;
my $window = shift;
open FILE,$input;
while(<FILE>){
        chomp;
        my @a = split;
        my $n = int($a[1]/$window);
        $count{$a[0]}{$n}++;
}
for $chr (keys %count){
        my @max = sort{$a<=>$b}keys %{$count{$chr}};
        for(my $i=0;$i<=$max[$#max];$i++){
                if($count{$chr}{$i}){
                        printf ("%s\t%s\t%s\t%d\n",$chr,$i*$window+1,($i+1)*$window,$count{$chr}{$i});
                }
                else{
                        printf ("%s\t%s\t%s\t0\n",$chr,$i*$window+1,($i+1)*$window);
                }
        }
}

