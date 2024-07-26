#!/usr/bin/perl
my (%hash,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
while(<FA>){
        if (/>PGA\_scaffold\_(\d+)\_\_(\d+)\_contigs\_\_length\_(\d+)/){
                 $scaffold=$1; 
		 $contig = $2;
		 $length = $3;
		print "\>Chr$1\n" if ($length > 10000000);
		print "\>Scaffold$1\n" if ($length < 10000000);      
        }else{
	print "$_";        
}
}
close FA;

