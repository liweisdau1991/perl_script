#!/usr/bin/perl
my (%hash,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
while(<FA>){
        chomp;
        if (/>/){
                 $name=$_;       
        }else{
        chomp($_);
        $hash{$name}.=$_;
        }
}
for (sort keys %hash) {
            $i=length($hash{$_});
        if ($i >=1) {
                print "$_\n";
                print "$hash{$_}\n";
        }
}
close FA;

