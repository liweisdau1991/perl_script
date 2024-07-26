#!/usr/bin/perl
my (%hash,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
my $length = $ARGV[1];
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
        if ($i < $length) {
                print "$_\n";
                print "$hash{$_}\n";
        }
}
close FA;

