#!/usr/bin/perl
my (%hash,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
while(<FA>){
        chomp;
        if ($_ =~ />(\S+)/){
                 $name=$1;       
        }else{
        chomp($_);
        $hash{$name}.=$_;
        }
}
open LIST, "<$ARGV[1]"  || die $!;
while(<LIST>){
	chomp;
	next if ($_ =~ /^#LTR_loc/);
	my @array = split/\t/,$_;
	my $ltrid = $array[0];
	my ($chr,$str) = split/\:/,$array[0];
	my ($ltr_st,$ltr_end) = split/\.\./,$str;
	$array[6] =~ s/IN\://g;
	my ($int_st,$int_end) = split/\.\./,$array[6];
	my $superfamily = $array[9];
	my $tetype = $array[10];
	if ($tetype eq "NA"){
	$tetype = "LTR";
	$superfamily = "unknown";
	}
	my $v1 = $int_st - 1;
	my $v2 = $int_end + 1;
	my $ltr1 = $chr . ":" . $ltr_st . ".." . $v1;
	my $int = $chr . ":" . $int_st . ".." . $int_end;
	my $ltr2 = $chr . ":" . $v2 . ".." . $ltr_end;
	print ">$ltrid\_LTR\#$tetype\/$superfamily\n$hash{$ltr1}\n";
	print ">$ltrid\_INT\#$tetype\/$superfamily\n$hash{$int}\n";
	#print ">$ltrid\_LTR\#$tetype\/$superfamily\n$hash{$ltr1}"."$hash{$int}\n";
}
