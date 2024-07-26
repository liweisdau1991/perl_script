#!/usr/bin/perl
use List::Util qw( min max );
my (%type,%number,%length,%percent);
open FA, "<$ARGV[0]"  || die $!;
while(<FA>){
	chomp;
	my @array = split/\t/,$_;
	my $type = $array[0];
	my $number = $array[1];
	my $length = $array[2];
	my $percent = $array[3];
#	my $percent =~ s/%//;
	if (!$type{$type}){
	$type{$type} = $type;
	$number{$type} = $number;
	$length{$type} = $length;
	$percent{$type} = $percent;
	} elsif ($type{$type}){
	$type{$type} = $type;
	$number{$type} = max($number{$type},$$number);
        $length{$type} = max($length{$type},$length);
        $percent{$type} = max($percent{$type},$percent);
	}
}
$type{"LTR/Other"}="LTR/Other";
$number{"LTR/Other"}=$number{"LTR"} - $number{"LTR/Gypsy"} - $number{"LTR/Copia"};
$length{"LTR/Other"}=$length{"LTR"} - $length{"LTR/Gypsy"} - $length{"LTR/Copia"};
$percent{"LTR/Other"}=($percent{"LTR"} - $percent{"LTR/Gypsy"} - $percent{"LTR/Copia"})."%";
$type{"ClassI"}="ClassIOther";
$number{"ClassI"} = $number{"SINE"}+$number{"LINE"};
$length{"ClassI"} = $length{"SINE"}+$length{"LINE"};
$percent{"ClassI"} = ($percent{"SINE"}+$percent{"LINE"})."%";
$type{"OtherRepeat"}="OtherRepeats";
$number{"OtherRepeat"}=$number{"Total"} - $number{"DNA"} - $number{"LTR"} - $number{"SINE"} - $number{"LINE"};
$length{"OtherRepeat"}=$length{"Total"} - $length{"DNA"} - $length{"LTR"} - $length{"SINE"} - $length{"LINE"};
$percent{"OtherRepeat"}=($percent{"Total"} - $percent{"DNA"} - $percent{"LTR"} - $percent{"SINE"} - $percent{"LINE"})."%";
my @new_array = ("LTR/Gypsy","LTR/Copia","LTR/Other","ClassI","DNA","OtherRepeat","Total");
foreach my $value (@new_array){
	print "$type{$value}\t$number{$value}\t$length{$value}\t$percent{$value}\n";
}
