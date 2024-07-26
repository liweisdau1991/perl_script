#!/usr/bin/perl
my (%hash,$name,$i);
open FA, "<$ARGV[0]"  || die $!;
print ("Chr\tRepeat_Start\tRepeat_End\tPeriod_Size\tCopy_No.\tAlignment_Score\tConsensus\tRepeat_Sequence\n");
while(<FA>){
	$cnt++;
  	next if /^\#/ || $cnt <= 8;
	if (/Sequence: (.*)/){
	chomp;
		$name = $1;
	}else{
	$hash{$name}.=$_;
	}
}
for (keys %hash) {
	my $chr = $_;
	my @array = split/\n/,$hash{$_};
	for ($i = 0; $i < @array; $i++){
	$line = $array[$i];
	if ($line =~ /^[0-9]/) {
	my @arr = split(" ", $line);
	$rep_start = $arr[0];
        $rep_end = $arr[1];
        $period_size = $arr[2];
        $copy_no = $arr[3];
        $pattern_size = $arr[4]; 
        $percent_match = $arr[5];
                        $percent_indel = $arr[6];
                        $align_score = $arr[7];
                        $a_percent = $arr[8];
                        $c_percent = $arr[9];
                        $g_percent = $arr[10];
                        $t_percent = $arr[11];
                        $entropy = $arr[12];
                        $consensus = $arr[13];
                        $repeat = $arr[14];
                        print ("$chr\t$rep_start\t$rep_end\t$period_size\t$copy_no\t$align_score\t$consensus\t$repeat\n");
	}
}
}
	
