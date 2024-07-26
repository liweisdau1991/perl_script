use Getopt::Long;


if(!(-e $ARGV[0])){
	&Usage(),exit(1);
}

$ARGV_len=@ARGV;
@ARGV_info;
for(my $i=1; $i<$ARGV_len; $i++){
	@ARGV_info[$i-1]=$ARGV[$i];
}

my $lineIN='';  #input file read
my $outf=$ARGV[0].'.out';

open(OUT,'>',"$outf");
open(IN,'<',"$ARGV[0]");

my @db;

#read in file to @db; 
print "read in data...\n";

	#use first value to cut
	#out
	#   897   10.3  7.3  3.2  scaffold1           1     151 (9561602) + rnd-3_family-127   Unknown               392    548   (382)      1 *
	#   596   25.4  1.1  0.6  scaffold1       26486   26663 (9535090) C rnd-5_family-1105  Unknown            (1967)    940     762     33 *
	#0-value 1-2-3- 4-chrom 5-st 6-ed 7- 8-c/+ 9-name 10-type


my $count=0;
my @row;
my $max=100000000000000;

$lineIN=<IN>;
$lineIN=<IN>;
$lineIN=<IN>;

#$out
#0-value 1-2-3- 4-chrom 5-st 6-ed 7- 8-c/+ 9-name 10-type
while($lineIN=<IN>){
	chomp($lineIN);
  my @tp= &Getheadinfo($lineIN);

  my $len=$tp[6]-$tp[5]+1;
  #mid #0-sortkey 1-chrom 2-st 3-ed 4-len 5-type 6-name 7-value   
	my $tpadd=$max+$tp[5]; 	
	#adjest sortkey DNA/RC/LINE/SINE score add 0 behand
	if(index($tp[10],'DNA/',0)>=0 ||index($tp[10],'RC/',0)>=0 ||index($tp[10],'LINE',0)>=0 ||index($tp[10],'SINE',0)>=0){
		$tp[0]=$tp[0]*10;
	}
	if(index($tp[10],'Unknown',0)>=0){
		$tp[0]=$tp[0]/10;
	}
	

	$row[$count]="$tp[4]-$tpadd\t$tp[4]\t$tp[5]\t$tp[6]\t$len\t$tp[10]\t$tp[9]\t$tp[0]";
	#print "$row[$count]\n";
	$count++;
	
}

#sort
my $otm=$ARGV[0].".mid";
open(OUTM,'>',"$otm");   
print "$count lines readed\nsort...\n";   
@row=sort(@row);
for(my $i=0;$i<$count;$i++){
	my @tp= &Getheadinfo($row[$i]);		
	print OUTM "$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]\n";
  	$db[$i]=[@tp];
  	#print "$db[$i][6]bbb$db[$i][7]aaa\n";
} 
close(OUTM);
print "sorted data save to $otm\n";
	  	
#loop clean
my %o;  #count overlap
my $ot=$ARGV[0].".detail";
open(OUTD,'>',"$ot");  
  my $d=9999;
  my $loop=0;
  while($d>0){
  	$d=0;
  	$loop++;
  	my $db=@db;    	
  	print "Loop $loop: $tp\t";
  	print OUTD "##Loop $loop: $tp-----------------------\n";
  	#mid #0-sortkey 1-chrom 2-st 3-ed 4-len 5-type 6-name 7-value   
  	#print $tp[0][0]."\n";
  	my @tpdb;  #used for sort
  	my $tpc=0;
  	for(my $i=1;$i<$db;$i++){
  		my $add=1;
  		if($db[$i][1] eq $db[$i-1][1]){	#same chrom
  			if($db[$i][2] < $db[$i-1][3]){   #overlap
  				$d++;
  				if($db[$i][5] eq $db[$i-1][5]){
  					#combin
  					print OUTD "n\t$db[$i-1][1]\t$db[$i-1][2]\t$db[$i-1][3]\t$db[$i-1][4]\t$db[$i-1][5]\t$db[$i-1][6]\t$db[$i-1][7]\n";
  					$add=0;	
  					$db[$i][2]=$db[$i-1][2];
  					$db[$i][3]=$db[$i][3]>$db[$i-1][3]?$db[$i][3]:$db[$i-1][3];  
  					$db[$i][4]=$db[$i][3]-$db[$i][2];
  					$db[$i][7]=$db[$i][7]>$db[$i-1][7]?$db[$i][7]:$db[$i-1][7]; 
  					$db[$i][6]=$db[$i][5]."_nest";    					
  				}else{
  					if($db[$i][2] == $db[$i-1][2]){ #leave big one to i, remove i-1
  						print OUTD "i\t$db[$i-1][1]\t$db[$i-1][2]\t$db[$i-1][3]\t$db[$i-1][4]\t$db[$i-1][5]\t$db[$i-1][6]\t$db[$i-1][7]\n";    					
						if($o{$db[$i-1][5]}>0){
							$o{$db[$i-1][5]}+=$db[$i-1][4];
						}else{
							$o{$db[$i-1][5]}=$db[$i-1][4];
						}  
  						$db[$i][3]=$db[$i][3]>$db[$i-1][3]?$db[$i][3]:$db[$i-1][3];   
  						$db[$i][4]=$db[$i][3]-$db[$i][2];
  						$add=0;								
  					}else{  					
  						if($db[$i][3] < $db[$i-1][3]){	#rm i
  							print OUTD "i\t$db[$i][1]\t$db[$i][2]\t$db[$i][3]\t$db[$i][4]\t$db[$i][5]\t$db[$i][6]\t$db[$i][7]\n";
							if($o{$db[$i][5]}>0){
								$o{$db[$i][5]}+=$db[$i][4];
							}	else{
								$o{$db[$i][5]}=$db[$i][4];
							} 
  							for(my $j=0;$j<8;$j++){
  								$db[$i][$j]=$db[$i-1][$j];
  							}
  							$db[$i][4]=$db[$i][3]-$db[$i][2];
  							$add=0;
  						}else{
  							#cut
  							if($db[$i][7] > $db[$i-1][7]){
  								print OUTD "d\t$db[$i-1][1]\t$db[$i-1][2]\t$db[$i-1][3]\t$db[$i-1][4]\t$db[$i-1][5]\t$db[$i-1][6]\t$db[$i-1][7]\n";
								if($o{$db[$i-1][5]}>0){
									$o{$db[$i-1][5]}+=$db[$i-1][4];
								}	else{
									$o{$db[$i-1][5]}=$db[$i-1][4];
								}  
  								$db[$i-1][3]=$db[$i][2]-1;
  								$db[$i-1][4]=$db[$i-1][3]-$db[$i-1][2]+1;   					
  							}else{
  								print OUTD "d\t$db[$i][1]\t$db[$i][2]\t$db[$i][3]\t$db[$i][4]\t$db[$i][5]\t$db[$i][6]\t$db[$i][7]\n"; 	
  								#print "d";
								if($o{$db[$i][5]}>0){
									$o{$db[$i][5]}+=$db[$i][4];
								}	else{
									$o{$db[$i][5]}=$db[$i][4];
								} 					
   								$db[$i][2]=$db[$i-1][3]+1;
  								$db[$i][4]=$db[$i][3]-$db[$i][2]+1;  
  								$add=2;
  								  							
  							}
  						}
  					}
  				}
  			}    			
  		}
  		if($add ==1){
			my $tpadd=$max+$db[$i-1][2]; 	
  			$tpdb[$tpc]="$db[$i-1][1]-$tpadd\t$db[$i-1][1]\t$db[$i-1][2]\t$db[$i-1][3]\t$db[$i-1][4]\t$db[$i-1][5]\t$db[$i-1][6]\t$db[$i-1][7]"; 	    			
  			#print "length $l last line $db[$i-1][1]-$tpadd\n";
  			$tpc++;  			
  		}
  		if($add ==2){
			my $tpadd=$max+$db[$i-1][2]; 	
  			$tpdb[$tpc]="$db[$i-1][1]-$tpadd\t$db[$i-1][1]\t$db[$i-1][2]\t$db[$i-1][3]\t$db[$i-1][4]\t$db[$i-1][5]\t$db[$i-1][6]\t$db[$i-1][7]"; 	    			
  			$tpc++;
			$tpadd=$max+$db[$i][2]; 	
  			$tpdb[$tpc]="$db[$i][1]-$tpadd\t$db[$i][1]\t$db[$i][2]\t$db[$i][3]\t$db[$i][4]\t$db[$i][5]\t$db[$i][6]\t$db[$i][7]"; 	    			
  			$tpc++;
  			$i++;
  		}
  	} 
	my $tpadd=$max+$db[$db-1][2]; 	
  	$tpdb[$tpc]="$db[$db-1][1]-$tpadd\t$db[$db-1][1]\t$db[$db-1][2]\t$db[$db-1][3]\t$db[$db-1][4]\t$db[$db-1][5]\t$db[$db-1][6]\t$db[$db-1][7]"; 	    			
  	$tpc++; 
  	
  	#print "$tpdb[tpc]\n";
  	   	
  	#print $tp[3][3]."\n";
  	@tpdb=sort(@tpdb);
  	
  	@db=();
   	for(my $i=0;$i<$tpc;$i++){
		my @tp= &Getheadinfo($tpdb[$i]);	
		$db[$i]=[@tp];
		#print "first $db[$i][0]\n"; 	
  	}    	
  	print "delete $d\n";	    	 
  }
  print "output...\n";
  my $db=@db;
  my %a;  #total count
  #mid #0-sortkey 1-chrom 2-st 3-ed 4-len 5-type 6-name 7-value   
  my $countA=0;
  for(my $i=0;$i<$db;$i++){
  	print OUT "$db[$i][1]\t$db[$i][2]\t$db[$i][3]\t$db[$i][4]\t$db[$i][5]\t$db[$i][6]\t$db[$i][7]\n";
  	$countA+=$db[$i][4];
  	if($a{$db[$i][5]}>0){
  		$a{$db[$i][5]}=$a{$db[$i][5]}+$db[$i][4];
  	}else{
  		$a{$db[$i][5]}=$db[$i][4];
  	}
  }
  print OUTD "\n##Overlap---------------\n";
  my @t=keys(%o);
  @t=sort(@t);
  my $t=@t;
  for(my $i=0;$i<$t;$i++){
  	print OUTD "$t[$i]\t$o{$t[$i]}\n";
  	$countD+=$o{$t[$i]};
  	
}  
	print OUTD "\n##Type Count---------------\n";
    my @tc=keys(%a);
    my $tc=@tc;
    @tc=sort(@tc);
    for(my $i=0;$i<$tc;$i++){
    	print OUTD "$tc[$i]\t$a{$tc[$i]}\n";
    	$countD+=$a{$tc[$i]};    	
	} 
	print OUTD "\n---------------\nTotal length\t$countA\nDelete Length\t$countD\n";    
	print "Total length\t$countA\nDelete Length\t$countD\n";      
    close(OUTD);


sub Getheadinfo{
	my $test = -1;
	my $li = 0;
	my @samble = (" ", "\t");
	$_=@_;
	if($_>1){
		for (my $i=1;$i<$_;$i++){
			$samble[$i+1]="$_[$i]";
			#print "$samble[$i+1]\t";
		}
	}
	my @headinfo;
	my $tp;
	
	#chomp($_[0]);
	
	for (my $i = 0; $i < length($_[0]); $i++)
	{
		$tp = substr($_[0],$i,1);
		if ($tp eq $samble[0] || $tp eq $samble[1] || $tp eq $_[1])
		{
			if ($test == 0){
				$li++;}
			$test =1;
		}
		else{
			my $canadd = 1;
	    	my $l = @_;
	    	for(my $j=1; $j<$l; $j++){
	    		if($tp eq $_[$j]){
	    	    	if ($test == 0){
	    	            $li++;}
	    	    		$test =1; 
	    	    		$canadd = 0;
	    	    		last;               			
	    		}
	    	}
	    	if ($canadd == 1){
	    		$headinfo[$li] = $headinfo[$li].$tp;
	    		$test = 0;
	    	}
		}       
	}
	return @headinfo;
}
