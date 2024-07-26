#get input information
use Getopt::Long;
#use Switch;
#use Math::GSL::Statistics qw /:all/;
use Math::Complex;
use Math::BigFloat;
use Data::Dumper qw(Dumper);

if(!(-e $ARGV[0])){
	&Usage(),exit(1);
}

$ARGV_len=@ARGV;
@ARGV_info;
for(my $i=1; $i<$ARGV_len; $i++){
	@ARGV_info[$i-1]=$ARGV[$i];
}
print "use @ARGV_info\t\n";	

my $lineIN='';  #input file read
my $outf=$ARGV[0].'.out';

open(OUT,'>',"$outf");
open(IN,'<',"$ARGV[0]");

my %db;
my @dba;
my @dbcount;
my $dblevel;
#read input file

#Do What?
#-----A--------------	
if($ARGV[1] eq "AddBlankToLine"){  #Add Blank To Line for align [AddBlankToLine]
	&AddBlankToLine(@ARGV_info);
}
if($ARGV[1] eq "Add_ID_for_vcf"){  #Add ID for vcf [Add_ID_for_vcf]
	&Add_ID_for_vcf(@ARGV_info);
}
elsif($ARGV[1] eq "AddLogforPl"){  #Add Log for perl [AddLogforPl]
	&AddLogforPl(@ARGV_info);
}
elsif($ARGV[1] eq "AddsubType"){  #add Sub Type of Name [AddsubType all w9 w8 r5 w10 r7]
	&addsubtype(@ARGV_info);	
}
elsif($ARGV[1] eq "AddLosstoLine"){ #add loss To Line (eg.fiste Line Is 1, 2, 5 Add the loss number 3, 4; Line Value )[AddLosstoLine 10 0]
	&AddLosstoLine(@ARGV_info);	
}		
elsif($ARGV[1] eq "Add0BehandNumber"){ #add 0 behand number (row Add To ? number of sites) [Add0BehandNumber 2 3]
	&add0behandnumber(@ARGV_info);
}	
elsif($ARGV[1] eq "Addlisttoinfile"){ #add list To infile (start digital frount after) [Addlisttoinfile 1 6 T n]
	&Addlisttoinfile(@ARGV_info);	
}	
elsif($ARGV[1] eq "AddfasTypeforRepeatMasker"){ #Add fas Type for RepeatMasker [AddfasTypeforRepeatMasker]
	&AddfasTypeforRepeatMasker(@ARGV_info);					
}
elsif($ARGV[1] eq "Addlisttoeachtype"){ #add list For each type (row start) only For sorted [Addlisttoeachtype 0 0]
	&Addlisttoeachtype(@ARGV_info);					
}
elsif($ARGV[1] eq "AddlistandReplacetoDuplicateID"){ #AddlistandReplacetoDuplicateID [AddlistandReplacetoDuplicateID]
	&AddlistandReplacetoDuplicateID(@ARGV_info);					
}
elsif($ARGV[1] eq "AddLineandMemberASymble"){ #add each line a symble for each memble (symble 4 digit)[AddLineandMemberASymble c 3]
	&AddLineandMemberASymble(@ARGV_info);					
}
elsif($ARGV[1] eq "AddLineASymble"){ #add each line a symble for all memble (symble 4 digit)[AddLineASymble CRTE 4]
	&AddLineASymble(@ARGV_info);					
}
elsif($ARGV[1] eq "AddLineSN"){ #add each line a SN [AddLineSN 1 BAR]
	&AddLineSN(@ARGV_info);					
}
elsif($ARGV[1] eq "AddLineSameSymble"){ #add each line same symble for all memble (f for head) (b for behand) ([AddLineSameSymble CRTE b]
	&AddLineSameSymble(@ARGV_info);					
}
elsif($ARGV[1] eq "AddNewClass"){ #Add New Class for blastclustal result [AddNewClass PGF -]
	&AddNewClass(@ARGV_info);					
}
elsif($ARGV[1] eq "AddNucleotide2fas"){ #Add Nucleotide to fas [AddNucleotide2fas GGG]
	&AddNucleotide2fas(@ARGV_info);					
}
elsif($ARGV[1] eq "AlignedFasFilter"){ #Aligned Fa sGet Percent (para defort 50) [AlignedFasFilter]
	&AlignedFasFilter(@ARGV_info);					
}
elsif($ARGV[1] eq "AlignedFasGetPercent"){ #Aligned Fas Get Percent to first [AlignedFasGetPercent]
	&AlignedFasGetPercent(@ARGV_info);					
}

#-----B--------------	
elsif($ARGV[1] eq "BatchAlntoFasbyList"){ #Batch Combin File bu Row Name [BatchAlntoFasbyList]
	&BatchAlntoFasbyList(@ARGV_info);
}
elsif($ARGV[1] eq "BatchAverageValueCount"){ #Batch Average Value Count (1-1 1-2) [BatchAverageValueCount]
	&BatchAverageValueCount(@ARGV_info);
}
elsif($ARGV[1] eq "BatchCombinFilebyRowName"){ #Batch Combin File bu Row Name [BatchCombinFilebyRowName]
	&BatchCombinFilebyRowName(@ARGV_info);
}
elsif($ARGV[1] eq "BatchComputeMedian"){ #Batch Compute Median [BatchComputeMedian]
	&BatchComputeMedian(@ARGV_info);
}
elsif($ARGV[1] eq "BatchClustalWbyList"){ #Batch ClustalW by list [BatchClustalwbyList]
	&BatchClustalWbyList(@ARGV_info);
}
elsif($ARGV[1] eq "Batch_PAML_codmel_mlc_Read"){ #Batch_PAML_codmel_mlc_Read £¨M0 branch model_A site£© [Batch_PAML_codmel_mlc_Read]
	&Batch_PAML_codmel_mlc_Read(@ARGV_info);
}
elsif($ARGV[1] eq "BatchTwoRowCompute"){ #batch compute row1 operator row2 ( 2 - 3 : 4 - 5 : 6 - 7) (+ - mult div) [BatchTwoRowCompute 2 - 3]
	&BatchTwoRowCompute(@ARGV_info);
}	
elsif($ARGV[1] eq "BatchMutiRowCompute"){ #batch compute muti row operator ( + 2 3 4 5 6 7) [BatchMutiRowCompute + 2 3 4 5 6 7 : + 4 5 6]
	&BatchMutiRowCompute(@ARGV_info);
}	
elsif($ARGV[1] eq "BatchSimpleCompute"){ #batch compute row operator number (+ - mult div) [BatchSimpleCompute 0 - 1] 
	&BatchSimpleCompute(@ARGV_info);
}		
elsif($ARGV[1] eq "BesInfoGet"){ #BES Get info from .fpc file [BesInfoGet 1]
	&GetBesInfo(@ARGV_info);
}	
elsif($ARGV[1] eq "BESctgRegulate"){ #BES Regulate ctg For CtgLink [BESctgRegulate 1]
	&RegulateBESctg(@ARGV_info);	
}	
elsif($ARGV[1] eq "BismarkExtractorOutputSlideCount"){  #Bismark Extractor Output Slide windows Count (window_size step_size)[BismarkExtractorOutputSlideCount 10000 1000]
	&BismarkExtractorOutputSlideCount(@ARGV_info);
}
elsif($ARGV[1] eq "BismarkExtractorOutCount"){  #only small genome can use:Bismark Extractor Output Count [BismarkExtractorOutCount]
	&BismarkExtractorOutCount(@ARGV_info);
}
elsif($ARGV[1] eq "BlatSearchFliter"){ #blat search fliter (range length_times) [BlatSearchFliter 0.8 2]
	&BlatSearchFliter(@ARGV_info);
}	
elsif($ARGV[1] eq "BlatFindPEmatch"){ #blat result Find Pair End match(region, file list Line must Add a useless record) [BlatFindPEmatch 200000 20]
	&FindPEmatch(@ARGV_info);
}	
elsif($ARGV[1] eq "Blat_2_loc"){ #blat result change to local >gengID|1|chrom:st..end [Blat_2_loc]
	&Blat_2_loc(@ARGV_info);
}
elsif($ARGV[1] eq "BlastGetFirstHit"){ #blast result get first hit [BlastGetFirstHit]
	&BlastGetFirstHit(@ARGV_info);
}
elsif($ARGV[1] eq "BlastGetDirect"){ #blast result get seq direct [BlastGetDirect]
	&BlastGetDirect(@ARGV_info);
}	
elsif($ARGV[1] eq "BlastHitLengthCount"){ #blast result musted cleaned, combin fregment information [BlastHitLengthCount]
	&BlastHitLengthCount(@ARGV_info);
}
elsif($ARGV[1] eq "BlastDeleteRedundance"){ #blast result delete Redundance hit [BlastDeleteRedundance]
	&BlastDeleteRedundance(@ARGV_info);
}
elsif($ARGV[1] eq "BlastResultCombin"){ #Blat Result Combin [BlastResultCombin]
	&BlastResultCombin(@ARGV_info);
}
elsif($ARGV[1] eq "BlastnChangeTargettoLocal"){ #Blastn Change to Local, add loss[BlastnChangeTargettoLocal]
	&BlastnChangeTargettoLocal(@ARGV_info);
}
elsif($ARGV[1] eq "BlastnChangeQuerrytoLocal"){ #Blastn Change to Local, add loss (p 20)[BlastnChangeQuerrytoLocal 20]
	&BlastnChangeQuerrytoLocal(@ARGV_info);
}
elsif($ARGV[1] eq "BreakSeqbyNNN"){ #[BreakSeqbyNNN]  ##
	&BreakSeqbyNNN(@ARGV_info);	
}
elsif($ARGV[1] eq "BreakSeqbyNNN_loc"){ #[BreakSeqbyNNN_loc]  ##give start and end of ctg in .loc file
	&BreakSeqbyNNN_loc(@ARGV_info);	
}
elsif($ARGV[1] eq "Bsmap_2_methylKit"){ #[Bsmap_2_methylKit]  
	&Bsmap_2_methylKit(@ARGV_info);	
}
elsif($ARGV[1] eq "bismark_2_methylKit"){ #[bismark_2_methylKit]  
	&bismark_2_methylKit(@ARGV_info);	
}
#-----C--------------	
elsif($ARGV[1] eq "CDHitClstrInfo"){ #get CD-Hit Clstr combined Info (bca.clstr file) [CDHitClstrInfo]
	&CDHitClstrInfo(@ARGV_info);	
}
elsif($ARGV[1] eq "Change_fas_to_complement"){ #Change_fas_to_complement [Change_fas_to_complement]
	&Change_fas_to_complement(@ARGV_info);	
}
elsif($ARGV[1] eq "ChangeSymbleToRef"){ #Change Symble To Ref [ChangeSymbleToRef 3 2 .]
	&ChangeSymbleToRef(@ARGV_info);	
}
elsif($ARGV[1] eq "ChangeRowbyCompair"){ #Change Row byCompair[ChangeRowbyCompair 8 9]
	&ChangeRowbyCompair(@ARGV_info);	
}
elsif($ARGV[1] eq "ChangeLTRinfoTolocinfo"){ #Change LTR info To loc info[ChangeLTRinfoTolocinfo]
	&ChangeLTRinfoTolocinfo(@ARGV_info);	
}	
elsif($ARGV[1] eq "ChangeFastaN"){ #Change Fasta N to -[ChangeFastaN N *]
	&ChangeFastaN(@ARGV_info);	
}
elsif($ARGV[1] eq "ChangeGffouttoLoc"){ #Change Gff out file to Loc file £¨CDS£©(gene) [ChangeGffouttoLoc CDS]
	&ChangeGffouttoLoc(@ARGV_info);	
}
elsif($ARGV[1] eq "ChangeStringtoCode"){ #Change String to Code (jap-w8-w9) [ChangeStringtoCode]
	&ChangeStringtoCode(@ARGV_info);	
}
elsif($ARGV[1] eq "ChangeSoapsnpConsensus_toFas"){ #Change Soapsnp Consensus to Fas [ChangeSoapsnpConsensus_toFas]
	&ChangeSoapsnpConsensus_toFas(@ARGV_info);	
}
elsif($ARGV[1] eq "Change_Soapsnp_Consensus_toFas"){ #Change Soapsnp Consensus to Fas (one by one) [Change_Soapsnp_Consensus_toFas]
	&Change_Soapsnp_Consensus_toFas(@ARGV_info);	
}
elsif($ARGV[1] eq "ChechIndelType"){ #Change Fasta N to -[ChechIndelType]
	&ChechIndelType(@ARGV_info);	
}
elsif($ARGV[1] eq "ChechPatternDiffer"){ #Chech Pattern Differ, find and change[ChechPatternDiffer]
	&ChechPatternDiffer(@ARGV_info);	
}
elsif($ARGV[1] eq "ChechSameFamilybyRow"){ #Chech Same Family by Row (row1 row2 cut_symble)[ChechSameFamilybyRow 0 1 _]
	&ChechSameFamilybyRow(@ARGV_info);	
}
elsif($ARGV[1] eq "CheckStatebyValue"){ #Check State (1>= 0<)by Value [CheckStatebyValue 0]
	&CheckStatebyValue(@ARGV_info);	
}
elsif($ARGV[1] eq "CombinFasbyID"){ #Combin Fas by ID [CombinFasbyID]
	&CombinFasbyID(@ARGV_info);	
}		
elsif($ARGV[1] eq "CombinbyLine"){ #combin lines by number [CombinbyLine 10]
	&CombinbyLine(@ARGV_info);	
}	
elsif($ARGV[1] eq "CombinLinestoOne"){ #combin lines To one Line [CombinLinestoOne 2]
	&CombinLinestoOne(@ARGV_info);	
}	
elsif($ARGV[1] eq "Compairstringrow"){ #compair tow String row accord operator,  eq  ne [Compairstringrow 5 eq 10]
	&compairstringrow(@ARGV_info);	
}	
elsif($ARGV[1] eq "compairintrow"){ #compair Int row accord operator(lt eq gt) [compairintrow 14 lt 9 10 11 12 13]
	&compairintrow(@ARGV_info);	
}
elsif($ARGV[1] eq "CombinCorrespondenceElements"){ #Combin Correspondence Elements [CombinCorrespondenceElements]
	&CombinCorrespondenceElements(@ARGV_info);
}
elsif($ARGV[1] eq "CombinTypetoLines"){ #Combin same Type to one Line[CombinTypetoLines _]
	&CombinTypetoLines(@ARGV_info);
}
elsif($ARGV[1] eq "Combinloc"){ #combin location by combin region (infile: Name start End ) [Combinloc 100]
	&Combinloc(@ARGV_info);	
}	
elsif($ARGV[1] eq "CombinAlternativeSplicingCDSloc"){ #Combin Alternative Splicing CDS loc [CombinAlternativeSplicingCDSloc]
	&CombinAlternativeSplicingCDSloc(@ARGV_info);	
}	
elsif($ARGV[1] eq "Combinclusterrecord"){ #combin cluster record (the Input must be corresponding) [Combinclusterrecord 1]
	&Combinclusterrecord(@ARGV_info);
}
elsif($ARGV[1] eq "CombinRepeatmaskerResult"){ #regulated repeatmasker result combin (combin results from littletool) [CombinRepeatmaskerResult 1]
	&CombinRepeatmaskerResult(@ARGV_info);	
}
elsif($ARGV[1] eq "Combinmutirow"){ #combin two row by Name [Combinmutirow 1 13 14]
	&Combinmutirow(@ARGV_info);		
}
elsif($ARGV[1] eq "CombinBlastclustalforBothLTR"){ #Combin Blastclustal result for Both LTR input [CombinBlastclustalforBothLTR]
	&CombinBlastclustalforBothLTR(@ARGV_info);		
}
elsif($ARGV[1] eq "CombinCDHitClstrToClass"){ #Combin CD-Hit Clstr To Blastclustal Class type [CombinCDHitClstrToClass]
	&CombinCDHitClstrToClass(@ARGV_info);		
}
elsif($ARGV[1] eq "CombinFastaSeqtoOne"){ #combin all seq in a fasta file to one, add NNN between seq [CombinFastaSeqtoOne 500 N]
	&CombinFastaSeqtoOne(@ARGV_info);		
}	
elsif($ARGV[1] eq "CombinFilesbyList"){ #Combin Files by List [CombinFilesbyList]
	&CombinFilesbyList(@ARGV_info);		
}	
elsif($ARGV[1] eq "CombinlBlastResult"){ #combin blast result (combin two record in same chromsome And Within cutoff distances)(Name chrom start End) [CombinlBlastResult 1]
	&CombinlBlastResult(@ARGV_info);	
}	
elsif($ARGV[1] eq "CombinltoBestQuerry"){ #combin same blast Loc To best querry (If some records was in same location,  keep the record With have largest e value) [CombinltoBestQuerry 1]
	&CombinltoBestQuerry(@ARGV_info);	
}	
elsif($ARGV[1] eq "CombinltoLargeSubType"){ #Combinl to Large Sub Type (type subtype value) [CombinltoLargeSubType 0 1 2]
	&CombinltoLargeSubType(@ARGV_info);	
}
elsif($ARGV[1] eq "CombinlSameType"){ #Combinl Same Type, subtype combin to a string (type) (other row combin) [CombinlSameType 0]
	&CombinlSameType(@ARGV_info);	
}
elsif($ARGV[1] eq "CompareValuestoCutoff"){ #Compare each Values to Cutoff to check number [CompareValuestoCutoff gt 0.66]
	&CompareValuestoCutoff(@ARGV_info);	
}
elsif($ARGV[1] eq "ConnectfasbyNNN"){ #Connect fasby NNN (50) [ConnectfasbyNNN 50]
	&ConnectfasbyNNN(@ARGV_info);	 
}
elsif($ARGV[1] eq "ConnectfasbyName"){ #Connect fasby Name, if need add ATGCATGCATGCATGC (>Name|detail) [ConnectfasbyName ATGCATGCATGCATGC]
	&ConnectfasbyName(@ARGV_info);	 
}
elsif($ARGV[1] eq "ConnectfasbyName_sorted"){ #Connect fasby same Name sorted [ConnectfasbyName_sorted]
	&ConnectfasbyName_sorted(@ARGV_info);	 
}
elsif($ARGV[1] eq "ConnectSortedfasbyName"){ #Connect sorted fasby Name [ConnectSortedfasbyName]
	&ConnectSortedfasbyName(@ARGV_info);	 
}
elsif($ARGV[1] eq "ConnectSortedLocal"){ #Connect sorted local (ID st ed) conn para [ConnectSortedLocal 200]
	&ConnectSortedLocal(@ARGV_info);	 
}
elsif($ARGV[1] eq "ConnectExcelEnter"){ #Connect Excel Enter [ConnectExcelEnter]
	&ConnectExcelEnter(@ARGV_info);	 
}
elsif($ARGV[1] eq "count_3_piont_trend"){ #count 3 piont trend [count_3_piont_trend]
	&count_3_piont_trend(@ARGV_info);
}
elsif($ARGV[1] eq "count01Matrix"){ #count 0 1 Matrix [count01Matrix]
	&count01Matrix(@ARGV_info);
}	
elsif($ARGV[1] eq "count_col_average_value"){ #count_col_average_avlue [count_col_average_value 2]
	&count_col_average_value(@ARGV_info);
}	
elsif($ARGV[1] eq "countFileLength"){ #count File length [countFileLength 1]
	&countFileLength(@ARGV_info);
}	
elsif($ARGV[1] eq "countGCcontant"){ #count GC contant [countGCcontant 1]
	&countGCcontant(@ARGV_info);
}
elsif($ARGV[1] eq "count_seq_Symble_No"){ #count seq Symble_No [count_seq_Symble_No]
	&count_seq_Symble_No(@ARGV_info);
}	
elsif($ARGV[1] eq "CountlineSymbleNo"){ #Count each line Symble No [CountlineSymbleNo RL (g) (c)]
	&CountlineSymbleNo(@ARGV_info);
}
elsif($ARGV[1] eq "Countline_value_No"){ #Count each line value No [Countline_value_No lt 0]
	&Countline_value_No(@ARGV_info);
}
elsif($ARGV[1] eq "CountlineLength"){ #Count each line tlength [CountlineLength]
	&CountlineLength(@ARGV_info);
}
elsif($ARGV[1] eq "CountlineTypeNo"){ #Count each line type No and list line types [CountlineTypeNo _]
	&CountlineTypeNo(@ARGV_info);
}
elsif($ARGV[1] eq "CountlineElementsNo"){ #Count each line Elements No [CountlineElementsNo]
	&CountlineElementsNo(@ARGV_info);
}
elsif($ARGV[1] eq "CountlineFamilyNo"){ #Count each line families No and list line families [CountlineFamilyNo]
	&CountlineFamilyNo(@ARGV_info);
}
elsif($ARGV[1] eq "CountSoapCoverageAvgDepth"){ #Count Soap Coverage Avg Depth [CountSoapCoverageAvgDepth]
	&CountSoapCoverageAvgDepth(@ARGV_info);	
}	
elsif($ARGV[1] eq "CountSoapCoverageWindowsAvgDepth"){ #Count Soap Coverage windows avg Depth [CountSoapCoverageWindowsAvgDepth 1000]
	&CountSoapCoverageWindowsAvgDepth(@ARGV_info);	
}	
elsif($ARGV[1] eq "Counttypeno"){ #count Type number by row  [Counttypeno 13]
	&Counttypeno(@ARGV_info);	
}
elsif($ARGV[1] eq "CountCombintypeno"){ #count 22 combin Type number by row  [CountCombintypeno 0 1]
	&CountCombintypeno(@ARGV_info);	
}
elsif($ARGV[1] eq "CountCombintypeValue"){ #count 22 combin Type value by row  [CountCombintypeValue 0 1 2]
	&CountCombintypeValue(@ARGV_info);	
}
elsif($ARGV[1] eq "CounttypeAllinfoNo"){ #Count type All info, then cout each sub type No  [CounttypeAllinfoNo 0]
	&CounttypeAllinfoNo(@ARGV_info);	
}	
elsif($ARGV[1] eq "Counttypecontinuationno"){ #Count type continuation no  [Counttypecontinuationno 0]
	&Counttypecontinuationno(@ARGV_info);	
}
elsif($ARGV[1] eq "CounttypePercentRegionValue"){ #count Type everage value, sorted (type value_1 value_2 value_3 ....)  [CounttypePercentRegionValue]
	&CounttypePercentRegionValue(@ARGV_info);	
}
elsif($ARGV[1] eq "CounttypeAverageValue"){ #count Type everage value, sorted (type value_1 value_2 value_3 ....)  [CounttypeAverageValue 0 1 2 3 4]
	&CounttypeAverageValue(@ARGV_info);	
}
elsif($ARGV[1] eq "CounttypeTotalValue"){ #count Type total value (type v1 v2 v3) [CounttypeTotalValue 0 1 2 3]
	&CounttypeTotalValue(@ARGV_info);	
}
elsif($ARGV[1] eq "Counttype_nonredundant_Length"){ #Count type_nonredundant_Length, input must sorted (type st ed) [Counttype_nonredundant_Length 0 1 2]
	&Counttype_nonredundant_Length(@ARGV_info);	
}
elsif($ARGV[1] eq "Countsubtypeno"){ #count Sub Type number by row (Type subtype) [Countsubtypeno 0 2]
	&Countsubtypeno(@ARGV_info);	
}
elsif($ARGV[1] eq "CounttypeRegionNumberAndTotalvalue"){ #count Type Region Number and total value (region size) [CounttypeRegionNumberAndTotalvalue 1]
	&CounttypeRegionNumberAndTotalvalue(@ARGV_info);	
}
elsif($ARGV[1] eq "Countsubloc"){ #count Type Loc no by row (Type range) (same range As same Loc) [Countsubloc 5 20]
	&Countsubloc(@ARGV_info);		
}
elsif($ARGV[1] eq "CountRegionNumber"){ #count region number (need change) [CountRegionNumber 100]
	&CountRegionNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "CountRegionLength"){ #count region number And length [CountRegionLength 10000]
	&CountRegionLength(@ARGV_info);	
}
elsif($ARGV[1] eq "CountRegionDetail"){ #count >para region number And length average etc (1000 10000 100000) [CountRegionDetail 0.5 1 2 5 20 50 500]
	&CountRegionDetail(@ARGV_info);	
}
elsif($ARGV[1] eq "CountfixedRegionNumber"){ #count region number (para 10000, fixed region) (region row1 row2 row3) [CountfixedRegionNumber 1000]
	&CountfixedRegionNumber(@ARGV_info);	
}	
elsif($ARGV[1] eq "CountSourtedRegionAvgValue"){ #count sorted region average value (inpute: 0-region 1-value,input must sorted) (para 10000, fixed region)  [CountSourtedRegionAvgValue 1000]
	&CountSourtedRegionAvgValue(@ARGV_info);	
}	
elsif($ARGV[1] eq "CountfixedRegionAvgValue"){ #count region average value (inpute: 0-region 1-value, start from 0) (para: 1-region 10000, 2-max 100000 )  [CountfixedRegionAvgValue 10000 500000]
	&CountfixedRegionAvgValue(@ARGV_info);	
}	
elsif($ARGV[1] eq "CountTypefixedRegionNumber"){ #count type region number (inpute: 0-type 1-region,input must sorted) (para 10000, fixed region)  [CountTypefixedRegionNumber 1000]
	&CountTypefixedRegionNumber(@ARGV_info);	
}	
elsif($ARGV[1] eq "CountchromRegionCoverage"){ #count chrom region  coverage (inpute: 0-chr 1-st 2-ed, start from 1)£¨not need sort£© (para: 1-regin )  [CountchromRegionCoverage 5000000]
	&CountchromRegionCoverage(@ARGV_info);	
}
elsif($ARGV[1] eq "CountchromRegionTotalValue"){ #count chrom region total value (inpute: 0-chr 1-loc 2-value, start from 1)£¨not need sort£© (para: 1-regin )  [CountchromRegionTotalValue 5000000]
	&CountchromRegionTotalValue(@ARGV_info);	
}
elsif($ARGV[1] eq "Counthitspan"){ # [Counthitspan ?]
	&counthitspan(@ARGV_info);
}
elsif($ARGV[1] eq "Counttypelength"){ #count Type length by Type row & length row [Counttypelength 14 10]
	&counttypelength(@ARGV_info);	
}	
elsif($ARGV[1] eq "countseqlength"){ #Count sequences length  [countseqlength 1]
	&countseqlength(@ARGV_info);	
}
elsif($ARGV[1] eq "countfastqtotallen_GC"){ #Count sequences length  [countfastqtotallen_GC]
	&countfastqtotallen_GC(@ARGV_info);	
}
elsif($ARGV[1] eq "countFastqLength"){ #Count sequences length (if had d, out put length detai;) [countFastqLength d]
	&countFastqLength(@ARGV_info);	
}
elsif($ARGV[1] eq "Countstringlength"){ #count String length [Countstringlength 1]
	&countstringlength(@ARGV_info);	
}
elsif($ARGV[1] eq "countAccumulationValue"){ #count Accumulation Value  [countAccumulationValue 1]
	&countAccumulationValue(@ARGV_info);	
}	
elsif($ARGV[1] eq "CPCgetCodingList"){ #CPC get CPCgetCodingList from one of forward and reverse seq score > 0 [CPCgetCodingList]
	&CPCgetCodingList(@ARGV_info);	
}
elsif($ARGV[1] eq "csv_2_tsv"){ #csv change to tsv [csv_2_tsv]
	&csv_2_tsv(@ARGV_info);	
}
elsif($ARGV[1] eq "CutfastqtoLength"){ #Count sequences length  [CutfastqtoLength 0 145]
	&CutfastqtoLength(@ARGV_info);	
}	
elsif($ARGV[1] eq "CutAlignFasbySeedSymble"){ #Cut aligned fas head and tail by Seed Symble  [CutAlignFasbySeedSymble AF_3357]
	&CutAlignFasbySeedSymble(@ARGV_info);	
}
elsif($ARGV[1] eq "CutAlignFasGapForProteinCoding"){ #Cut aligned fas gap for protein coding sequence at 3 code [CutAlignFasGapForProteinCoding]
	&CutAlignFasGapForProteinCoding(@ARGV_info);	
}
elsif($ARGV[1] eq "CutSequence_by_length"){ #Cut fasta by length  [CutSequence_by_length 100]
	&CutSequence_by_length(@ARGV_info);	
}	
elsif($ARGV[1] eq "CutfastqfromStart"){ #Cut fastq from Start  [CutfastqfromStart 20]
	&CutfastqfromStart(@ARGV_info);	
}	
elsif($ARGV[1] eq "CutfastqfromEnd"){ #Cut fastq from End  [CutfastqfromEnd 3]
	&CutfastqfromEnd(@ARGV_info);	
}	
elsif($ARGV[1] eq "CutfastafromEnd"){ #Cut fasta from End  [CutfastafromEnd 3]
	&CutfastafromEnd(@ARGV_info);	
}
elsif($ARGV[1] eq "CutLenToCode"){ #Cut Len To Code [CutLenToCode]
	&CutLenToCode(@ARGV_info);	
}
elsif($ARGV[1] eq "Cut_to_Max_value"){ #change value higher to max to max value [Cut_to_Max_value 0 40]
	&Cut_to_Max_value(@ARGV_info);	
}
elsif($ARGV[1] eq "CutSortedDataToMeanValue"){ #Cut Sorted Data To Mean Value (meav_value row) [CutSortedDataToMeanValue 4177 1]
	&CutSortedDataToMeanValue(@ARGV_info);	
}	
#-----D--------------	
elsif($ARGV[1] eq "Deletesameline"){ #delete same Line (sorted)  [Deletesameline 1]
	&Deletesameline(@ARGV_info);
}	
elsif($ARGV[1] eq "Deletenestloc"){ #delete nest Loc in same Type (Type	chrom start End)(Input seq must sorted) [Deletenestloc 1 1 2 3]
	&deletenestloc(@ARGV_info);	
}
elsif($ARGV[1] eq "Deleteshortline"){ #delete short Line (row number,  cut short Line)  [Deleteshortline 1]
	&Deleteshortline(@ARGV_info);
}
elsif($ARGV[1] eq "DEseq2_input_from_HTseq"){ #DEseq2_input_regulate from_HTseq by list  [DEseq2_input_from_HTseq]
	&DEseq2_input_from_HTseq(@ARGV_info);
}	
elsif($ARGV[1] eq "DEseq_2_Goseq_input"){ #DEseq output change to Goseq input)  [DEseq_2_Goseq_input 0.05]
	&DEseq_2_Goseq_input(@ARGV_info);
}		
elsif($ARGV[1] eq "DndToTrees"){ #Clustalx dnd change to trees for paml (name must same length) [DndToTrees]
	&DndToTrees(@ARGV_info);	
}
elsif($ARGV[1] eq "DrawSaturationCurve"){ #Draw Saturation Curve us all state [DrawSaturationCurve]
	&DrawSaturationCurve(@ARGV_info);	
}
#-----E--------------	
elsif($ARGV[1] eq "Embl_2_geneloc"){ #Embl_2_geneloc (doing) [Embl_2_geneloc]
	&Embl_2_geneloc(@ARGV_info);	
}
elsif($ARGV[1] eq "EMMAX_ps_2_Manhattan_input"){ #EMMAX_ps_2_Manhattan_input [EMMAX_ps_2_Manhattan_input]
	&EMMAX_ps_2_Manhattan_input(@ARGV_info);	
}
elsif($ARGV[1] eq "Extract_KEGG_compound_info"){ #Extract_KEGG_compound_info  [Extract_KEGG_compound_info]
	&Extract_KEGG_compound_info(@ARGV_info);	
}
#-----F--------------	
elsif($ARGV[1] eq "FasToPhy"){ #Fasta change to phy for paml (name must same length) cut to 10 [FasToPhy 10]
	&FasToPhy(@ARGV_info);	
}
elsif($ARGV[1] eq "FastqToFasta"){ #Fastq change to Fasta [FastqToFasta]
	&FastqToFasta(@ARGV_info);	
}
elsif($ARGV[1] eq "Findnestloc"){ #find nest Loc (nest: distance With in 400)£¨Input file: Name chrom START End£© [Findnestloc 400]
	&Findnestloc(@ARGV_info);	
}
elsif($ARGV[1] eq "Find6alignbigestRegion"){ #Find 6align local bigest Region, need no gap local [Find6alignbigestRegion]
	&Find6alignbigestRegion(@ARGV_info);
}
elsif($ARGV[1] eq "Find6alignhit"){ #Find and seperate 6align hit in 6 species [Find6alignhit]
	&Find6alignhit(@ARGV_info);
}
elsif($ARGV[1] eq "FindbestORF"){ #find the best translate ORF And Count Stop coden [FindbestORF 1]
	&FindbestORF(@ARGV_info);
}
elsif($ARGV[1] eq "FindAllORF"){ # List all orf [FindAllORF]
	&FindAllORF(@ARGV_info);
}
elsif($ARGV[1] eq "FindlinesOverlap"){ #Find lines Overlap to frevious one [FindlinesOverlap 1]
	&FindlinesOverlap(@ARGV_info);
}
elsif($ARGV[1] eq "FindlocalOverlap"){ #Find tow local Overlap to frevious one (chr1 st1 ed1 chr2 st2 ed2) [FindlocalOverlap 1 3 4 6 8 9]
	&FindlocalOverlap(@ARGV_info);
}
elsif($ARGV[1] eq "FindLonglinebyCutoff"){ #Find Long line by Cutoff [FindLonglinebyCutoff 200]
	&FindLonglinebyCutoff(@ARGV_info);
}
elsif($ARGV[1] eq "FindPairdata"){ #Find Pair  data [FindPairdata]
	&FindPairdata(@ARGV_info);
}
elsif($ARGV[1] eq "Find_NNN_MoreThanCutoff"){ #Find continue NNN More Than Cut off [Find_NNN_MoreThanCutoff 10]
	&Find_NNN_MoreThanCutoff(@ARGV_info);
}
elsif($ARGV[1] eq "FindSameModule"){ #find the Right hit from module (module module1 info module2 info module3 info ... ...)(numbrt means info length, e.g. 2 means 2 group info, one group info use module's n value)  [FindSameModule 5]
	&FindSameModule(@ARGV_info);	
}
elsif($ARGV[1] eq "Floot_2_int_up"){ # floot to int, (0.1 -1) to 1  (para line) [Floot_2_int_up 12]
	&Floot_2_int_up(@ARGV_info);	
}	
#-----G--------------
elsif($ARGV[1] eq "gb_gbk_2_genelocal"){ #gb_gbk_2_genelocal [gb_gbk_2_genelocal]
	&gb_gbk_2_genelocal(@ARGV_info);	
}	
elsif($ARGV[1] eq "gb_gbk_2_pro_fas"){ #gb_gbk_2_pro_fas [gb_gbk_2_pro_fas]
	&gb_gbk_2_pro_fas(@ARGV_info);	
}
elsif($ARGV[1] eq "gff_trans_2_pro_fas"){ #gff_trans_2_pro_fas [gff_trans_2_pro_fas]
	&gff_trans_2_pro_fas(@ARGV_info);	
}	
elsif($ARGV[1] eq "gff_2_gtf_for_tophat"){ #gff_2_gtf_for_tophat [gff_2_gtf_for_tophat]
	&gff_2_gtf_for_tophat(@ARGV_info);	
}	
elsif($ARGV[1] eq "getBlastRelationship"){ #get Blast Relationship for LTR classify [getBlastRelationship]
	&getBlastRelationship(@ARGV_info);	
}	
elsif($ARGV[1] eq "getGOannofromInterproFile"){ #get row from number And info behand symble (all.interpro) [getGOannofromInterproFile]
	&getGOannofromInterproFile(@ARGV_info);	
}
elsif($ARGV[1] eq "GetbigSmallone"){ #get small/big one from rows (big / small) [GetbigSmallone big]
	&GetbigSmallone(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGapInfo"){ #get gap loc and length from seq [GetGapInfo N]
	&GetGapInfo(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGeneCodebyTECode"){ #Get Gene Code by TE Code (Gene must sorted) [GetGeneCodebyTECode]
	&GetGeneCodebyTECode(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGeneExpfromSoftFile"){ #Get Gene Expression Data from Soft File download from NCBI [GetGeneExpfromSoftFile]
	&GetGeneExpfromSoftFile(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGeneExpfromSoftFile"){ #Get Gene Expression Data from Soft File download from NCBI [GetGeneExpfromSoftFile]
	&GetGeneExpfromSoftFile(@ARGV_info);	
}
elsif($ARGV[1] eq "GetEnsemblGeneAnnCDS"){ #Get Ensembl Gene Annotation Info [GetEnsemblGeneAnnCDS]
	&GetEnsemblGeneAnnCDS(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGffanno_by_name"){ #Get Gff anno by name [GetGffanno_by_name ID Parent product Note]
	&GetGffanno_by_name(@ARGV_info);	
}	
elsif($ARGV[1] eq "GetGffanno_by_index"){ #Get Gff anno by index (index) [GetGffanno_by_index ID=cds- gene= Parent=rna-]
	&GetGffanno_by_index(@ARGV_info);	
}		
elsif($ARGV[1] eq "GetGff3GeneInfo_our"){ #Get Gff3 Gene Info [GetGff3GeneInfo_our]
	&GetGff3GeneInfo_our(@ARGV_info);	
}	
elsif($ARGV[1] eq "GetGff3Gene_exon_length_jv6"){ #Get Gff3 Gene exon length [GetGff3Gene_exon_length_jv6]
	&GetGff3Gene_exon_length_jv6(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGff3Gene_mRNA_CDS_name"){ #Get Gff3 Gene Info [GetGff3Gene_mRNA_CDS_name]
	&GetGff3Gene_mRNA_CDS_name(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGff3GeneInfo"){ #Get Gff3 Gene Info [GetGff3GeneInfo]
	&GetGff3GeneInfo(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGff3CDSInfo_jv7"){ #Get Gff3 Gene CDS Info [GetGff3CDSInfo_jv7]
	&GetGff3CDSInfo_jv7(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGff3out_to_gdb"){ #change Gff3 out to a file for gdb sqlit3 read [GetGff3out_to_gdb]
	&GetGff3out_to_gdb(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGff3GeneInfo_jv6"){ #Get Gff3 Gene Info [GetGff3GeneInfo_jv6]
	&GetGff3GeneInfo_jv6(@ARGV_info);	
}	
elsif($ARGV[1] eq "GetGff3GeneInfo_detail_our"){ #***both our and jv6 can use  Get Gff3 Gene Info  [GetGff3GeneInfo_detail_our]
	&GetGff3GeneInfo_detail_our(@ARGV_info);	
}
elsif($ARGV[1] eq "GetGff3GeneInfo_AlternativeSplicing_jv6"){ #Get Gff3 Gene Info [GetGff3GeneInfo_AlternativeSplicing_jv6]
	&GetGff3GeneInfo_AlternativeSplicing_jv6(@ARGV_info);	
}
elsif($ARGV[1] eq "GetMin"){ #GetMin [GetMin]
	&GetMin(@ARGV_info);	
}	
elsif($ARGV[1] eq "GetMax"){ #GetMin [GetMax]
	&GetMax(@ARGV_info);	
}	
elsif($ARGV[1] eq "GetMafAlignedTable"){ #Get lastz/last Maf Aligned Table [GetMafAlignedTable]
	&GetMafAlignedTable(@ARGV_info);	
}	
elsif($ARGV[1] eq "GetMafAlignedQualifiedRegion"){ #Get lastz/last Maf Aligned Qualified Region by region cut length  [GetMafAlignedQualifiedRegion 50]
	&GetMafAlignedQualifiedRegion(@ARGV_info);	
}
elsif($ARGV[1] eq "GetNRfirstname"){ #Get N Rfirst name  [GetNRfirstname]
	&GetNRfirstname(@ARGV_info);	
}	
elsif($ARGV[1] eq "getLinefromNumber"){ #get Line from Number (start how_many) [getLinefromNumber 2549401 100000 ]
	&getLinefromNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "getLineFromSymble"){ #get Line from symble  [getLineFromSymble symble ]
	&getLineFromSymble(@ARGV_info);	
}
elsif($ARGV[1] eq "getLineWithoutSymble"){ #get Line without symble [getLineWithoutSymble symble ]
	&getLineWithoutSymble(@ARGV_info);	
}
elsif($ARGV[1] eq "getLinfromName"){ #get Line from name (name row) [getLinfromName name 2 ]
	&getLinfromName(@ARGV_info);	
}
elsif($ARGV[1] eq "getLinWithoutName"){ #get Line without name (name row) [getLinWithoutName name 2 ]
	&getLinWithoutName(@ARGV_info);	
}
elsif($ARGV[1] eq "GetLineinEverageRegion"){  #Get lines in clase everage region (class value region) [GetLineinEverageRegion 0 2 3000]
	&GetLineinEverageRegion(@ARGV_info);	
}
elsif($ARGV[1] eq "getlocMapEverageDepth"){ #get loc Map Everage Depth from getlocMapDepth_lowMem result [getlocMapEverageDepth]
	&getlocMapEverageDepth(@ARGV_info);	
}	
elsif($ARGV[1] eq "getMaskedRegionLocalfromFasta"){ #ge tMasked Region Local from Fasta [getMaskedRegionLocalfromFasta]
	&getMaskedRegionLocalfromFasta(@ARGV_info);	
}
elsif($ARGV[1] eq "GetRowMatrixfromFilelist"){ #get row from number by file list, out put matrix ( key = 0 ) [GetRowMatrixfromFilelist 9]
	&GetRowMatrixfromFilelist(@ARGV_info);	
}
elsif($ARGV[1] eq "GetRowMatrixfromFilelist_sort_muti"){ #get row from number by file list, out put matrix,file must be same length and sorted ( key = 0 ) [GetRowMatrixfromFilelist_sort_muti 1 2]
	&GetRowMatrixfromFilelist_sort_muti(@ARGV_info);	
}
elsif($ARGV[1] eq "GetRowfromNumber"){ #get row from number,  can Add more than one row ) [GetRowfromNumber 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32]
	&GetRowfromNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "GetRowfromNumber_splite"){ #get row from number,  can Add more than one row (> | : ..) [GetRowfromNumber_splite 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32]
	&GetRowfromNumber_splite(@ARGV_info);	
}
elsif($ARGV[1] eq "GetRowWithoutNumber"){ #get row without number,  can Add more than one row ) [GetRowWithoutNumber 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32]
	&GetRowWithoutNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "GetfasNamefromNumber"){ #get fas name from number,  sequence no change) [GetfasNamefromNumber 1]
	&GetfasNamefromNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "GetRowtoLinesfromNumber"){ #get row from number,  can Add more than one row, add ^n between rows ) [GetRowtoLinesfromNumber 3 6]
	&GetRowtoLinesfromNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "getAllSeqAfterName"){ #get seq have both side hit [getAllSeqAfterName name]
	&getAllSeqAfterName(@ARGV_info);	
}
elsif($ARGV[1] eq "getseqfromsymble"){ #get seq from symble [getseqfromsymble MIKCC]
	&getseqfromsymble(@ARGV_info);	
}
elsif($ARGV[1] eq "getseqwithoutsymble"){ #get seq without symble [getseqwithoutsymble MIKCC]
	&getseqwithoutsymble(@ARGV_info);	
}
elsif($ARGV[1] eq "Getseqwithbothside"){ #get seq have both side hit [Getseqwithbothside 1]
	&getseqwithbothside(@ARGV_info);	
}
elsif($ARGV[1] eq "GetPairLTRbyDirect"){  #Get Pair LTR by Direct(name direct loc) [GetPairLTRbyDirect 13]
	&GetPairLTRbyDirect(@ARGV_info);	
}
elsif($ARGV[1] eq "GetFirstSubtype"){  #Get First Sub type (Sorted) [GetFirstSubtype 0 1]
	&GetFirstSubtype(@ARGV_info);	
}	
elsif($ARGV[1] eq "GetSequencefromEstdata"){ #Get Sequence from Est data file down load from NCBI [GetSequencefromEstdata]
	&GetSequencefromEstdata(@ARGV_info);	
}
elsif($ARGV[1] eq "GetStringAfterString"){  #Get String After String [GetStringAfterString]
	&GetStringAfterString(@ARGV_info);	
}	
elsif($ARGV[1] eq "Gettypefirstrow"){  #Get Type first row (sorted) [Gettypefirstrow 0 1]
	&Gettypefirstrow(@ARGV_info);	
}
elsif($ARGV[1] eq "GettypeToprows"){  #Get Type first row (sorted) (line 0 top 5) [GettypeToprows 0 5]
	&GettypeToprows(@ARGV_info);	
}
elsif($ARGV[1] eq "Gettypefirstsequence"){  #Get Type first sequenced (sorted) [Gettypefirstsequence]
	&Gettypefirstsequence(@ARGV_info);	
}
elsif($ARGV[1] eq "GetTypeLongestfas"){  #Get Type longest sequenced (sorted) [GetTypeLongestfas]
	&GetTypeLongestfas(@ARGV_info);	
}
elsif($ARGV[1] eq "GetUnigeneFromTrinityFasta"){  #Get Unigene From Trinity.Fasta[GetUnigeneFromTrinityFasta]
	&GetUnigeneFromTrinityFasta(@ARGV_info);	
}
elsif($ARGV[1] eq "Gff3_2_gtf"){  #change Gff3 to gtf [Gff3_2_gtf]
	&Gff3_2_gtf(@ARGV_info);	
}
elsif($ARGV[1] eq "Gtf_2_gene"){  #change Gtf to gene [Gtf_2_gene]
	&Gtf_2_gene(@ARGV_info);	
}
elsif($ARGV[1] eq "Gtf_2_join_loc"){  #change Gtf to Gtf_txt [Gtf_2_join_loc]
	&Gtf_2_join_loc(@ARGV_info);	
}
elsif($ARGV[1] eq "Gtf_txt_2_join_loc"){  #change Gtf.txt to joined gene local [Gtf_txt_2_join_loc]
	&Gtf_txt_2_join_loc(@ARGV_info);	
}
elsif($ARGV[1] eq "GiveEachLineMemberASymble"){  #Give Each Line Member A Symble (1, 2, 3...)[GiveEachLineMemberASymble 0]
	&GiveEachLineMemberASymble(@ARGV_info);	
}
#-----H--------------
elsif($ARGV[1] eq "HTseq_to_DESeq"){ #HTseq output change to DESeq input  [HTseq_to_DESeq]
	&HTseq_to_DESeq(@ARGV_info);
}
elsif($ARGV[1] eq "Hmp_2_NN"){ #hmp file change heterozygosis to NN [Hmp_2_NN]
	&Hmp_2_NN(@ARGV_info);
}
#-----K--------------
elsif($ARGV[1] eq "Kobas_file_to_annolist"){  #Kobas anno file to table [Kobas_file_to_annolist]
	&Kobas_file_to_annolist(@ARGV_info);	
}
#-----L--------------	
elsif($ARGV[1] eq "Listallelement"){ #list all element (^t And " " set apart element) [Listallelement 1]
	&Listallelement(@ARGV_info);	
}
elsif($ARGV[1] eq "ListsubtypeMember"){  #List subtype Member (Type subtype) [ListsubtypeMember 0 1]
	&ListsubtypeMember(@ARGV_info);	
}
elsif($ARGV[1] eq "ListsubtypeNumber"){  #List subtype number by list (Type subtype) [ListsubtypeNumber 0 2 + -]
	&ListsubtypeNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "ListandFindMaxSubtypeNumber"){  #List subtype number by list and find morest subtype (Type subtype) [ListandFindMaxSubtypeNumber 0 2 + -]
	&ListandFindMaxSubtypeNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "Listsubtypevalue"){  #List subtype value by list (Type subtype value) [Listsubtypevalue 0 1 5 SAT NIV GLA BAR GLU MER]
	&Listsubtypevalue(@ARGV_info);	
}
elsif($ARGV[1] eq "ListLineElementNo"){ #ListLineElementNo [ListLineElementNo B2N B4N BDL ETL HVL JAP OTL SBL SIL TAL ZJL ZML]
	&ListLineElementNo(@ARGV_info);	
}
elsif($ARGV[1] eq "ListLinetypevalue"){  #List Line type value (name type1 value1 type2 value2) [ListLinetypevalue 0M 1M 2M 3M 4M 5M]
	&ListLinetypevalue(@ARGV_info);	
}
elsif($ARGV[1] eq "LocalCombinbyDirectandDistancd"){  #Local Combin by Direct and Distancd (name chrom st ed tst ted direct) [LocalCombinbyDirectandDistancd]
	&LocalCombinbyDirectandDistancd(@ARGV_info);	
}
elsif($ARGV[1] eq "LTR_retriever_get_seq_loc"){  #LTR_retriever_get_seq_loc, full,LTR5,LTR3,IN [LTR_retriever_get_seq_loc]
	&LTR_retriever_get_seq_loc(@ARGV_info);	
}
elsif($ARGV[1] eq "LTR_RT_rename_2021"){  #LTR_RT_rename_2021[LTR_RT_rename_2021]
	&LTR_RT_rename_2021(@ARGV_info);	
}
#-----M--------------
elsif($ARGV[1] eq "meg_2_fas"){  #meg_2_fas  [meg_2_fas]
	&meg_2_fas(@ARGV_info);	
}
elsif($ARGV[1] eq "MakerAnnotationRename"){  #Maker Annotation Rename  [MakerAnnotationRename WR]
	&MakerAnnotationRename(@ARGV_info);	
}
elsif($ARGV[1] eq "MakerAnnotationRenameConn"){  #Maker Annotation Rename (mRNA connected) [MakerAnnotationRenameConn WR]
	&MakerAnnotationRenameConn(@ARGV_info);	
}
elsif($ARGV[1] eq "MatrixRowLineChange"){  #Matrix Row Line Change [MatrixRowLineChange]
	&MatrixRowLineChange(@ARGV_info);	
}
elsif($ARGV[1] eq "MatrixToFrame"){  #Matrix to Frame (in R) [MatrixToFrame head]
	&MatrixToFrame(@ARGV_info);	
}
elsif($ARGV[1] eq "MatrixToList"){  #Matrix to Frame (with head) [MatrixToList]
	&MatrixToList(@ARGV_info);	
}
elsif($ARGV[1] eq "MergedGtfGetTranscriptsAndGene"){  #Merged Gtf Get Transcripts and Gene [MergedGtfGetTranscriptsAndGene]
	&MergedGtfGetTranscriptsAndGene(@ARGV_info);	
}
elsif($ARGV[1] eq "MethylationTypeIdentityINGenome"){  #Check the Methylation Type Identity In input Genome [MethylationTypeIdentityINGenome]
	&MethylationTypeIdentityINGenome(@ARGV_info);	
}
elsif($ARGV[1] eq "Methy_gene_Local_count"){  #count gene up down in meth % [Methy_gene_Local_count 50]
	&Methy_gene_Local_count(@ARGV_info);	
}
elsif($ARGV[1] eq "MisaResultSeperate"){  #Misa Result Seperate to lines [MisaResultSeperate]
	&MisaResultSeperate(@ARGV_info);	
}
elsif($ARGV[1] eq "MisaResult_samedirect_Seperate"){  #Misa Result Seperate to lines [MisaResult_samedirect_Seperate]
	&MisaResult_samedirect_Seperate(@ARGV_info);	
}
elsif($ARGV[1] eq "MummerResultCombin"){  #Mummer Result Combin [MummerResultCombin]
	&MummerResultCombin(@ARGV_info);	
}	
elsif($ARGV[1] eq "MutationCountForCode123"){  #Mutation Count For Code 1 2 3 (input must be aligned fas) [MutationCountForCode123 seedname]
	&MutationCountForCode123(@ARGV_info);	
}	
#-----P--------------
elsif($ARGV[1] eq "PrankforCombinedFile"){ #Prank for Combined File, e.g. 6_align.fas [PrankforCombinedFile]
	&PrankforCombinedFile(@ARGV_info);
}
elsif($ARGV[1] eq "PranknexandEventsMutateCount"){ #Prank out .best.nex Count indel [PranknexandEventsMutateCount]
	&PranknexandEventsMutateCount(@ARGV_info);
}
elsif($ARGV[1] eq "PrankEventsCount"){ #Prank out .events Count indel [PrankEventsCount]
	&PrankEventsCount(@ARGV_info);
}
elsif($ARGV[1] eq "PPI_Link_Clean"){ #Link Clean for STRING [PPI_Link_Clean]
	&PPI_Link_Clean(@ARGV_info);	
}
elsif($ARGV[1] eq "PvalueFoldchangeCheck"){ #Pvalue Fold changeCheck (ID p f p f ...) [PvalueFoldchangeCheck]
	&PvalueFoldchangeCheck(@ARGV_info);
}
#-----R--------------
elsif($ARGV[1] eq "ReadCLUSTALfilebyList"){ #ReadCLUSTALfilebyList  [ReadCLUSTALfilebyList]
	&ReadCLUSTALfilebyList(@ARGV_info);
}
elsif($ARGV[1] eq "ReadinfilebyList"){ #ReadinfilebyList and add type by file name  [ReadinfilebyList]
	&ReadinfilebyList(@ARGV_info);
}
elsif($ARGV[1] eq "ReversedOrderForFasta"){ #Reversed Order For Fasta file  [ReversedOrderForFasta 1]
	&ReversedOrderForFasta(@ARGV_info);
}
elsif($ARGV[1] eq "RreturnLocTorows"){ #return Loc info To rows  [RreturnLocTorows 1]
	&RreturnLocTorows(@ARGV_info);
}	
elsif ($ARGV[1] eq "RegulateFastaFile"){ #Regulate Fasta File [RegulateFastaFile 70]
	&RegulateFastaFile(@ARGV_info);	
}
elsif($ARGV[1] eq "Regulateblasttoloc"){ #regulate blast result To Loc info (use target)(If have parametre "+",  amend the result£©  [Regulateblasttoloc 1 +]
	&regulateblasttoloc(@ARGV_info);	
}
elsif($ARGV[1] eq "RegulateblastQuerrytoloc"){ #regulate blast result To Loc info (use querry)(If have parametre "+",  amend the result£©  [RegulateblastQuerrytoloc x +]
	&RegulateblastQuerrytoloc(@ARGV_info);	
}	
elsif($ARGV[1] eq "Regulatestsfile"){ #regulate sts file (ID L-primer R-primer) Name are same in two side [Regulatestsfile 1]
	&regulatestsfile(@ARGV_info);
}
elsif($ARGV[1] eq "RemoveAllSameLineExceptRows"){ #Remove All Same Line Except para Rows [RemoveAllSameLineExceptRows 0 2 3 4]
	&RemoveAllSameLineExceptRows(@ARGV_info);	
}
elsif($ARGV[1] eq "Remove_0_row_percentage"){ #Remove row had many 0, more than cutoff (defult 0.99) [Remove_0_row_percentage 0.99]
	&Remove_0_row_percentage(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveLineWtihAllSameRow"){ #Remove line with All Same row Except first [RemoveLineWtihAllSameRow]
	&RemoveLineWtihAllSameRow(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveLinesAll0"){ #Remove line with All 0 row Except first [RemoveLinesAll0]
	&RemoveLinesAll0(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveLinesBytotalValue"){ #Remove Lines By total Value [RemoveLinesBytotalValue 10]
	&RemoveLinesBytotalValue(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveLinefromNumber"){ #Remove line from number [RemoveLinefromNumber 1 3 4]
	&RemoveLinefromNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "RemvoeChar"){ #remvoe Char [RemvoeChar -]
	&RemvoeChar(@ARGV_info);
}
elsif($ARGV[1] eq "RemoveSameNameFas"){ #Remove Same Name Fas [RemoveSameNameFas]
	&RemoveSameNameFas(@ARGV_info);
}
elsif($ARGV[1] eq "RemoveInfoAfterSymble"){ #Remove Info After Symble [RemoveInfoAfterSymble _]
	&RemoveInfoAfterSymble(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveNestingCodeinGene"){ #Remove Nesting Code in same Gene [RemoveNestingCodeinGene 500]
	&RemoveNestingCodeinGene(@ARGV_info);	
}	
elsif($ARGV[1] eq "Remove0behandnumber"){ #remove 0 behand number (row) [Remove0behandnumber 2]
	&remove0behandnumber(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveExtremeValuebyType"){ #Remove Extreme Value by Type(Name Type Value) [RemoveExtremeValuebyType 2]
	&RemoveExtremeValuebyType(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveElementSameRow"){ #Remove Element Same Row [RemoveElementSameRow 1 2 3 4 5 6]
	&RemoveElementSameRow(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveAllSameRowbyType"){ #Remove Element Same Row (Type other_same_row)[RemoveAllSameRowbyType 4 9 10 11 12 13 14 15 16 17]
	&RemoveAllSameRowbyType(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveSameCodeforAlignedFas"){ #Remove Same Code for Aligned Fas (1 for single remove, 3 for 3 code remove)[RemoveSameCodeforAlignedFas]
	&RemoveSameCodeforAlignedFas(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveSequenceBy_NNN_MoreThanCutoff"){ #Remove Sequence By NNN More Than Cutoff 20[RemoveSequenceBy_NNN_MoreThanCutoff 20]
	&RemoveSequenceBy_NNN_MoreThanCutoff(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveSequencWithSpecialCharacters"){ #Remove Sequenc With Special Characters other than ATGC[RemoveSequencWithSpecialCharacters]
	&RemoveSequencWithSpecialCharacters(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatnoCombinCount"){ #repeatmasker cat Count (never combin, family cut by input) [RepeatMaskerCatnoCombinCount _]
	&RepeatMaskerCatnoCombinCount(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatCount"){ #repeatmasker cat Count (use this 1) [RepeatMaskerCatCount]
	&RepeatMaskerCatCount(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatCombinCount"){ #repeatmasker cat Count (use this 1, combin to @@@) [RepeatMaskerCatCombinCount]
	&RepeatMaskerCatCombinCount(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatOUTfamilyCombin"){ #repeatmasker cat out same family combin  (use this 2)  [RepeatMaskerCatOUTfamilyCombin 0 0]
	&RepeatMaskerCatOUTfamilyCombin(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatOUTTypeSimpleCombin"){ #repeatmasker cat sample combin  (use this only one) (type row) [RepeatMaskerCatOUTTypeSimpleCombin out]
	&RepeatMaskerCatOUTTypeSimpleCombin(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatOUTTypeLoopCombin"){ #repeatmasker cat sample combin  (**use this only one) (type row) [RepeatMaskerCatOUTTypeLoopCombin out]
	&RepeatMaskerCatOUTTypeLoopCombin(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerOUTTypeLoopCombin"){ #repeatmasker out sample combin  (**use this only one) (type row) [RepeatMaskerOUTTypeLoopCombin]
	&RepeatMaskerOUTTypeLoopCombin(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatOUTtypetotalCombin"){ #repeatmasker cat out same sub type combin  (use this combin)  [RepeatMaskerCatOUTtypetotalCombin 50]
	&RepeatMaskerCatOUTtypetotalCombin(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerCatLTRCombin"){ #RepeatMasker Cat count result .out.LTR.detail result Combin (LTR use this 3) [RepeatMaskerCatLTRCombin 10000 100 S]
	&RepeatMaskerCatLTRCombin(@ARGV_info);	
}	
elsif ($ARGV[1] eq "RepeatMaskerCatDetailCount"){ #repeatmasker cat Count [RepeatMaskerCatDetailCount _]
	&RepeatMaskerCatDetailCount(@ARGV_info);	
}	
elsif ($ARGV[1] eq "RepeatMaskerCat_2_out"){ #repeatmasker cat to out [RepeatMaskerCat_2_out]
	&RepeatMaskerCat_2_out(@ARGV_info);	
}	
elsif ($ARGV[1] eq "RepeatMaskerOutPercent"){ #RepeatMasker Out Percent [RepeatMaskerOutPercent]
	&RepeatMaskerOutPercent(@ARGV_info);	
}
elsif ($ARGV[1] eq "RepeatMaskerOutReadsType"){ #RepeatMasker Out Reads Type [RepeatMaskerOutReadsType]
	&RepeatMaskerOutReadsType(@ARGV_info);	
}
elsif ($ARGV[1] eq "REannotatedCount"){ #Reannotated result  Count [REannotatedCount 1]
	&REannotatedCount(@ARGV_info);	
}	
elsif ($ARGV[1] eq "Replace_Low_High_to_cutoff"){ #Replace_Low_High_to_digit (gt/lt cutoff)[Replace_Low_High_to_cutoff lt 0.0001]
	&Replace_Low_High_to_cutoff(@ARGV_info);	
}
elsif ($ARGV[1] eq "ReplaceSymble1by2"){ #replace symble1 by symble2 (Change in code) [ReplaceSymble1by2 2]
	&ReplaceSymble1by2(@ARGV_info);	
}
elsif ($ARGV[1] eq "ReplaceFastAOtherLetter"){ #replace other letter (except atgc)fasta formate seq [ReplaceFastAOtherLetter 1 _]
	&ReplaceFastAOtherLetter(@ARGV_info);	
}
elsif($ARGV[1] eq "Removefilehead"){ #Remove file head rows  [Removefilehead 3]
	&Removefilehead(@ARGV_info);	
}
elsif($ARGV[1] eq "RemoveTotalIncludeData"){ #Remove Total Include Data (chr start length) [RemoveTotalIncludeData 2 3 4]
	&RemoveTotalIncludeData(@ARGV_info);	
}	
elsif($ARGV[1] eq "Renameclass"){ #rename class (name start_from length) [Renameclass TL 1 3]
	&Renameclass(@ARGV_info);	
}
elsif($ARGV[1] eq "RenamebyTypelist"){ #rename sequence by rule (type row, digital) [RenamebyTypelist 0 3]
	&RenamebyTypelist(@ARGV_info);	
}	
elsif($ARGV[1] eq "Renamesequence"){ #rename sequence by rule (first_number symble Add_to_4_digit) [Renamesequence 1 GRA 5]
	&renamesequence(@ARGV_info);	
}	
elsif($ARGV[1] eq "RenameFamilysequence"){ #rename sequence by rule hold family name [RenameFamilysequence 1 ruf 3]
	&RenameFamilysequence(@ARGV_info);	
}	
elsif($ARGV[1] eq "ReverseFilebyRow"){ #Reverse File by Row [ReverseFilebyRow]
	&ReverseFilebyRow(@ARGV_info);	
}
elsif($ARGV[1] eq "RNAhybridOutputCheck"){ #RNAhybrid Output Check [RNAhybridOutputCheck]
	&RNAhybridOutputCheck(@ARGV_info);	
}
#-----S--------------
elsif($ARGV[1] eq "SamRemoveSinglehit"){ #Sam file Remove Single hit in one place [SamRemoveSinglehit]
	&SamRemoveSinglehit(@ARGV_info);	
}	
elsif($ARGV[1] eq "SearchAlignPercentfromPhyFile"){ #1 Search Align Percent from Phy File from window (now only 10) [SearchAlignPercentfromPhyFile]
	&SearchAlignPercentfromPhyFile(@ARGV_info);	
}	
elsif($ARGV[1] eq "SearchAlignPefectRegion"){ #2 Search Align Pefect Region (6 means >60%) [SearchAlignPefectRegion 6]
	&SearchAlignPefectRegion(@ARGV_info);	
}	
elsif($ARGV[1] eq "SearchAlignRegionGetFas"){ #3 Search Align RegionGet Fas [SearchAlignPefectRegion 6]
	&SearchAlignRegionGetFas(@ARGV_info);	
}
elsif($ARGV[1] eq "SlideWindowCount"){ #Count date by Slide Window (window size, step size)  [SlideWindowCount 10 1]
	&SlideWindowCount(@ARGV_info);	
}	
elsif($ARGV[1] eq "SlideWindowCount_sorttype"){ #Count date by Slide Window for sorted type (0-type 1-loc 2-value) (window size, step size)  [SlideWindowCount_sorttype 10 1]
	&SlideWindowCount_sorttype(@ARGV_info);	
}
elsif($ARGV[1] eq "SlideWindow_continue_step"){ #Count Slide Window by continue sorted steps (0-chrom 1-st 2-ed 3-value) (window = 10 step)  [SlideWindow_continue_step 10]
	&SlideWindow_continue_step(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateFastAbySeq"){ #seperate fas seq to single file by sequence  (fa for .fa) [SeperateFastAbySeq aatgc]
	&SeperateFastAbySeq(@ARGV_info);	
}  
elsif($ARGV[1] eq "SeperateGroupFastAtoFile"){ #seperate each seq to a single file by name *** caution [SeperateGroupFastAtoFile]
	&SeperateGroupFastAtoFile(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperatFastAbyNumber"){ #seperate seq to file by number [SeperatFastAbyNumber 100]
	&SeperatFastAbyNumber(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateFastAbyLength"){ #seperate seq to file by fas length [SeperateFastAbyLength 500000000]
	&SeperateFastAbyLength(@ARGV_info);	
}
elsif($ARGV[1] eq "Seperate_cut_FastAbyLength"){ #cut seq to file by fas length [Seperate_cut_FastAbyLength 5000]
	&Seperate_cut_FastAbyLength(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateEachLineValue"){ #Seperate Each Line as a file use title [SeperateEachLineValue]
	&SeperateEachLineValue(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateFilebyRowType"){ #seperate each type to a single file by row [SeperateFilebyRowType 0]
	&SeperateFilebyRowType(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateFilebyLine"){ #Seperate File by each 2000000 Line [SeperateFilebyLine 2000000]
	&SeperateFilebyLine(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateFilebySortedRowType"){ #seperate each type to a single file by row #sorted [SeperateFilebySortedRowType 0]
	&SeperateFilebySortedRowType(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateLinebyType"){ #Seperate each name a Line by Type (Type row) [SeperateLinebyType 0]
	&SeperateLinebyType(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateGOannoResult"){ #seperate to each line a gene [SeperateGOannoResult]
	&SeperateGOannoResult(@ARGV_info);
}
elsif($ARGV[1] eq "SeperateGroupPair"){ #Seperate seperate group pair [SeperateGroupPair 5]
	&SeperateGroupPair(@ARGV_info);	
}
elsif($ARGV[1] eq "Seperategff3locToregionloc"){ #seperate each group seq to a single file by name (name mast be )[Seperategff3locToregionloc]
	&Seperategff3locToregionloc(@ARGV_info);	
}
elsif($ARGV[1] eq "Seperate_Scaffold_2_contig"){ #Seperate_Scaffold_2_contig, cut from NNN [Seperate_Scaffold_2_contig]
	&Seperate_Scaffold_2_contig(@ARGV_info);	
}
elsif($ARGV[1] eq "SeperateTreesReadfromfilelist"){ #Seperate Trees Read and clean from filelist (.dnd)[SeperateTreesReadfromfilelist]
	&SeperateTreesReadfromfilelist(@ARGV_info);	
}
elsif($ARGV[1] eq "Seperate_or_Type"){ #Seperate or Type such as a|b|c [Seperate_or_Type 5]
	&Seperate_or_Type(@ARGV_info);	
}
elsif($ARGV[1] eq "SequenceRemoveGapLocCorrespondence"){ #Sequence Remove Gap Local Correspondence 1,4,5,6,8,9 (loss number is gap) [SequenceRemoveGapLocCorrespondence]
	&SequenceRemoveGapLocCorrespondence(@ARGV_info);	
}
elsif($ARGV[1] eq "SelectBuyRowValue"){ #select buy all row value (lt = gt eq ne) [SelectBuyRowValue 8 gt 1]
	&SelectBuyRowValue(@ARGV_info);	
}
elsif($ARGV[1] eq "SelectBuyAtleastOneRowValue"){ #select buy at least one row value (lt = gt eq ne) [SelectBuyAtleastOneRowValue 1 2 3 4 gt 1]
	&SelectBuyAtleastOneRowValue(@ARGV_info);	
}
elsif($ARGV[1] eq "SelectbyRowString"){ #select by row String [SelectbyRowString 25 ne NA]
	&SelectbyRowString(@ARGV_info);
}	
elsif($ARGV[1] eq "SelectbyRowStringLength"){ #select by row String length  (lt eq gt) [Rowlengthcutoff 6 gt 95]
	&rowlengthcutoff(@ARGV_info);
}	
elsif($ARGV[1] eq "SelectMutiRowbyValue"){ #select buy row value (lt = gt eq ne) [SelectMutiRowbyValue 1 gt 0 2 gt 0 2 gt 0 2 gt 0]
	&SelectMutiRowbyValue(@ARGV_info);	
}
elsif ($ARGV[1] eq "Selectsequencebylength"){ #select sequence by length  [Selectsequencebylength 500]
	&selectsequencebylength(@ARGV_info);
}
elsif ($ARGV[1] eq "Selectfastqbylength"){ #select fastq seq by length  [Selectfastqbylength lt 80]
	&Selectfastqbylength(@ARGV_info);
}
elsif ($ARGV[1] eq "SelectLargestValueLine"){ #Select  Line with Largest Value (type value) [SelectLargestValueLine 0 1]
	&SelectLargestValueLine(@ARGV_info);
}
elsif ($ARGV[1] eq "SelectbyLargestSubtypeValue"){ #Select  Line with Largest subtype Value £¨cut off must > cutoff such as 20%, if not mark to unknown type£© (input must be: type subtype value) [SelectbyLargestSubtypeValue 40]
	&SelectbyLargestSubtypeValue(@ARGV_info);
}
elsif ($ARGV[1] eq "SignlocalLargeThanMedian"){ #Sign 1 to local Large Than ? times of Median, other sign 0 (median_row gt/lt  times) [SignlocalLargeThanMedian 6 gt 1]
	&SignlocalLargeThanMedian(@ARGV_info);
}		
elsif($ARGV[1] eq "SoapDirectCount"){ #soap result .pe, count direct  mate<- -> and base-><- [SoapDirectCount]
	&SoapDirectCount(@ARGV_info);	
}
elsif($ARGV[1] eq "SoapDepthCount"){ #(soap result .pe + .se, remove seq info, sort by 6,7) then count windows depth [SoapDepthCount]
	&SoapDepthCount(@ARGV_info);	
}
elsif($ARGV[1] eq "SoapDepthFastCount"){ #(soap result .pe + .se, remove seq info, sort by 6,7) then count windows depth [SoapDepthFastCount]
	&SoapDepthFastCount(@ARGV_info);	
}
elsif($ARGV[1] eq "SortbyrownumberNR"){ #sort String by row number (+:A-Z -:Z-A ) //fast,  but can't handle repeat [SortbyrownumberNR 1 +]
	&Sortbyrownumber(@ARGV_info);	
}
elsif($ARGV[1] eq "SortbyrownumberRepeat"){ #sort String by row number (+:A-Z -:Z-A ) [SortbyrownumberRepeat 1 +]
	&SortbyrownumberRepeat(@ARGV_info);	
}
elsif($ARGV[1] eq "sort2IntinEachRow"){ #sort 2 Int in Each Row (1 < 2) [sort2IntinEachRow 1 2]
	&sort2IntinEachRow(@ARGV_info);	
}
elsif($ARGV[1] eq "sortbyIntrow"){ #sort by number row (row direct)+:from small to large  -:from large to small [sortbyIntrow 8 +]
	&sortbyIntrow(@ARGV_info);	
}
elsif($ARGV[1] eq "sortbyIntandStrrow"){ #sort by first int and second string row (row direct)+:from small to large  -:from large to small [sortbyIntandStrrow 8 3 +]
	&sortbyIntandStrrow(@ARGV_info);	
}
elsif($ARGV[1] eq "sortbyIntOrStrrow"){ #sort by more than two row row  [sortbyIntOrStrrow str 12 str 3 int 4]
	&sortbyIntOrStrrow(@ARGV_info);	
}        
elsif($ARGV[1] eq "sortbyStrandIntrow"){ #sort by first String and second int row (row direct)+:from small to large  -:from large to small [sortbyStrandIntrow 8 2 +]
	&sortbyStrandIntrow(@ARGV_info);	
}
elsif($ARGV[1] eq "sortby2Introw"){ #sort by 2 int row (row direct)+:from small to large  -:from large to small [sortby2Introw 1 4 +]
	&sortby2Introw(@ARGV_info);	
}
elsif($ARGV[1] eq "sortbyStringrow"){ #sort by String row (row direct)+:from small to large  -:from large to small [sortbyStringrow 1 -]
	&sortbyStringrow(@ARGV_info);
}
elsif($ARGV[1] eq "sortbyChr"){ #sort by chr row (loc sorted) [sortbyChr 0]
	&sortbyChr(@ARGV_info);
}
elsif($ARGV[1] eq "sortRow"){ #sort Row [sortRow]
	&sortRow(@ARGV_info);	
}
elsif($ARGV[1] eq "sortSequencebyName"){ #sort Sequence by Name [sortSequencebyName]
	&sortSequencebyName(@ARGV_info);	
}
elsif($ARGV[1] eq "sortSequencebyLength"){ #sort Sequence by length [sortSequencebyLength]
	&sortSequencebyLength(@ARGV_info);	
}
elsif($ARGV[1] eq "SSRCountConsensusNumber"){ #SSR count Consensus sequence number, eg 2-4 3-10  TA->AT, TCT->AGA [SSRCountConsensusNumber 4]
	&SSRCountConsensusNumber(@ARGV_info);
}
elsif($ARGV[1] eq "SSRFindConsensusSequence"){ #SSR find sort Consensus sequence, eg TA->AT, TCT->AGA [SSRFindConsensusSequence]
	&SSRFindConsensusSequence(@ARGV_info);
}
elsif($ARGV[1] eq "StasticEachRowDivideMean"){ #Stastic Each Row Divide by Mean [StasticEachRowDivideMean]
	&StasticEachRowDivideMean(@ARGV_info);
}
elsif($ARGV[1] eq "StasticCountAlternativeSplicingTEloc"){ #Stastic Count Alternative Splicing TE loc [StasticCountAlternativeSplicingTEloc]
	&StasticCountAlternativeSplicingTEloc(@ARGV_info);
}
elsif($ARGV[1] eq "StasticCountAverageBy0or1"){ #find Stop coden loc and count [StasticCountAverageBy0or1]
	&StasticCountAverageBy0or1(@ARGV_info);
}
elsif($ARGV[1] eq "StasticGetRandlines"){ #Stastic Get Rand lines by number [StasticGetRandlines 100]
	&StasticGetRandlines(@ARGV_info);
}	
elsif($ARGV[1] eq "StopCodenCount"){ #find Stop coden loc and count [StopCodenCount]
	&StopCodenCount(@ARGV_info);
}		
#-----T--------------
elsif($ARGV[1] eq "TableCombinTo01Matrix"){ #Table combin To 01 Matrix by cut off  [TableCombinTo01Matrix 0.5]
	&TableCombinTo01Matrix(@ARGV_info);
}	
elsif($ARGV[1] eq "TableChangeTo01ByCutoff"){ #Table Change To 0 1 By Cutoff  [TableChangeTo01ByCutoff 60]
	&TableChangeTo01ByCutoff(@ARGV_info);
}	
elsif($ARGV[1] eq "TableChangeTo01ByMutiIf"){ #Table Change To 0 1 By Muti If  [TableChangeTo01ByMutiIf]
	&TableChangeTo01ByMutiIf(@ARGV_info);
}
elsif($ARGV[1] eq "TASSEL_distance_2_phylip"){ #change TASSEL_distance_2_phylip  [TASSEL_distance_2_phylip]
	&TASSEL_distance_2_phylip(@ARGV_info);
}
elsif($ARGV[1] eq "TASSEL_glm_2_Manhattan_input"){ #TASSEL_glm_2_Manhattan_input  [TASSEL_glm_2_Manhattan_input]
	&TASSEL_glm_2_Manhattan_input(@ARGV_info);
}
elsif($ARGV[1] eq "TASSEL_mlm_2_Manhattan_input"){ #TASSEL_mlm_2_Manhattan_input  [TASSEL_mlm_2_Manhattan_input]
	&TASSEL_mlm_2_Manhattan_input(@ARGV_info);
}
elsif($ARGV[1] eq "ThreeSamplesUpDownCompare"){ #Three Samples Up Down Compare  [ThreeSamplesUpDownCompare]
	&ThreeSamplesUpDownCompare(@ARGV_info);
}
elsif($ARGV[1] eq "TrimFastaHeadinfo"){ #Trim Fasta Headinfo, remain some row  [TrimFastaHeadinfo 1]
	&TrimFastaHeadinfo(@ARGV_info);
}
elsif($ARGV[1] eq "Typelowvaluemask"){ #add symble To Type low value (Input£ºType value)  [Typelowvaluemask 1]
	&Typelowvaluemask(@ARGV_info);
}
elsif($ARGV[1] eq "TransposonsNestCombin"){ #Transposons Nest Combin, DNA transposons priority level higher than LTR transposons  [TransposonsNestCombin]
	&TransposonsNestCombin(@ARGV_info);
}
elsif($ARGV[1] eq "trinotate_annotation_venn"){ #trinotate annotation venn [trinotate_annotation_venn]
	&trinotate_annotation_venn(@ARGV_info);
}
elsif($ARGV[1] eq "truncFloatToInt"){ #truncFloatToInt [truncFloatToInt]
	&truncFloatToInt(@ARGV_info);
}
#-----V--------------
elsif($ARGV[1] eq "vcf_filter"){ #vcf filter (depth cut, percentage cut) [vcf_filter 8 1]
	&vcf_filter(@ARGV_info);
}
elsif($ARGV[1] eq "vcf_2_MatrixEQTL_SNP"){ #vcf_2_MatrixEQTL_SNP [vcf_2_MatrixEQTL_SNP]
	&vcf_2_MatrixEQTL_SNP(@ARGV_info);
}
elsif($ARGV[1] eq "venn_state_by_filelist"){ #make state table by filelist [venn_state_by_filelist]
	&venn_state_by_filelist(@ARGV_info);
}
#--------------------
elsif($ARGV[1] eq "test_test"){ #test [test]
	&test(@ARGV_info);
}
else{
	print "Nothing had been done!!?????\n";
}

sub test{
	#for(my $locn=1;$locn<=20;$locn++){
		#my $y=0.5-atan2($x,1)/3.14;
		#my $aa = 1-(-0.0002*($locn**3) + 0.007*($locn**2) - 0.0141*$locn + 0.0051);
		#my $y=log($x)/(log(10)*1.5);
		#my $y=atan2($x,1);
		#my $y = -0.0002*($x**3) + 0.007*($x**2) - 0.0141*$x + 0.0051;
		#my $y = 0.0005*($x**3) - 0.0652*($x**2) + 0.8561*$x + 97.812;
		#print "$locn $aa\n";
	#}
	my @a=(1,2,3,4);
	@b=splice(@a,2,1);
	
	print "@a\n@b\n";
}



#------------------


sub LTR_retriever_get_seq_loc{
	
	##LTR_loc        Category        Motif   TSD     5_TSD 3_TSD       Internal        Identity      Strand  SuperFamily  TE_type     Insertion_Time
	#000258F|arrow_u:523762..536026  pass    motif:TGCA      TSD:AAAGC       523757..523761  536027..536031  IN:526819..532967       0.9768  ?       unknown NA      907709
	#000258F|arrow_u:731472..738885  pass    motif:TGCA      TSD:CATAG       731467..731471  738886..738890  IN:732813..737544       0.9933  +       Gypsy   LTR     259293
	#000258F|arrow_u:751808..756531  pass    motif:TGCA      TSD:CTAAT       751803..751807  756532..756536  IN:752129..756210       0.9875  -       Copia   LTR     483297

	my $name=$_[1];
	my $n=0;
	$lineIN=<IN>;
	
	open(OUTA,'>',$ARGV[0].".all");
	open(OUTL5,'>',$ARGV[0].".LTR5");
	open(OUTL3,'>',$ARGV[0].".LTR3");
	open(OUTIN,'>',$ARGV[0].".IN");
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my @all=&Getheadinfo($tp[0],":",".");
		my @in=&Getheadinfo($tp[6],":",".");
		if($tp[8] eq "?"){
			$tp[8]="+";
		}
		print OUTA "$name\_$n\_$tp[9]\t$all[0]\t$all[1]\t$all[2]\t$tp[8]\n";
		print OUTIN "$name\_$n\_$tp[9]\t$all[0]\t$in[1]\t$in[2]\t$tp[8]\n";
		if($tp[8] eq "+"){
			$ed=$in[1]-1;
			print OUTL5 "$name\_$n\_$tp[9]\t$all[0]\t$all[1]\t$ed\t$tp[8]\n";
			$st=$in[2]+1;
			print OUTL3 "$name\_$n\_$tp[9]\t$all[0]\t$st\t$all[2]\t$tp[8]\n";
		}else{
			$ed=$in[1]-1;
			print OUTL3 "$name\_$n\_$tp[9]\t$all[0]\t$all[1]\t$ed\t$tp[8]\n";
			$st=$in[2]+1;
			print OUTL5 "$name\_$n\_$tp[9]\t$all[0]\t$st\t$all[2]\t$tp[8]\n";			
		}
		$n++;
	}
	print "$n done\n";
	close(OUTA);
	close(OUTL5);
	close(OUTL3);
	close(OUTIN);
}





sub LTR_RT_rename_2021{
	
	#>2G_112_Copia|copia=Maximus=Mtr12.5     ORF=3   StopCodon=3     6,19,29,
	
	open(OUTR,'>',$ARGV[0].".rename");
	%remame = ('Maximus','Max','Tat', 'Tat', 'Tekay', 'Tky', 'Athlia', 'Ath', 'Reina', 'Rei', 'CRM', 'CRM', 'Angela', 'Ang', 'GMR', 'GMR', 'TAR', 'TAR', 'Bianca', 'Bia', 'Ale', 'Ale', 'Ivana', 'Iva');

	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp = &Getheadinfo($lineIN,"_","=","|");
			#my $o=">$tp[1]_$remame{$tp[3]}";
			my $o=">$tp[3]_$tp[1]_$tp[4]";
			print OUT "$o\n";
			print OUTR "$lineIN\t$o\n";
		}else{
			print OUT "$lineIN\n";
		}		
	}
	print "done\n";
	close(OUTR);
}


sub ConnectExcelEnter{
	
	#"Fatty Acid Elongation & Wax Biosynthesis
	#Cutin Synthesis & Transport 1
	#Suberin Synthesis & Transport 1"

	my $add=0;
	my $o="";
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/\"/){
			if($add==0){
				$add=1;
				$o=$lineIN;
			}else{
				$add=0;
				$o=$o." / ".$lineIN;
				print OUT "$o\n";
			}
		}else{
			if($add==0){
				print OUT "$lineIN\n";
			}else{
				$o=$o.$lineIN;
			}
		}
		
	}
	print "done\n";
}



sub bismark_2_methylKit{
	#in:Chromosome	Site	Direct	meth	un_meth	type	context	
	#Chr1_RaGOO      9       +       1       0       CHH     CTA
	
	print OUT "chrBase	chr	base	strand	coverage	freqC	freqT\n";
	
	while($lineIN=<IN>){
		if($lineIN=~/rev_G_count/){
			print "..";
			$lineIN=<IN>;
		}
		chomp($lineIN);
		my @tp = split(/\t/,$lineIN);
		my $chrBase="$tp[0].$tp[1]";
		my $depth=$tp[3]+$tp[4];
		my $freqC=sprintf("%.2f",(100*$tp[3]/$depth));
		my $freqT=sprintf("%.2f",(100-$freqC));	
		my $strand="F";
		if($tp[2] eq "-"){
			$strand="R";
		}
		print OUT "$chrBase\t$tp[0]\t$tp[1]\t$strand\t$depth\t$freqC\t$freqT\n";
	}
	print "done\n";
}



sub Bsmap_2_methylKit{
	#in
	#chr     pos     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower        CI_upper
	#Chr10_RaGOO     2688    +       CG      0.020   51.00   1       51      35      35      0.003   0.103
	#Chr10_RaGOO     2689    -       CG      0.639   36.00   23      36      50      50      0.476   0.775	
	#chrBase	chr	base	strand	coverage	freqC	freqT
	#Chr10_RaGOO.2688	Chr10_RaGOO	2688	F	51.00	2.0	98.0
	
	print OUT "chrBase	chr	base	strand	coverage	freqC	freqT\n";
	
	while($lineIN=<IN>){
		if($lineIN=~/rev_G_count/){
			print "..";
			$lineIN=<IN>;
		}
		chomp($lineIN);
		my @tp = split(/\t/,$lineIN);
		my $chrBase="$tp[0].$tp[1]";
		my $freqC=sprintf("%.2f",($tp[4]*100));
		my $freqT=sprintf("%.2f",(100-$freqC));	
		my $strand="F";
		if($tp[2] eq "-"){
			$strand="R";
		}
		print OUT "$chrBase\t$tp[0]\t$tp[1]\t$strand\t$tp[5]\t$freqC\t$freqT\n";
	}
	print "done\n";
}



sub GetMax{
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp = split(/\t/,$lineIN);
		my $o=$lineIN[$_[1]]>$lineIN[$_[2]]?$lineIN[$_[1]]:$lineIN[$_[2]];
		print OUT "$lineIN\t$o\n";
	}
	print "done\n";
}

sub GetMin{
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp = split(/\t/,$lineIN);
		my $o=$lineIN[$_[1]]<$lineIN[$_[2]]?$lineIN[$_[1]]:$lineIN[$_[2]];
		print OUT "$lineIN\t$o\n";
	}
	print "done\n";
}


sub meg_2_fas{

	$lineIN=<IN>;
	$lineIN=<IN>;
	$lineIN=<IN>;
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		if(length($lineIN)>0){
			if ($lineIN=~/#/){
				$lineIN =~ s/#/>/;
			}
			print OUT "$lineIN\n";
		}
	}
	print "done\n";
}



sub SeperateFastAbySeq{
	
	my $name="NA";
	my $cut=$_[1];
	my $seq="";
	my $count=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		if ($lineIN=~/>/){
			if($name ne "NA"){
				my $loc=index($seq,$cut,);
				if($loc>0){
					my $tpseq=substr($seq,0,$loc);
					print OUT "$name"."_"."1\n$tpseq\n";
					$tpseq=substr($seq,$loc,length($seq));
					print OUT "$name"."_"."2\n$tpseq\n";
				}
			}
			$name=$lineIN;
			$seq="";
			$count=0;
		}else{
			$seq=$seq.$lineIN;			
		}

	}
	#last
	my $loc=index($seq,$cut,);
	if($loc>0){
		my $tpseq=substr($seq,0,$loc);
		print OUT "$name"."_"."1\n$tpseq\n";
		$tpseq=substr($seq,$loc,length($seq));
		print OUT "$name"."_"."2\n$tpseq\n";
	}
	print "done\n";
}



sub Seperate_Scaffold_2_contig{
	
	my $name="NA";
	my $cut=10;
	my $seq="";
	my $count=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		if ($lineIN=~/>/){
			if($name ne "NA"){
				my @arr = split('N', $seq);
				my $arr=@arr;
				#print "$name $arr\n";
				for(my $i=0;$i<$arr;$i++){
					if(length($arr[$i])>=$cut){
						print OUT "$name"."_"."$count\n$arr[$i]\n";
						$count++;
					}
				}
			}
			$name=$lineIN;
			$seq="";
			$count=0;
		}else{
			$seq=$seq.$lineIN;			
		}

	}
	#last
	my @arr = split('N', $seq);
	my $arr=@arr;
	#print "$name $arr\n";
	for(my $i=0;$i<$arr;$i++){
		if(length($arr[$i])>=$cut){
			print OUT "$name"."_"."$count\n$arr[$i]\n";
			$count++;
		}
	}
	print "done\n";
}


sub truncFloatToInt{

	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my $l=index($tp[$_[1]],".",);
		#print "$_[1] $tp[$_[1]] $l\n";
		if($l>=0){
			my $a1=substr($tp[$_[1]],0,$l);
			my $a2=substr($tp[$_[1]],($l+1),1);
			#print "$a2\n";
			if($a2>=5){
				$tp[$_[1]]=$a1+1;
			}else{
				$tp[$_[1]]=$a1;
			}
		}
		print OUT "$tp[0]";
		for(my $i=1;$i<$tp;$i++){
			print OUT "\t$tp[$i]";
		}
		print OUT "\n";
	}
	print "done\n";
}

sub EMMAX_ps_2_Manhattan_input{

### for EMMAX.ps
#
#	[out_prefix].ps : Each line consist of
#	1. SNP ID
#	2. Beta (1 is effect allele)
#	3. SE(beta)
#	4. p-value.

### for Manhattan_inpu for  R
#
#	SNP CHR BP         P
# rs1   1  1 0.9148060
# rs2   1  2 0.9370754
#

	print OUT "\tSNP\tCHR\tBP\tP\n";	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp=split(/\t/,$lineIN);
		my @snp=split(/\_/,$tp[0]);
		my $chr=substr($snp[0],1,length($snp[0]));
		print OUT "$tp[0]\t$chr\t$snp[1]\t$tp[3]\n";
	}
	print "done\n";
}


sub TASSEL_glm_2_Manhattan_input{

### for glm_?_1.txt
#
#Trait   Marker  Chr     Pos     marker_F        p       marker_Rsq      add_F   add_p   dom_F   dom_p   marker_df       marker_MS       error_df        error_MS        model_df        model_MS        minorObs
#M216    S1_1742 1       1742    0.57105 0.56565 0.00444 0.24878 0.61836 0.89342 0.34544 2       27.73206        256     48.56339        2       27.73206        53


### for Manhattan_inpu for  R
#
#	SNP CHR BP         P
# rs1   1  1 0.9148060
# rs2   1  2 0.9370754
#

	my $trait="NA";
	my $tpfname="";

	
	$lineIN=<IN>;
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp=split(/\t/,$lineIN);
		if($tp[0] eq $trait){
			if($tp[1] ne "None"){
				print OUTF "$tp[1]\t$tp[2]\t$tp[3]\t$tp[5]\n";
			}
		}else{
			if($trait ne "NA"){
				close(OUTF);
			}
			#creat new file
			$trait = $tp[0];
			$tpfname="./manhattan/$trait.txt";
			open(OUTF,'>',"$tpfname");	
			print OUTF "\tSNP\tCHR\tBP\tP\n";
		}
	}
	close(OUTF);
	print "done\n";
}


sub TASSEL_mlm_2_Manhattan_input{

### for mlm_?_9.txt
#
#Trait   Marker  Chr     Pos     df      F       p       add_effect      add_F   add_p   dom_effect      dom_F   dom_p   errordf MarkerR2        Genetic Var     Residual Var    -2LnLikelihood
#M168    None                    0       NaN     NaN     NaN     NaN     NaN     NaN     NaN     NaN     259     NaN     0.01082 0.00359 -2.9294E2
#M168    S1_1742 1       1742    2       0.06917 0.93318 0.01785 0.07226 0.78829 0.01251 0.03271 0.85662 259     5.3556E-4       0.01082 0.00359 -2.9294E2


### for Manhattan_inpu for  R
#
#	SNP CHR BP         P
# rs1   1  1 0.9148060
# rs2   1  2 0.9370754
#

	my $trait="NA";
	my $tpfname="";

	
	$lineIN=<IN>;
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp=split(/\t/,$lineIN);
		if($tp[0] eq $trait){
			if($tp[1] ne "None"){
				print OUTF "$tp[1]\t$tp[2]\t$tp[3]\t$tp[6]\n";
			}
		}else{
			if($trait ne "NA"){
				close(OUTF);
			}
			#creat new file
			$trait = $tp[0];
			$tpfname="./manhattan/$trait.txt";
			open(OUTF,'>',"$tpfname");	
			print OUTF "\tSNP\tCHR\tBP\tP\n";
			print OUT "$trait\n";
		}
	}
	close(OUTF);
	print "done\n";
}




sub Extract_KEGG_compound_info{

#	ENTRY       C00001                      Compound
#	NAME        H2O;
#	            Water
#	FORMULA     H2O
#	EXACT_MASS  18.0106
#	MOL_WEIGHT  18.0153

	my $outt=$ARGV[0].".list";
	open(OUTT,'>',"$outt");
	
	my $ENTRY="";
	my @NAME;		#start from 12
	my $count=0;
	my $FORMULA="";
	my $MOL_WEIGHT=0;
	my $add=0;
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($tp[0] eq "ENTRY"){
			#output
			if($count>0){			
				print OUT "$ENTRY\t@NAME\t$FORMULA\t$MOL_WEIGHT\n";
				for(my $i=0;$i<$count;$i++){
					print OUTT "$NAME[$i]\t$ENTRY\t$FORMULA\t$MOL_WEIGHT\n";
				}
			}
			#add new
			$add=0;		
			$ENTRY=$tp[1];
			@NAME=();
			$count=0;
		}elsif($tp[0] eq "NAME"){
			$add=1;
			$NAME[$count]=substr($lineIN,12,length($lineIN));	
			#print "$NAME[$count]\n";
			$count++;
		}elsif($tp[0] eq "FORMULA"){
			$add=0;	
			$FORMULA=$tp[1];	
		}elsif($tp[0] eq "MOL_WEIGHT"){
			$add=0;	
			$MOL_WEIGHT=$tp[1];	
		}elsif($add==1){
			$NAME[$count]=substr($lineIN,12,length($lineIN));
			#print "$NAME[$count]\n";
			$count++;
		}else{
		}		
	}
	print OUT "$ENTRY\t@NAME\t$FORMULA\t$MOL_WEIGHT";
	for(my $i=0;$i<$count;$i++){
		print OUTT "$NAME[$i]\t$ENTRY\t$FORMULA\t$MOL_WEIGHT\n";
	}
	print "\ndone\n";
	close(OUTT);
}


sub vcf_2_MatrixEQTL_SNP{

# 0/0	0
# 0/1	1
# 1/1	2

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  C100    C101 
	

	my $outt=$ARGV[0].".snploc";
	open(OUTT,'>',"$outt");
	print OUTT "snpid\tchr\tpos\n";
	
	my $head=1;
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp=split(/\t/,$lineIN);
		my $tp=@tp;
		if($head==1){
			if($lineIN=~/CHROM/){
				#head
				print OUT "snpid";
				for(my $i=9;$i<$tp;$i++){
					print OUT "\t$tp[$i]";
				}
				print OUT "\n";
				$head=0;
			}	
		}else{
			#body
			my $name="S$tp[0]_$tp[1]";
			print OUT "$name";		
			for(my $i=9;$i<$tp;$i++){
				if($tp[$i] eq "0/0"){
					print OUT "\t0";
				}elsif($tp[$i] eq "1/1"){
					print OUT "\t2";
				}else{
					print OUT "\t1";
				}
			}
			print OUT "\n";	
			#loc	
			print OUTT "$name\t$tp[0]\t$tp[1]\n";	
		}				
	}

	close(OUTT);
	print "done\n";
}




sub TASSEL_distance_2_phylip{
#c2b04056    0.000000  0.000010  0.000010  0.000010  0.000010  0.000010

	my $first=1;
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/#/){
			#print "$lineIN\n";
		}else{
			if($first==1){
				print OUT "$lineIN\n";
				print "$lineIN samples\n"; 
				$first=0;
			}else{
				my @tp=&Getheadinfo($lineIN);
				my $tp=@tp;
				my $count=0;
				for(my $i=0; $i<$tp; $i++){
					$count++;
					if($count>7){
						$count=1;
						print OUT "\n";
					}
					my $aa=	$tp[$i];
					my $len=length($aa);		
					if($i==0){						
						if($len>9){
							print "Length of $tp[0] > 9\n";
							last;
						}else{
							for(my $j=$len; $j<10;$j++){
								$aa=$aa." ";
							}							
						}
					}else{
						if($len>8){
							$aa="  ".substr($aa,0,8);
						}else{
							for(my $j=$len; $j<8;$j++){
								$aa=$aa."0";
							}	
							$aa="  ".$aa;						
						}						
					}
					print OUT $aa;					
				}
				print OUT "\n";
			}
		}		
	}	
	print "done\n";
}


sub MakerAnnotationRenameConn{


#Chr02	maker	gene	23184331	23187775	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4;Name=maker-Chr02-exonerate_est2genome-gene-231.4
#Chr02	maker	mRNA	23184331	23187775	2957	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1;Parent=maker-Chr02-exonerate_est2genome-gene-231.4;Name=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1;_AED=0.02;_eAED=0.02;_QI=145|1|1|1|0|0|5|1219|530
#Chr02	maker	mRNA	23184331	23187775	15157	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2;Parent=maker-Chr02-exonerate_est2genome-gene-231.4;Name=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2;_AED=0.00;_eAED=0.00;_QI=145|1|1|1|0|0|4|1325|530
#Chr02	maker	exon	23185410	23185728	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2:exon:12645;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	exon	23184331	23184985	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:exon:12644;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1,maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	exon	23185169	23185303	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:exon:12643;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1,maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	exon	23185410	23185538	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:exon:12642;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	exon	23185645	23185728	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:exon:12641;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	exon	23185822	23187775	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:exon:12640;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1,maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	five_prime_UTR	23187631	23187775	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:five_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	CDS	23186038	23187630	.	-	0	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:cds;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	three_prime_UTR	23185822	23186037	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	three_prime_UTR	23185645	23185728	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	three_prime_UTR	23185410	23185538	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	three_prime_UTR	23185169	23185303	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	three_prime_UTR	23184331	23184985	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-1
#Chr02	maker	five_prime_UTR	23187631	23187775	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2:five_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	CDS	23186038	23187630	.	-	0	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2:cds;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	three_prime_UTR	23185822	23186037	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	three_prime_UTR	23185410	23185728	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	three_prime_UTR	23185169	23185303	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2
#Chr02	maker	three_prime_UTR	23184331	23184985	.	-	.	ID=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2:three_prime_utr;Parent=maker-Chr02-exonerate_est2genome-gene-231.4-mRNA-2


	
	$sp=$_[1];
	$count=0;
	$gID="";
	@tID;	#0-old name; 1-new name ; 2-D
	$tcount=0;
	@set;	#the number for 0-exon,1-u5,2-CDS,3-u3 for each mRNA; 2-D
	@mRNA;
	@exon;
	@CDS;
	@u5;	
	@u3;			
	$direct;

	my $ot=$ARGV[0].".gene";
	open(OUTG,'>',"$ot");
	my $ot=$ARGV[0].".mRNA";
	open(OUTD,'>',"$ot");	
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN,"=",";");
		my $tp=@tp;
		my $o="$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]";	
			
		if($tp[2] eq "contig"){
			print OUT "$lineIN\n";
		}
		elsif($tp[2] eq "gene"){		
			#output last one if had
			if($tcount>0){
				&printmRNAs();
			}
			
			#new gene	
			$count++;
			$tcount=0;
			$gID=$sp.substr($tp[0],3,2)."g";
			my $tpl=length($count);
			for(my $i=$tpl;$i<6;$i++){
				$gID=$gID."0";
			}
			$gID=$gID.$count;
			$direct=$tp[6];
			#print "$gID\n";
			print OUT "$o\tID=$gID;Name=$gID\n";
			print OUTG "$tp[9]\t$gID\n";						
		}
		elsif($tp[2] eq "mRNA"){
			$tcount++;
			$tID[$tcount-1][0]=$tp[9];	#old name
			$tID[$tcount-1][1]=$gID.".".$tcount;	#new name
			$mRNA[$tcount-1]= $o;
			print OUTD "$tp[9]\t$tID[$tcount-1][1]\n";
			$set[$tcount-1]=[0,0,0,0];	#0-exon,1-u5,2-CDS,3-u3
		} 
		elsif($tp[2] eq "exon"){
			for(my $i=0;$i<$tcount;$i++){
				if(index($lineIN,$tID[$i][0],)>0){
					#add a head for sort
					my $head="";
					for(my $j=length($tp[3]);$j<10;$j++){
						$head=$head."0";
					}
					$head=$head.$tp[3];
					#add exon
					$exon[$i][$set[$i][0]]= $head.$o;
					#print "$set[$i][0]\t$exon[$i][$set[$i][0]]\n";
					$set[$i][0]++;
				}
			}		
		}
		elsif($tp[2] eq "five_prime_UTR"){
			for(my $i=0;$i<$tcount;$i++){
				if(index($lineIN,$tID[$i][0],)>0){
					$u5[$i][$set[$i][1]]= $o;
					$set[$i][1]++;
				}
			}	 	
		}				 	
		elsif($tp[2] eq "CDS"){
			for(my $i=0;$i<$tcount;$i++){
				if(index($lineIN,$tID[$i][0],)>0){
					$CDS[$i][$set[$i][2]]= $o;
					$set[$i][2]++;
				}
			}	  	
		}
		elsif($tp[2] eq "three_prime_UTR"){
			for(my $i=0;$i<$tcount;$i++){
				if(index($lineIN,$tID[$i][0],)>0){
					$u3[$i][$set[$i][3]]= $o;
					$set[$i][3]++;
				}
			}   	
		}else{
			print "$lineIN\n";
		}	
		
	}
	&printmRNAs();
	close(OUTD);
	close(OUTG);
	print "done\n";
}

sub printmRNAs{
	#print out
	for(my $i=0;$i<$tcount;$i++){
		#mRNA
		print OUT $mRNA[$i]."\tID=$tID[$i][1];Parent=$gID;Name=$tID[$i][1]\n";
		
		#exon sort and out put
		my @tpexon;
		for(my $j=0;$j<$set[$i][0];$j++){
			$tpexon[$j]=$exon[$i][$j]
		}
			
		@tpexon=sort(@tpexon);

		if($direct eq "+"){
			for(my $j=0;$j<$set[$i][0];$j++){
				my $tpid=$j+1;	
				my $tpo=substr($tpexon[$j],10,length($tpexon[$j]));	
				print OUT "$tpo\tID=$tID[$i][1]:exon_$tpid;Parent=$tID[$i][1]\n";
			}
		}else{
			for(my $j=$set[$i][0]-1;$j>=0;$j--){
				my $tpid=$set[$i][0]-$j;
				my $tpo=substr($tpexon[$j],10,length($tpexon[$j]));	
				#print "$tpo\n";
				print OUT "$tpo\tID=$tID[$i][1]:exon_$tpid;Parent=$tID[$i][1]\n";
			}			
		}
		
		#utr5
		for(my $j=0;$j<$set[$i][1];$j++){
			my $tpid=$j+1;
			print OUT $u5[$i][$j]."\tID=$tID[$i][1]:utr5_$tpid;Parent=$tID[$i][1]\n";
		}	
			
		#CDS
		for(my $j=0;$j<$set[$i][2];$j++){
			my $tpid=$j+1;
			print OUT $CDS[$i][$j]."\tID=$tID[$i][1]:CDS_$tpid;Parent=$tID[$i][1]\n";
		}	
			
		#utr3
		for(my $j=0;$j<$set[$i][3];$j++){
			my $tpid=$j+1;
			print OUT $u3[$i][$j]."\tID=$tID[$i][1]:utr3_$tpid;Parent=$tID[$i][1]\n";
		}		
	}
	
	#clean
	$tcount=0;	
	@tID=();
	@set=();		
	@mRNA=();
	@exon=();
	@CDS=();
	@u5=();
	@u3=();		
	
}

sub MakerAnnotationRename{

#Chr00   maker   gene    10058   18557   .       -       .       ID=maker-Chr00-snap-gene-0.43;Name=maker-Chr00-snap-gene-0.43
#Chr00   maker   mRNA    10058   18557   .       -       .       ID=maker-Chr00-snap-gene-0.43-mRNA-1;Parent=maker-Chr00-snap-gene-0.43;Name=maker-Chr00-snap-gene-0.43-mRNA-1;_AED=0.02;_eAED=0.02;_QI=0|0|0|1|0.5|0.33|3|0|2060
#Chr00   maker   exon    16476   18557   .       -       .       ID=maker-Chr00-snap-gene-0.43-mRNA-1:exon:6;Parent=maker-Chr00-snap-gene-0.43-mRNA-1
#Chr00   maker   CDS     16476   18557   .       -       0       ID=maker-Chr00-snap-gene-0.43-mRNA-1:cds;Parent=maker-Chr00-snap-gene-0.43-mRNA-1
#Chr12	maker	three_prime_UTR	11599103	11599421	.	+	.	ID=maker-Chr12-augustus-gene-116.60-mRNA-1:three_prime_utr;Parent=maker-Chr12-augustus-gene-116.60-mRNA-1
	
	my $sp=$_[1];
	my $count=0;
	my $gID="";
	my $tID="";
	my $tcount=0;
	my @exon;
	my $CDScount=0;
	my $exoncount=0;
	my $UTRcount=0;
	my $direct;

	my $ot=$ARGV[0].".rename";
	open(OUTD,'>',"$ot");
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN,"=",";");
		my $tp=@tp;
		my $o="$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]";
		#print $lineIN."\n";
		if($tp[2] eq "exon"){			
			$exon[$exoncount]=[@tp];
			$exoncount++;
		}else{
			#exon output
			if($exoncount>0){		
				
				if($direct eq "+"){					
					for (my $i=0;$i<$exoncount;$i++){
						my $o="$exon[$i][0]\t$exon[$i][1]\t$exon[$i][2]\t$exon[$i][3]\t$exon[$i][4]\t$exon[$i][5]\t$exon[$i][6]\t$exon[$i][7]";
						$tpno=$i+1;
						print OUT $o."\tID=$tID:exon_$tpno;Parent=$tID\n"
					}
				}
				else{  # -		
					my $dis=1;			
					if($exoncount>1){
						$dis= $exon[1][3]-$exon[0][3];
						if ($dis<0 && $direct eq "-"){
							#print "$gID $direct $dis\n";
						}					
					}
					if($dis >0){
						for (my $i=$exoncount-1;$i>=0;$i--){
							my $o="$exon[$i][0]\t$exon[$i][1]\t$exon[$i][2]\t$exon[$i][3]\t$exon[$i][4]\t$exon[$i][5]\t$exon[$i][6]\t$exon[$i][7]";
							my $tpno=$exoncount-$i;
							print OUT $o."\tID=$tID:exon_$tpno;Parent=$tID\n"
						}
					}else{
						for (my $i=0;$i<$exoncount;$i++){
							my $o="$exon[$i][0]\t$exon[$i][1]\t$exon[$i][2]\t$exon[$i][3]\t$exon[$i][4]\t$exon[$i][5]\t$exon[$i][6]\t$exon[$i][7]";
							my $tpno=$i+1;
							print OUT $o."\tID=$tID:exon_$tpno;Parent=$tID\n"
						}					
					}
				}				
			}
			#clean
			$exoncount=0;
			@exon=();
			
			
			if($tp[2] eq "contig"){
				print OUT "$lineIN\n";
			}
			elsif($tp[2] eq "gene"){					
					$count++;
					$tcount=0;
					$gID=$sp.substr($tp[0],3,2)."g";
					my $tpl=length($count);
					for(my $i=$tpl;$i<6;$i++){
						$gID=$gID."0";
					}
					$gID=$gID.$count;
					$direct=$tp[6];
					#print "$gID\n";
					print OUT "$o\tID=$gID;Name=$gID\n";
					#print OUTD "$tp[9]\t$gID";
			}
			elsif($tp[2] eq "mRNA"){
				$tcount++;
				$CDScount=0;
				$exoncount=0;
				$UTRcount=0;			
				$tID=$gID.".".$tcount;
				print OUT "$o\tID=$tID;Parent=$gID;Name=$tID\n";
				print OUTD "$tp[9]\t$tID\n";
			}  	
			elsif($tp[2] eq "CDS"){
					$CDScount++;
					print OUT "$o\tID=$tID:CDS_$CDScount;Parent=$tID\n";    	
			}
			elsif($tp[2] eq "five_prime_UTR" || $tp[2] eq "three_prime_UTR"){
					$UTRcount++;
					print OUT "$o\tID=$tID:utr_$UTRcount;Parent=$tID\n";     	
			}else{
				print "$lineIN\n";
			}	
		}			
	}
	close(OUTD);
	print "done\n";
}


sub Batch_PAML_codmel_mlc_Read{
	#M0 branch model_A site
	
    #input 
	#kappa (ts/tv) = 23.12459
	
	#lnL(ntime: 34  np: 36):   -756.985546      +0.000000
	

	my $model=$_[1];
	print "Model = $model\n";
	if($model eq "branch"){
		#branch
		#w (dN/dS) for branches:  0.10108 1.53449 2.06752	
		print "GeneID	lnL	ts/tv	w1	w2	w3\n";		
		print OUT "GeneID	lnL	ts/tv	w1	w2	w3\n";
	}elsif($model eq "model_A"){
		#model_A
		#MLEs of dN/dS (w) for site classes (K=4)
		#	
		#site class             0        1       2a       2b
		#proportion       0.89157  0.08624  0.02023  0.00196
		#background w     0.00000  1.00000  0.00000  1.00000
		#foreground w     0.00000  1.00000 43.89416 43.89416		
		print "GeneID	lnL	ts/tv	p_0	p_1	p_2a	p_2b	w_0	w_1	w_2a	w_2b\n";
		print OUT "GeneID	lnL	ts/tv	p_0	p_1	p_2a	p_2b	w_0	w_1	w_2a	w_2b\n";
		my $outfname=$ARGV[0].".BEB";
		open(OUTF,'>',"$outfname");				
	}elsif($model eq "M0"){
		#M0
		#omega (dN/dS) =  0.36753	
		print "GeneID	lnL	ts/tv	dN/dS\n";		
		print OUT "GeneID	lnL	ts/tv	dN/dS\n";
	}elsif($model eq "site"){
		#site  
		#here two model in one mlc file
		#Model 1: NearlyNeutral (2 categories)
		#MLEs of dN/dS (w) for site classes (K=2)
		#
		#p:   0.92023  0.07977
		#w:   0.00000  1.00000
		
		#Model 2: PositiveSelection
		#MLEs of dN/dS (w) for site classes (K=3)
		#
		#p:   0.92660  0.05810  0.01530
		#w:   0.00000  1.00000 25.94057	
		print "site\n";
	}
	my $count=0;
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my $fname=$lineIN.".mlc";
		my $out=$lineIN;
		print ".";
		open(TPIN,'<',"$fname");		
		my $lineINF="";	
		my $lnL    ="NA";
		my $kappa  ="NA";	
		my $w	   ="NA";
		$count++;
		
		my $state=1;
		if(length($_[2])>0){
			$state=0;
		}		
		while($lineINF=<TPIN>){
			if($lineINF=~/$_[2]/){
				$state=1;
			}
			if($state==1){
				chomp($lineINF);
		    	
				if($lineINF=~/lnL\(ntime/){
					my @tp = &Getheadinfo($lineINF);
					$lnL=$tp[4];
				}			
				if($lineINF=~/kappa/){
					my @tp = &Getheadinfo($lineINF);
					$kappa=$tp[3];
				}
				if($lineINF=~/branches:/){
					my @tp = &Getheadinfo($lineINF);
					$w="$tp[4]\t$tp[5]\t$tp[6]";
					last;
				}
				if($lineINF=~/omega/){
					my @tp = &Getheadinfo($lineINF);
					$w="$tp[3]";
					last;
				}			
				
				if($lineINF=~/MLEs/){
					$lineINF=<TPIN>;
					$lineINF=<TPIN>;
					if($model eq "site"){
						chomp($lineINF);
						my @tp = &Getheadinfo($lineINF);
						#print "@tp\n";
						$w="$tp[1]\t$tp[2]\t$tp[3]";
						$lineINF=<TPIN>;
						chomp($lineINF);
						@tp = &Getheadinfo($lineINF);
						$w=$w."\t$tp[1]\t$tp[2]\t$tp[3]";
						#last;
					}elsif($model eq "model_A"){
						$lineINF=<TPIN>;
						chomp($lineINF);				
						my @tp = &Getheadinfo($lineINF);
						$w="$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]";
						$lineINF=<TPIN>;
						$lineINF=<TPIN>;
						chomp($lineINF);
						@tp = &Getheadinfo($lineINF);
						$w=$w."\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]";
						#last;
					}
				}
				
				if($lineINF=~/Nielsen/ && $model eq "model_A"){
				#if($lineINF=~/Nielsen/){
					print "+";
					#Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
					#Positive sites for foreground lineages Prob(w>1):
					#    67 D 0.955*
					#   136 R 0.945
					#   137 R 0.944						
					$lineINF=<TPIN>;
					while($lineINF=<TPIN>){
						chomp($lineINF);
						my @tp = &Getheadinfo($lineINF);
						my $tp=@tp;
						if($tp<3){
							last;
						}	
						#read and write
						print OUTF "$lineIN\t$tp[0]\t$tp[1]\t$tp[2]\n";
					}					
					last;								
				}
				
			}	
		}
		print OUT "$lineIN\t$lnL\t$kappa\t$w\n";	
	}
	close(TPIN);
	if($model eq "model_A"){
		close(OUTF);
	}	
	print "done\n";
}


sub RemoveSameCodeforAlignedFas{
    #input must be aligned
	my $n=3;  #remove unit
	if($_[1]>0){
		$n=$_[1];
	}
	
	#readin
	my @name=();
	my @seq=();
	my $count=-1;
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/>/){
			$count++;
			$name[$count]=$lineIN;
			$seq[$count]="";
		}else{
			$seq[$count]=$seq[$count].$lineIN;
		}
	}
	$count++;
	
	#cut
	my @cut=();
	for(my $j=0;$j<$count;$j++){
		$cut[$j]="";
	}	
	my $len=length($seq[0]);
	for(my $i=0;$i<$len;$i+=$n){
		my $add=0;
		my $char=substr($seq[0],$i,$n);
		for(my $j=1;$j<$count;$j++){
			if($char ne substr($seq[$j],$i,$n)){
				$add=1;
				last;
			}
		}
		if($add>0){
			for(my $j=0;$j<$count;$j++){
				$cut[$j]=$cut[$j].substr($seq[$j],$i,$n)
			}
		}			
	}

	#out
	for(my $j=0;$j<$count;$j++){
		print OUT "$name[$j]\n$cut[$j]\n";
	}	
	print "done\n";
}


sub CutAlignFasGapForProteinCoding{
      
	my @seq=();
	my @name=();
	my $count=-1;
	
	#read in seq
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp = &Getheadinfo($lineIN);
			$count++;
			$name[$count]=$tp[0];
			$seq[$count]="";
		}else{
			$seq[$count]=$seq[$count].$lineIN;
		}
	}
	$count++;
	print "$count seq readed\n";
	
	#cut
	my @out=();
	for(my $i=0;$i<$count;$i++){
		$out[$i]="";
	}	
	my $len=length($seq[0]);	
	for(my $k=0;$k<$len;$k+=3){
		my $add=1;
		for(my $i=0;$i<$count;$i++){
			my $tpseq=substr($seq[$i],$k,3);
			if($tpseq=~/-/){
				$add=0;
				last;
			}
		}
		#add
		if($add==1){
			for(my $i=0;$i<$count;$i++){
				my $tpseq=substr($seq[$i],$k,3);
				$out[$i]=$out[$i].$tpseq;
			}			
		}
	}
	#out
	for(my $i=0;$i<$count;$i++){
		print OUT "$name[$i]\n$out[$i]\n";
	}	
	print "done\n";
}


sub MutationCountForCode123{
     
	my @seq=();
	my @name=();
	my $count=-1;
	my @m=();	#m[0]-code 1 m[1]-code 2 m[2]-code 3
	
	my $seed="";
	
	#read in seq
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp = &Getheadinfo($lineIN);
			$count++;
			$name[$count]=$tp[0];
			$seq[$count]="";
		}else{
			$seq[$count]=$seq[$count].$lineIN;
		}
	}
	$count++;
	print "$count seq readed\n";
	#find seed
	my $seed=0;
	if($_[1]){
		for(my $i=0;$i<$count;$i++){
			if($name[$i]=~/$_[1]/){
				$seed=$i;
				print "seed is $name[$i]\n";
			}
		}
	}
	
	#check	
	my @out=();
	for(my $i=0;$i<$count;$i++){
		$m[0][$i]=0;
		$m[1][$i]=0;
		$m[2][$i]=0;
	}	
	my $len=length($seq[0]);
	for(my $k=0;$k<$len;$k+=3){
		for(my $j=0;$j<3;$j++){
			my $seedchr=substr($seq[$seed],($k+$j),1);
			for(my $i=0;$i<$count;$i++){
				if(substr($seq[$i],($k+$j),1) ne $seedchr){
					$m[$j][$i]++;
				}
			}
		}
	}
	#out
	for(my $i=0;$i<$count;$i++){
		print OUT "$name[$i]\t$m[0][$i]\t$m[1][$i]\t$m[2][$i]\t$len\n";
	}	
	print "done\n";
}


sub CutAlignFasbySeedSymble{
      
	my $seedsymble=$_[1];

	my @seq=();
	my @name=();
	my $count=-1;
	
	my $seed="";
	
	#read in seq
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp = &Getheadinfo($lineIN);
			$count++;
			$name[$count]=$tp[0];
			$seq[$count]="";
		}else{
			$seq[$count]=$seq[$count].$lineIN;
		}
	}
	$count++;
	print "$count seq readed\n";
	
	#find seed
	my $len=length($seq[0]);
	for(my $i=0;$i<$count;$i++){
		if($name[$i]=~/$seedsymble/){
			$seed=$i;
		}
	}

	my $loc=index($seq[$seed],"-",);
	while($loc>=0){		
		print ".";		
		for(my $i=0;$i<$count;$i++){
			my $tpseq=substr($seq[$i],0,$loc).substr($seq[$i],($loc+1),($len-$loc+1));
			$seq[$i]=$tpseq;
		}
		$loc=index($seq[$seed],"-",);
	}
	#out
	for(my $i=0;$i<$count;$i++){
		print OUT "$name[$i]\n$seq[$i]\n";
	}	
	print "done\n";
}



sub Blat_2_loc{
      
	#input
	#match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
	#     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------
	#992	4	0	0	0	0	4	237	-	XM_002372115.1	996	0	996	JZDT01000032.1	16086	7215	8448	5	238,85,91,112,470,	0,238,323,414,526,	7215,7517,7653,7789,7978,
   	  	
   	#out
   	#>filename=ID|exonno|chrom:complement(join(2354845..2355351,2355357..2355806))
   	
   	my $filename=substr($ARGV[0],0,index($ARGV[0],".",));
   	print "$filename\n";

	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp         = &Getheadinfo($lineIN);
		my @blockSizes = &Getheadinfo($tp[18],",");
		my @tStarts    = &Getheadinfo($tp[20],",");
				
		my $out=">$filename=$tp[9]|$tp[17]|$tp[13]:";
		if($tp[8] eq "-"){
			$out=$out."complement(";
		}
		
		if($tp[17]==1){
			my $ed=$tStarts[0]+$blockSizes[0]-1;
			#print "ed=$ed\n";
			$out=$out."$tStarts[0]..$ed";
		}else{
			my $ed=$tStarts[0]+$blockSizes[0]-1;
			#print "ed=$ed\n";
			$out=$out."join($tStarts[0]..$ed";
			for(my $i=1;$i<$tp[17];$i++){
				$ed=$tStarts[$i]+$blockSizes[$i]-1;
				$out=$out.",$tStarts[$i]..$ed";
			}
			$out=$out.")";
		}

		if($tp[8] eq "-"){
			$out=$out.")";
		}		
		print OUT "$out\n";
	}
	
	print "done\n";
}



sub CDHitClstrInfo{
      
	#input
#	>Cluster 0
#	0	6667nt, >3=B2N_0844... *
#	1	6609nt, >5=B2N_0844... at +/97.28%
#	>Cluster 1
#	0	334nt, >3=B4N_0765... at -/88.02%
#	1	162nt, >3=B4N_0828... at +/85.80%
#	2	853nt, >3=B4N_0928... at -/83.00%
#	3	4234nt, >3=B4N_1149... at +/88.14%
#	4	333nt, >5=B4N_0765... at -/88.29%
#	5	161nt, >5=B4N_0828... at +/89.44%
#	6	851nt, >5=B4N_0928... at -/82.61%
#	7	4236nt, >5=B4N_1149... *


   	
   	my $cluster=-1;
   	my @detail=();
   	my @mark=();
	while($lineIN=<IN>){	
		chomp($lineIN);
		if(substr($lineIN,0,1) eq ">"){
			$cluster++;
			$detail[$cluster]=" ";
		}else{
			my @tp= &Getheadinfo($lineIN);
			if($tp[3] eq "*"){
				$mark[$cluster]=substr($tp[2],1,(length($tp[2])-4))."\t".$tp[1];		
			}else{
				$detail[$cluster]=$detail[$cluster]." ".substr($tp[2],1,(length($tp[2])-4));
			}
		}
	}
	
	for(my $i=0;$i<$cluster;$i++){
		print OUT "$mark[$i]$detail[$i]\n";
	}

	print "done\n";
}



sub Add_ID_for_vcf{
      
	#input
#	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  C100    C101    C102    C103
#	1       1742    .       C       T       .       .       PR      GT      0/0     0/0     0/1     0/0 
 
#OUTPUT ID is S1_1742 
   	
	my $add=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($add == 1){
			#add ID
			my @tp = split(/\t/,$lineIN);
			my $tp=@tp;
			my $ID="S$tp[0]_$tp[1]";
			print OUT "$tp[0]\t$tp[1]\t$ID";
			for(my $i=3;$i<$tp;$i++){
				print OUT "\t$tp[$i]";
			}
			print OUT "\n";
		}else{
			print OUT "$lineIN\n";
		}
		
		if($lineIN=~/#CHROM/){
			$add=1;
		}
		
	}

	print "done\n";
}


sub AddBlankToLine{
      
	#input
	#aaaaa	aaaaa
	#aaa aaaaaa
	#output
	#aaaaa	aaaaa
	#aaa    aaaaaa	
   	
   	my @linelength=();
   	
   	my @array=();
   	my $count=0;
   		
   	$lineIN=<IN>;
   	chomp($lineIN);
   	my @tp= &Getheadinfo($lineIN);
   	$array[$count]=[@tp];
   	$count++;
 
   	my $lineNo=@tp;
	for(my $i=0;$i<$lineNo;$i++){
		$linelength[$i]=length($tp[$i]);
	}   	
   	
	while($lineIN=<IN>){	
		chomp($lineIN);
		@tp= &Getheadinfo($lineIN);
		$array[$count]=[@tp];
		
		for(my $i=0;$i<$lineNo;$i++){
			my $tplen=length($tp[$i]);
			if($linelength[$i]<$tplen){
				$linelength[$i]=$tplen;
				print ".";
			}
		}	
		$count++;
	}

	print "Line No = $lineNo\t@linelength\nCount = $count\n";
	#output
	for(my $i=0;$i<$count;$i++){
		my $o="";
		for(my $j=0;$j<($lineNo-1);$j++){
			$o=$o.$array[$i][$j];
			my $addlen=$linelength[$j]-length($array[$i][$j])+1;
			for(my $k=0;$k<=$addlen ;$k++){
				$o=$o." ";
			}
		}
		$o=$o.$array[$i][$lineNo-1];
		print OUT "$o\n";
	}
	print "done\n";
}


sub AddfasTypeforRepeatMasker{
      
	#input
	#>OL005_JAP_0000
	#out
	#>OL005_JAP_0000#LTR/OL005
   	
   	my $first=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp= &Getheadinfo($lineIN, "_");
			my $type=substr($tp[0],1,length($tp[0]));
			$lineIN=$lineIN."#LTR/$type";
		}
		print OUT "$lineIN\n";
		
	}

	print "done\n";
}

sub venn_state_by_filelist{
      
	#input
	#list file 1
	#list file 2
   	
   	my $head="ID";
   	my @data=();
   	my %ID=();		#ID hash
   	my $count=0;	#@data length
   	my $fno=0;  #read in file no
	while($lineIN=<IN>){	
		chomp($lineIN);
		
		#read in file
		open(INF,'<',"$lineIN");	
		my $tplineIN;
		my %tpdata=();
		while($tplineIN=<INF>){
			chomp($tplineIN);
			$tpdata{$tplineIN}=1;
			#cheak hav or not
			if($ID{$tplineIN}!=1){
				$ID{$tplineIN}=1;
				$data[$count][0]=$tplineIN;
				my $add="";
				for(my $i=0;$i<$fno;$i++){
					$add=$add."\t0";
				}
				$data[$count][1]=$add;
				$count++;
			}
		}	
		close(INF);
		$fno++;
		
		#add name
		$head=$head."\t$lineIN";
		#add state		
		for(my $i=0;$i<$count;$i++){
			if($tpdata{$data[$i][0]}==1){
				$data[$i][1]=$data[$i][1]."\t1";
			}else{
				$data[$i][1]=$data[$i][1]."\t0";
			}
		}		
	}

	#print out
	print OUT "$head\n";
	for(my $i=0;$i<$count;$i++){
		print OUT "$data[$i][0]$data[$i][1]\n";
	}

	print "done\n";
}



sub count_3_piont_trend{
	#if 1 and 3 had 2 fold change, then 1,2 and 2,3 only need 1.5 fold change
      
	#input
	#0-name	1-1 2-2 3-3
   	#out: 0=nuchange, 1-up, 2=down
   	
   	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $c1v3=1;
		if($tp[1]==0){
			$c1v3=$tp[3];			
		}elsif($tp[3]==0){	$c1v3=1/$tp[1];
		}else{				$c1v3=$tp[3]/$tp[1];
		}
		#1v2
		my $c1v2=1;
		if($tp[1]==0){		$c1v2=$tp[2];
		}elsif($tp[2]==0){	$c1v2=1/$tp[1];
		}else{				$c1v2=$tp[2]/$tp[1];
		}
			
		my $o1=0;
		if($c1v2>=2){		$o1=1;}	
		if($c1v2<=0.5){		$o1=2;}
		
		#2v3		
		my $c2v3=1;
		if($tp[2]==0){		$c2v3=$tp[3];		
		}elsif($tp[3]==0){	$c2v3=1/$tp[2];
		}else{				$c2v3=$tp[3]/$tp[2];
		}
			
		my $o2=0;
		if($c2v3>=2){		$o2=1;}	
		if($c2v3<=0.5){		$o2=2;}		
		
		#compair
		if($c1v3>2){
			if($c1v2>1 && $c2v3 >1){
				if($c1v2>=1.5){		$o1=1;}	
				if($c2v3>=1.5){		$o2=1;}							
			}
		}
		if($c1v3<0.5){
			if($c1v2<1 && $c2v3 <1){
				if($c1v2<=0.66){	$o1=2;}	
				if($c2v3<=0.66){	$o2=2;}							
			}
		}	
		print OUT "$lineIN\t$o1$o2\n";
	}

	print "done\n";
}

sub csv_2_tsv{
      
	#input
	#"","baseMean","log2FoldChange","lfcMLE","lfcSE","stat","pvalue","padj"
	#"XLOC_000001",0.655621431295059,0.17959214042933,0.830789744811771,0.570941635816685,0.314554289200574,0.753100090413155,NA
	#output
   	
   	my $first=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		$lineIN=~ tr/\"//d;
		my @tp=split(/\,/,$lineIN);
		my $tp=@tp;
		if($first==0){
			$tp[0]="ID";
			$first=1;
		}
		print OUT "$tp[0]";
		for(my $i=1;$i<$tp;$i++){
			print OUT "\t$tp[$i]";
		}
		print OUT "\n";
	}

	print "done\n";
}



sub ThreeSamplesUpDownCompare{
      
	#input
	#0-ID	1-CK_24.avg	2-CK_24.sd	3-CK_24.sd%	4-CK_48.avg	5-CK_48.sd	6-CK_48.sd%	7-CK_72.avg	8-CK_72.sd	9-CK_72.sd%
  	#outpput
   	
   	
   	$lineIN=<IN>;
   	chomp($lineIN);
   	my @tp= &Getheadinfo($lineIN);
   	print OUT "$tp[0]\t$tp[1]\t$tp[4]\t$tp[7]\tAvg\tlog2(2v1)\tlog2(3v2)\tCode\tSig_Code\tlog(1)\tlog(2)\tlog(3)\n";
   	
	while($lineIN=<IN>){	
		chomp($lineIN);
		@tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		$tp[1]=$tp[1]>0.1?$tp[1]:0.1;
		$tp[4]=$tp[4]>0.1?$tp[4]:0.1;
		$tp[7]=$tp[7]>0.1?$tp[7]:0.1;
		my $avg=($tp[1]+$tp[4]+$tp[7])/3;
		my $max=$tp[1]>$tp[4]?$tp[1]:$tp[4];
		$max=$max>$tp[7]?$max:$tp[7];
		my $log2v1=0;
		my $log3v2=0;
		my $code1="=";
		my $scode1="=";	
		my $code2="=";			
		my $scode2="=";	
		my $p1=0;			
		my $p2=0;
		my $p3=0;
		
		if($tp[1]>0){
			$p1=log($tp[1])/log(2);
		}else{
			$p1=-1;
		}
		if($tp[4]>0){
			$p2=log($tp[4])/log(2);
		}else{
			$p2=-1;
		}
		if($tp[7]>0){
			$p3=log($tp[7])/log(2);
		}else{
			$p3=-1;
		}						

		if($avg>0.3 || $max >1){
			if($tp[1]==0){
				if($tp[4]>1){
					$log2v1=log($tp[4])/log(2);
				}
			}else{
				if($tp[4]==0){
					if($tp[1]>1){
						$log2v1=-1* log($tp[1])/log(2);
					}
				}else{
					$log2v1=log($tp[4]/$tp[1])/log(2);
				}
			}
			
			if($tp[4]==0){
				if($tp[7]>1){
					$log3v2=log($tp[7])/log(2);
				}
			}else{
				if($tp[7]==0){
					if($tp[4]>1){
						$log3v2=-1* log($tp[4])/log(2);
					}
				}else{
					$log3v2=log($tp[7]/$tp[4])/log(2);
				}
			}
			if($log2v1>=1){
				$code1="<";
			}elsif($log2v1<=-1){
				$code1=">";			
			}
			if($log3v2>=1){
				$code2="<";			
			}elsif($log3v2<=-1){
				$code2=">";			
			}
			if($tp[1]>0 && $tp[4]>0 && $tp[7]>0){
				if(($tp[4]/$tp[1])>($tp[3]*2) && ($tp[4]/$tp[1])>($tp[6]*2)){
					$scode1=$code1;
				}
				if(($tp[7]/$tp[4])>($tp[6]*2) && ($tp[7]/$tp[4])>($tp[9]*2)){
					$scode2=$code2;
				}
			}					
		}	
		print OUT "$tp[0]\t$tp[1]\t$tp[4]\t$tp[7]\t$avg\t$log2v1\t$log3v2\t$code1$code2\t$scode1$scode2\t$p1\t$p2\t$p3\n";
	}

	print "done\n";
}

sub SeperateEachLineValue{
      
 	my @info;
 	my $infono=0;
 	my $n=0;
 	my @name;
 	
 	$lineIN=<IN>;
 	chomp($lineIN);
 	my @tp= &Getheadinfo($lineIN);
 	@name=@tp;
 	$n=@tp;
   	
	while(
	$lineIN=<IN>){	
		chomp($lineIN);
		@tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		for(my $i=0;$i<$tp;$i++){
			$info[$infono][$i]=$tp[$i];
		}
		
		#print "$info[$infono][0]\n";
		$infono++;
	}
	print "n=$n\tinfono=$infono\n";
	#print "aa=$info[0][0]\t$info[0][1]\n";
	
	for(my $i=1;$i<=$n;$i++){
		my $outfname=$name[$i].".txt";
		print "$outfname\n";
		open(OUTF,'>',"$outfname");		
		for(my $j=0;$j<$infono;$j++){
			print OUTF "$info[$j][0]\t$info[$j][$i]\n";
			
		}
		close(OUTF);	

	}
	print "done\n";
}




sub SeperateGroupPair{
      
	#input
	#ID	W	W	W	W	W	LN	LN	LN	LN	LN	HN	HN	HN	HN	HN
	#F0983	0.999394755	0.999798252	0.999899126	0.999764627	0.99971419	0.999899126	0.999739408	0.999831876	0.999731002	0.999495629	0.999529254	0.999831876	0.99936113	0.999764627	0.998117014

   	#outpput
   	#	G1	G1	G1	G1	G1	G2	G2	G2	G2	G2	
   	#W=LN=F0983	0.999394755	0.999798252	0.999899126	0.999764627	0.99971419	0.999899126	0.999739408	0.999831876	0.999731002	0.999495629
   	
   	my $n=$_[1];
   	if($n<=0){
   		$n=5;
   	}
  	
   	$lineIN=<IN>;
   	my @tp= &Getheadinfo($lineIN);
   	my $tp=@tp;
   	my $count=0;
   	my @group=();
   	for(my $i=1;$i<$tp;$i+=$n){
   		$group[$count]=$tp[$i];
   		$count++;
   	}
   	
 	my @info;
 	my $infono=0;
   	
	while($lineIN=<IN>){	
		chomp($lineIN);
		@tp= &Getheadinfo($lineIN);
		$info[$infono]=@tp;
		$infono++;
	}
	
	for(my $i=0;$i<$group;$i++){
	my $outfname=$name.".up";
	open(OUTF,'>',"$outfname");	
	
	close(OUTF);

	}
	print "done\n";
}




sub Gtf_txt_2_join_loc{
      
	#input
	#XLOC_000024     TCONS_00000031  Chr1_A_oryzae_RIB40     +       177787  178235  177787  178235  2       177787,178131,  178059,178235,
   	#outpput
   	#>ID|1|Chr:complement(join(2354845..2355351,2355357..2355806))
   	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my @st=split(/\,/,$tp[9]);
		my @ed=split(/\,/,$tp[10]);
		my $st=@st;
		my $name=">$tp[1]=$tp[0]|1|";
		my $loc="";
		if($st>1){
			$loc="join($st[0]..$ed[0]";
			for(my $i=1;$i<$st;$i++){
				$loc=$loc.",$st[$i]..$ed[$i]";
			}
			$loc=$loc.")";
		}else{
			$loc="$st[0]..$ed[0]";
		}
		if($tp[3] eq "+"){
			$name=$name."$tp[2]:$loc\n";
		}else{
			$name=$name."$tp[2]:complement($loc)\n";
		}
		print OUT "$name";
	}

	print "done\n";
}



sub DEseq2_input_from_HTseq{
	#input  compair list
	#CK_24	CK_48
	
	my $n=3;
	my $dirname="";
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		$dirname=$tp[0]."_vs_".$tp[1];
		system("mkdir $dirname");
		for(my $i=1;$i<=$n;$i++){
			my $tpname1="$tp[0]"."_R$i.reads.count.out";
			my $tpname2="$tp[1]"."_R$i.reads.count.out";
			system("cp $tpname1 $dirname/treated$i.txt");
			system("cp $tpname2 $dirname/untreated$i.txt");
		}
	}
	print "done\n";
}

sub DEseq_2_Goseq_input{
	#input  
	#GeneID			baseMean	log2FoldChange	lfcMLE		lfcSE		stat		pvalue		padj   
    #XLOC_005236	1411.047135	3.461955606		3.56448419	0.224711403	15.40623024	1.49E-53	1.03E-49
 	#output
 	#CK3_LA3_up	gene5862
 	
    my $qvalue=0.05;    
    my $name=substr($ARGV[0],0,index($ARGV[0],".",));  
    
   	my $up="";
   	my $down="";
   	my $uplist="";
	my $downlist="";
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN,"\"",",");
		my $tp=@tp;
		my $tpq=new Math::BigFloat $tp[7];   
		#print "$tpq qvalue=$qvalue\n";
		if($tpq<=$qvalue){
			if($tp[2] >= 1 ){
				$up="$up$name"."_up\t$tp[0]\n";	
				$uplist="$up$tp[0]\n";	
			}elsif($tp[2] <= -1){
				$down="$down$name"."_down\t$tp[0]\n";	
				$downlist="$downlist$tp[0]\n";
			}
		}		
	}
	print OUT "$up$down";
	#
	my $outfname=$name.".up";
	open(OUTF,'>',"$outfname");	
	print OUTF $uplist;
	close(OUTF);
	#
	my $outfname=$name.".down";
	open(OUTF,'>',"$outfname");	
	print OUTF $downlist;
	close(OUTF);
	
	print "done\n";
}

sub GetStringAfterString{
      
    my $symble=$_[1];
   	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my $out="NA";
		for(my $i=0;$i<$tp;$i++){
			if($tp[$i] eq $symble){
				print OUT "$tp[$i+1]\n";
				last;
			}
		}
		print $out "\n";
	}

	print "done\n";
}


sub ChangeSymbleToRef{
      
    my $symble=".";
    if($_[3]){
    	$symble=$_[3];
    }     	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($tp[$_[1]] eq $symble){
			$tp[$_[1]] = $tp[$_[2]];
		}
		print OUT "$tp[0]";
		for(my $i=1;$i<$tp;$i++){
			print OUT "\t$tp[$i]";
		}
		print OUT "\n";
	}

	print "done\n";
}


sub PvalueFoldchangeCheck{
		#ID	CT-W-T	CT-W-log2fold	CT-LN-T	CT-LN-log2fold	CT-HN-T	CT-HN-log2fold
		#OTU_368	1.44E-06	-1.622930351	0.027003967	-1.111031312	0.234609829	-0.641546029
		

    $pcut=0.05;
    #0=nochange,1-up,2=down
    #out put single score and Composite score
    
    #tittle
    $lineIN=<IN>;
		chomp($lineIN);
		my @name= &Getheadinfo($lineIN);
		my $name=@name;
		my $groupno=0;
		my @group;
		for(my $i=1;$i<$name;$i+=2){
			$group[$groupno]=$name[$i];
			$groupno++;
		}
		
		my @id;
		my $count=0;
		  	
		while($lineIN=<IN>){	
				chomp($lineIN);
				my @tp= &Getheadinfo($lineIN);
				my $tp=@tp;
				$id[$count][0]=$tp[0];
				my $tpg=1;
				for(my $i=1;$i<$tp;$i+=2){					
						if($tp[$i] eq "#DIV/0!"){							
							$id[$count][$tpg]=0;
						}elsif($tp[$i] <= $pcut){
								if($tp[$i+1] eq "#DIV/0!"){
									$id[$count][$tpg]=2;
									 print ".";
								}elsif($tp[$i+1] eq "#NUM!" ){
									$id[$count][$tpg]=1;
									print ".";
								}elsif($tp[$i+1]>=1){
									$id[$count][$tpg]=1;
									print ".";
								}elsif($tp[$i+1]<=-1){
									$id[$count][$tpg]=2;
									print ".";
								}else{
									$id[$count][$tpg]=0;
								}
						}else{
							$id[$count][$tpg]=0;
						}
						$tpg++;
				}
				$count++;
		}
		print "$groupno group and $count OTUs\n";
		
		
		#output single
		#head
		print OUT "ID";
		for(my $j=0;$j<$groupno;$j++){
			print OUT "\t$group[$j]";
		}		
		print OUT "\n";
		#score
		for(my $i=0;$i<$count;$i++){
			print OUT "$id[$i][0]";		
			for(my $j=0;$j<$groupno;$j++){
				print OUT "\t$id[$i][$j+1]";
			}
			print OUT "\n";
		}
		
		#output pair
		my $outfpair=$ARGV[0].'.pair';
		open(OUTP,'>',"$outfpair");	
		#head
		print OUTP "ID";
		for(my $i=0;$i<$groupno;$i++){
			for(my $j=$i+1;$j<$groupno;$j++){
					print OUTP "\t$group[$i]-vs-$group[$j]";
			}
		}		
		print OUTP "\n";
		#score
		for(my $k=0;$k<$count;$k++){
			print OUTP "$id[$k][0]";
			for(my $i=0;$i<$groupno;$i++){
				for(my $j=$i+1;$j<$groupno;$j++){
					print OUTP "\t$id[$k][$i+1]$id[$k][$j+1]";
				}
			}
			print OUTP "\n";
		}		
		
		close(OUTP);
		
		
	print "done\n";
}



sub PPI_Link_Clean{
	
	#gene6908	pp	gene12641
	#gene6969	pp	gene11559
      
    my %links;
        	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		
		my $add="$tp[0]\tpp\t$tp[2]";
		if($tp[0] gt $tp[2]){
			#print "do\n";
			$add="$tp[2]\tpp\t$tp[0]";
		}else{
		}
		if($links{$add}>0){
			$links{$add}++;
			print ".";
		}else{
			$links{$add} =1;
			#print "+";
		}

	}

	my @outs=sort(keys(%links));
	my $outs=@outs;
	for(my $i =0;$i<$outs;$i++){
		print OUT "$outs[$i]\n";
	}

	print "done\n";
}

sub AddNucleotide2fas{
      
    my $seq="";     	
	while($lineIN=<IN>){	
		chomp($lineIN);
		if($lineIN=~/>/){
			print OUT "$lineIN\n";
			if(length($seq)>0){
				print OUT "$seq$_[1]\n";
				$seq="";
			}
		}else{
			$seq=$seq.$lineIN;
		}	
	}
	print OUT "$seq$_[1]\n";
	print "done\n";
}



sub Hmp_2_NN{
#	rs#     alleles chrom   pos     strand  assembly#       center  protLNID        assayLNID       panelLNID       QCcode  C100 
#	S1_1742 C/T     1       1742    +       NA      NA      NA      NA      NA      NA      C       C       Y       C  
#	C/T	T/C	Y
#	G/C	C/G	S
#	T/G	G/T	K
#	A/C	C/A	M
#	A/G	G/A	R
#	A/T	T/A	W

	$lineIN=<IN>;
	print OUT "$lineIN";

	while($lineIN=<IN>){	
		chomp($lineIN);
		my $o=substr($lineIN,0,index($lineIN,"\t",));
		my $tp=substr($lineIN,index($lineIN,"\t",),length($lineIN));
		$tp =~ tr/YSKMRW/N/;
		print OUT "$o$tp\n";

		### double
		#my @tp= split(/\t/,$lineIN);
		#my $tp=@tp;
		#print OUT "$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]\t$tp[8]\t$tp[9]\t$tp[10]";
		#for(my $i=11;$i<$tp;$i++){
		#	if($tp[$i] eq "A"){$tp[$i] = "AA";}
		#	if($tp[$i] eq "T"){$tp[$i] = "TT";}
		#	if($tp[$i] eq "G"){$tp[$i] = "GG";}
		#	if($tp[$i] eq "C"){$tp[$i] = "CC";}
		#	if($tp[$i] eq "Y"){$tp[$i] = "NN";}
		#	if($tp[$i] eq "S"){$tp[$i] = "NN";}
		#	if($tp[$i] eq "K"){$tp[$i] = "NN";}
		#	if($tp[$i] eq "M"){$tp[$i] = "NN";}
		#	if($tp[$i] eq "R"){$tp[$i] = "NN";}
		#	if($tp[$i] eq "W"){$tp[$i] = "NN";}
		#	print OUT "\t$tp[$i]";
		#}
		#print OUT "\n";

	}
	print "done\n";
}

sub HTseq_to_DESeq{
	#remove head and tail
	 
	#35812 GFF lines processed.
	#XLOC_009440	0
	#__no_feature	7712272
	   	
	while($lineIN=<IN>){	
		chomp($lineIN);
		if(index($lineIN,"processed",)>=0){
			#print ".";
		}elsif(index($lineIN,"__",)>=0){
			#print "+";
		}else{
			print OUT "$lineIN\n";
		}
	}
	print "done\n";
}

sub Kobas_file_to_annolist{

	#////
	#Query:              	gene3
	#Gene:               	afv:AFLA_122970	
	#Entrez Gene ID:      	7920912
	#Pathway:            	Metabolic pathways	KEGG PATHWAY	afv01100
	#                    	Biosynthesis of secondary metabolites	KEGG PATHWAY	afv01110
	#                    	Starch and sucrose metabolism	KEGG PATHWAY	afv00500
	#                   	Cyanoamino acid metabolism	KEGG PATHWAY	afv00460
     
    my $Query="";
    my $Gene="";
    my $Entrez="";
    my $Pathway="";
    my $add=0;               	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		#print "$tp[0]\n";
		if($tp[0] eq "Query:"){
			$Query=$tp[1];
		}elsif($tp[0] eq "Gene:"){
			$Gene=$tp[1];
		}elsif($tp[0] eq "Entrez"){
			$Entrez=$tp[3];
		}elsif($tp[0] eq "Pathway:"){
			my $Pathway_term=$tp[1];
			for(my $i=2;$i<($tp-3);$i++){
				$Pathway_term="$Pathway_term\_$tp[$i]";
			}
			print OUT "$Query\t$Gene\t$Entrez\t$tp[$tp-1]\t$Pathway_term\n";		
			$add=1;
			
		}elsif($tp[0] eq "////"){
			$add=0;
			
		}else{
		
			if($add==1){
				my $Pathway_term=$tp[0];
				for(my $i=1;$i<($tp-3);$i++){
					$Pathway_term="$Pathway_term\_$tp[$i]";
				}
				print OUT "$Query\t$Gene\t$Entrez\t$tp[$tp-1]\t$Pathway_term\n";
			}
		}
	
	}
	print "done\n";
}


sub vcf_filter{
#	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  DSW11769.sorted.bam
#	chr01   1004    .       A       .       32.9906 .       DP=1;MQ0F=0;AF1=0;AC1=0;DP4=1,0,0,0;MQ=29;FQ=-29.9956   GT:PL   0/0:0
#	chr01   1005    .       A       .       32.9874 .       DP=1;MQ0F=0;AF1=0;AC1=0;DP4=1,0,0,0;MQ=29;FQ=-29.9987   GT:PL   0/0:0
#	chr01   1037    .       TAAA    .       22.5538 .       INDEL;IDV=1;IMF=0.2;DP=5;SGB=-0.379885;MQ0F=0;AF1=1;AC1=2;DP4=0,0,1,0;MQ=29;FQ=-37.6609 GT:PL   0/0:15
	
	my $depth_cut=8;
	my $p_cut=1;
	
	if($_[1] >0){
		$depth_cut=$_[1];
		$p_cut=$_[2];		
	}
		
	my $name=substr($ARGV[0],0,(index($ARGV[0],".",)));
	my @namearray=split(/_/,$name);
	my $namearray=@namearray;
	
	my $fname=substr($ARGV[0],0,(length($ARGV[0])-4));
	my $vcf_indel=$fname.'.Indel.vcf';
	open(OUTVI,'>',"$vcf_indel");
	my $vcf_SNP=$fname.'.clean.SNP.vcf';
	open(OUTVS,'>',"$vcf_SNP");
	my $out_SNP=$fname.'.clean.SNP.txt';
	open(OUTS,'>',"$out_SNP");
	my $out_v=$fname.'.clean.v.txt';
	open(OUTV,'>',"$out_v");	

	my $add=0;
	while($lineIN=<IN>){
		if(substr($lineIN,0,1) ne "#"){
			#print "$lineIN";
			if($lineIN=~/INDEL/ ){
				print OUTVI "$lineIN";
				print ".";
			}else{
				chomp($lineIN);
				my @tp= &Getheadinfo($lineIN,";");
				my $tp=@tp;
            	
				my $dp=substr($tp[7],3,length($tp[7]));
				my @dp4=();				
				for(my $i=7;$i<$tp;$i++){
					if($tp[$i]=~/DP4/){
						@dp4=split(/,/,substr(($tp[$i],4,length($tp[$i]))));
						last;
					}
				}
				#print "$dp=@dp4\n";
				#check
				
				if($dp >= $depth_cut && ( ($dp4[0]+$dp4[1])/$dp >= $p_cut || ($dp4[2]+$dp4[3])/$dp >= $p_cut)){
					my $v=0;
					my $output_snp="";
					print OUTVS "$lineIN\n";
					$output_snp= "$tp[0]\t$tp[1]\t$tp[3]\t";
					if($tp[4] eq "."){
						$output_snp=$output_snp."$tp[3]";
					}else{
						$output_snp=$output_snp."$tp[4]";
						$v=1;
					}
					for(my $i=0;$i<$namearray;$i++){
						$output_snp=$output_snp."\t$namearray[$i]";
					}
					
					print OUTS $output_snp."\n";
					if($v==1){
						print OUTV $output_snp."\n";
					}
		
				}else{
					#print "-";
				}
						
			}
		}
	}
	print "\ndown\n";
	close(OUTVI);
	close(OUTVS);
	close(OUTS);
	close(OUTV);
}

sub CPCgetCodingList{
	
	#SAT_000028=TCONS_00000075	354	noncoding	-1.13398
	#(revcomp)SAT_000028=TCONS_00000075	354	noncoding	-1.114
	
	my %coding=();

	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $score=2;
		if($tp[3]>0){
			$score=2;
		}else{
			$score=1;
		}		
		if($lineIN=~/(revcomp)/){
			my $name=substr($tp[0],9,length($tp[0]));
			if($coding{$name}<$score){
				$coding{$name}=$score;				
			}
		}else{
			$coding{$tp[0]}=$score;
		}		
	}
	#print
	my @outkeys=sort(keys(%coding));
	my $outkeys=@outkeys;
	for(my $i=0;$i<$outkeys;$i++){
		if($coding{$outkeys[$i]}==2){
			print OUT $outkeys[$i]."\n";
		}
	}

	print "\ndown\n";
}


sub DrawSaturationCurve{

	#my @dba;   	#put input data	
	$data_len=0;
	
	#read data
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($data_len==0){
			$data_len=@tp;
		}
		$dba[$dblevel]=[@tp];
		$dblevel++;	
	}
	
	#at all satuation
	#get list
	print "data length = $data_len; data row = $dblevel\n";
	for(my $i=1;$i<=$data_len;$i++){	
		$loop_len=$i;
		print "loop length = $loop_len\n";
		print OUT "Lv$i";
		@loop_arr;
		my $loop_st=0;
		my $loop_floor=0;
		&GroupLoop($loop_st,$loop_floor);
		print OUT "\n";
	}
	#count
	
	print "done\n";
}

sub GroupLoop{
	
	my $st=$_[0];
	my $floor=$_[1];
	#print "$st $floor\n";
	for(my $i=$st;$i<$data_len;$i++){
		$loop_arr[$floor]=$i;		
		if($floor==$loop_len-1){
			print "@loop_arr\n";
			CoutLineTotalState($loop_len,@loop_arr);
			
		}else{
			&GroupLoop($i + 1, $floor + 1);
		}
		
	}
}

sub CoutLineTotalState{
	my $r=0;
	
	for(my $i=0;$i<$dblevel;$i++){
		my $tpadd=0;
		for(my $j=0;$j<$_[0];$j++){
			$tpadd+=$dba[$i][$_[$j+1]];
			#print "=$dba[$i][$_[$j+1]].";
		}
		if($tpadd>0){
			$r++;
		}
	}
	print OUT "\t$r";
}


sub CheckStatebyValue{
	my $cut=0;
	if($_[1]){
		$cut=$_[1];
	}
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		print OUT "$lineIN";
		for(my $i=0;$i<$tp;$i++){
			if($tp[$i]>=$cut){
				print OUT "\t1";
			}else{
				print OUT "\t0";
			}
		}
		print OUT "\n";
	}
	print "done\n";
}



sub Removefilehead{
	my $rm=1;
	if($_[1]>1){
		$rm=$_[1];
	}
	my $count=0;
	while($lineIN=<IN>){	
		chomp($lineIN);	
		if($count<=$rm){
			$count++;
		}else{
			print OUT "$lineIN\n";
		}
	}
}


sub FindlocalOverlap{
	@loc1=(0,1,2);
	@loc2=(3,4,5);
	if($_[1]>=0){
		@loc1=($_[1],$_[2],$_[3]);
		@loc2=($_[4],$_[5],$_[6]);
	}	
	print "@loc1 vs @loc2\n";
	my $count=0;
	while($lineIN=<IN>){	
		my @tp= &Getheadinfo($lineIN);
		my $o=-1;
		if($tp[$loc1[0]] eq $tp[$loc2[0]]){
			
			$o=&FindOverlaplength($tp[$loc1[1]],$tp[$loc1[2]],$tp[$loc2[1]],$tp[$loc2[2]]);	
			#print "$o\n";	
		}
		print OUT "$lineIN\t$o\n";
	}

}

sub FindOverlaplength{ #s1 sorted
	#0-s1, 1-e1, 2-s2, 3-e2
	#return overlap length
	#print "@_\n";
	my $r=0;
	if($_[1]>=$_[2] && $_[1]<=$_[3]){
		if($_[0]>=$_[2]){
			$r=$_[1]-$_[2]+1;
		}else{
			$r=$_[1]-$_[0]+1;   #include
		}
		return $r;
	}
	if($_[2]>=$_[0] && $_[3]<=$_[1]){
		$r=$_[3]-$_[2]+1;
		return $r;
	}
	if($_[2]<=$_[0] && $_[3]>=$_[1]){
		$r=$_[3]-$_[0]+1;
		return $r;
	}	
	if($_[0]>=$_[2] && $_[0]<=$_[3]){
		if($_[1]>=$_[3]){
			$r=$_[1]-$_[0]+1;
		}else{
			$r=$_[1]-$_[0]+1;  #inclued
		}
		return $r;
	}
	return $r;
}

sub GetSequencefromEstdata{
	#||
	#dbEST Id:       19864234
	#EST name:       BpHi 045
	#GenBank Acc:    CF568773
	#GenBank gi:     34994856	
	#SEQUENCE
	#                TACTAACTGGAGTGAGCAGGGGCTGCATTACTTCTTTCTTCTGTCAGCGCAGCAATAGCC

	my $name="";
	my $seq="";
	my $count=0;
	my $canadd=0;
	while($lineIN=<IN>){	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($tp[0] eq "||"){
			$canadd=0;	
			$name="";
			$seq="";				
		}
		elsif($tp[0] eq "dbEST" && $canadd==0){
			$name = $tp[2]; 
			$canadd=1;				
		}		
		elsif($tp[0] eq "GenBank" && $canadd==1){
			$name = $name."|$tp[1]$tp[2]"; 		
		}
		elsif($tp[0] eq "SEQUENCE" && $canadd==1){
			$canadd=2;	
		}	
		elsif($canadd==2){
			if(length($tp[0])>0){
				$seq=$seq.$tp[0];
			}else{
				print OUT ">$name\n$seq\n";
				$count++;
				$name="";
				$seq="";
				$canadd=3;
			}	
		}	
	}
	print "$count sequence readed\n";
}



sub SlideWindow_continue_step{
	
#input  0-chrom 1-st 2-ed 3-value
#out avg_value
	my $n = 10;
	if($_[1]>0 ){
		$n = $_[1];
	}

	my @data;  
	my $count=0;
	my $chrom="";
	while($lineIN=<IN>){	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		
		if($chrom ne $tp[0]){
			if($count >0){	
				my $add=1;	
				for(my $i=0;$i<$count;$i++){
					my $o=0;
					
					my $ied=$i+$n-1;
					if($ied>=($count-1)){
						$ied=($count-1);
						$add=0;
					}
					#print "$count $i $ied $data[$ied][2]\n";
					for(my $j=$i;$j<=$ied;$j++){
						$o=$o+$data[$j][3];
					}
					$o=$o/($ied-$i+1);
					print OUT "$data[$i][0]\t$data[$i][1]\t$data[$ied][2]\t$o\n";
					if($add==0){
						last;
					}
				}
				#new
				@data=();
				$count=0;	
				
			}				
			$chrom = $tp[0];	
		}
		$data[$count]=[@tp];
		$count++;
	}
	#lst one
				my $add=1;	
				for(my $i=0;$i<$count;$i++){
					my $o=0;
					
					my $ied=$i+$n-1;
					if($ied>=($count-1)){
						$ied=($count-1);
						$add=0;
					}
					#print "$count $i $ied $data[$ied][2]\n";
					for(my $j=$i;$j<=$ied;$j++){
						$o=$o+$data[$j][3];
					}
					$o=$o/($ied-$i+1);
					print OUT "$data[$i][0]\t$data[$i][1]\t$data[$ied][2]\t$o\n";
					if($add==0){
						last;
					}	
				}

}



sub SlideWindowCount_sorttype{
	
#input  0-type 1-loc(sorted) 2-value 
	my $winsize = 5;
	if($_[1]>0 ){
		$winsize = $_[1];
	}
	my $stepsize = 1;
	if($_[2]>0 ){
		$stepsize = $_[2];
	}	
	print "winsize=$winsize stepsize=$stepsize\n";
	my @data;
	my $count=0;
	my $name="";
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($tp[0] ne $name){
			$name=$tp[0];
			$count=0;
			@data=();
		}
		
		#add
		$data[$count]=[@tp];
		
		#sum
		my $sum=0;
		my $denominator = $winsize;
		for(my $i=0;$i<$winsize;$i++){
			if($count-$i>=0){
				$sum=$sum+$data[$count-$i][2];
			}else{
				$denominator=$i;
				last;
			}
		}
		#avg
		my $avg=$sum/$denominator;
		
##### cut
		my $sum_10=0;
		my $denominator_10 = 10;
		for(my $i=0;$i<10;$i++){
			if($count-$i>=0){
				$sum_10=$sum_10+$data[$count-$i][2];
			}else{
				$denominator_10=$i;
				last;
			}
		}
		my $avg_10=$sum_10/$denominator_10;		
		
		my $sum_20=0;
		my $denominator_20 = 20;
		for(my $i=0;$i<20;$i++){
			if($count-$i>=0){
				$sum_20=$sum_20+$data[$count-$i][2];
			}else{
				$denominator_20=$i;
				last;
			}
		}
		my $avg_20=$sum_20/$denominator_20;				
		
		print OUT $lineIN."\t$avg\t$avg_10\t$avg_20\n";
		
		$count++;
	}

}


sub SlideWindowCount{
	
#input  0-loc(sorted) 2-value 
	my $winsize = 100;
	if($_[1]>0 ){
		$winsize = $_[1];
	}
	my $stepsize = 10;
	if($_[2]>0 ){
		$stepsize = $_[2];
	}	
	print "winsize=$winsize stepsize=$stepsize\n";
	my @data;
	my $count=0;
	while($lineIN=<IN>){	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		$data[$count]=$tp[1];	
		$count++;
	}
		
	for(my $i=0;$i<$count;$i+=$stepsize){
		my $o=0;
		for(my $j=0;$j<$winsize;$j++){
			$o=$o+$data[$i+$j];
		}
		print OUT "$i\t$o\n";
	}
}



sub MergedGtfGetTranscriptsAndGene{
	#Chr1	Cufflinks	exon	2895	3255	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "LOC_Os01g01010"; oId "CUFF.10.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	
	my $count=0;
	my $info="";	
	my $name="NA";
	my $st="";
	my $ed=0;
	my $tail="";
	
	my $g_name="NA";
	my $g_st="";
	my $g_ed=0;
	my $g_tail="";	
	my $g_info="";
		
	my $outt=$ARGV[0].".transcripts";
	open(OUTS,'>',"$outt");	
	my $outg=$ARGV[0].".gene";
	open(OUTG,'>',"$outg");		
			
	while($lineIN=<IN>){	
		chomp($lineIN);		
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($tp[11] eq $name){
			#add
			$ed=$tp[4];
			$info=$info."$lineIN\n";
		}else{			
			if($name ne "NA"){
				#output
				print OUT "$st$ed$tail$info";
				print OUTS "$st$ed$tail";
			}
			#add new
			$name = $tp[11];
			$st="$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t";
			$ed=$tp[4];
			$tail="\t$tp[5]\t$tp[6]\t$tp[7]\t$tp[8]\t$tp[9]\t$tp[10]\t$tp[11]";	
			for(my $i=12; $i<$tp;$i+=2){
				if($tp[$i] ne "exon_number" ){
					$tail=$tail."\t$tp[$i]\t$tp[$i+1]";
				}
			}		
			$tail=$tail."\n";
			
			#new
			$info="$lineIN\n";
		}
		
		if($tp[9] eq $g_name){
			#add
			$g_ed=$tp[4];
		}else{			
			if($g_name ne "NA"){
				#output
				print OUTG "$g_st$g_ed$g_tail";
			}
			#add new
			$g_name = $tp[9];
			$g_st="$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t";
			$g_ed=$tp[4];
			$g_tail="\t$tp[5]\t$tp[6]\t$tp[7]\t$tp[8]\t$tp[9]\t$tp[10]\t$tp[11]";	
			for(my $i=12; $i<$tp;$i+=2){
				if($tp[$i] ne "exon_number" ){
					$g_tail=$g_tail."\t$tp[$i]\t$tp[$i+1]";
				}
			}		
			$g_tail=$g_tail."\n";
		}		
		
	}
	close(OUTS);
	close(OUTG);
}


sub Gtf_2_join_loc{
	
	#gtf
	#Chr1	Cufflinks	exon	2895	3255	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr1	Cufflinks	exon	3354	3616	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "2"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr1	Cufflinks	exon	4357	4455	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "3"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr1	Cufflinks	exon	5457	5763	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "4"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr11	Cufflinks	exon	2087097	2087678	.	-	.	gene_id "XLOC_013958"; transcript_id "TCONS_00022146"; exon_number "1"; gene_name "LOC_Os11g04910"; oId "CUFF.5957.4"; nearest_ref "LOC_Os11g04910.1"; class_code "j"; tss_id "TSS16278";
	#Chr11	Cufflinks	exon	2091028	2091107	.	-	.	gene_id "XLOC_013958"; transcript_id "TCONS_00022146"; exon_number "2"; gene_name "LOC_Os11g04910"; oId "CUFF.5957.4"; nearest_ref "LOC_Os11g04910.1"; class_code "j"; tss_id "TSS16278";
	#Chr11	Cufflinks	exon	2091178	2091321	.	-	.	gene_id "XLOC_013958"; transcript_id "TCONS_00022146"; exon_number "3"; gene_name "LOC_Os11g04910"; oId "CUFF.5957.4"; nearest_ref "LOC_Os11g04910.1"; class_code "j"; tss_id "TSS16278";

   	#outpput
   	#>ID|1|Chr:complement(join(2354845..2355351,2355357..2355806))
	
	my $count=0;
	my $gene_id="NA";
	my $gene_dir="NA";
	my $chrom=0;
	my $transcript_no=0;	
	my $transcript_id="NA";
	my $exon_no=0;
	my @transcript_st=0;
	my @transcript_ed=0;	
	
	open(OUTG,'>',"$outfg");		
	
	while($lineIN=<IN>){	
		my @tp= &Getheadinfo($lineIN, "=",";",":","\"");	
		my $tp=@tp;
		if($tp[11] ne $transcript_id){
			#print transcript 
			if($transcript_id ne "NA"){
				#output
				my $out= "$transcript_id=$gene_id|$transcript_no|$chrom:";
				my $loc=0;
				if($exon_no>1){
					$loc="join($transcript_st[0]..$transcript_ed[0]";
					for(my $i=2;$i<$exon_no;$i++){
						$loc=$loc.",$transcript_st[$i]..$transcript_ed[$i]";
					}
					$loc=$loc.")";
				}else{
					$loc="$transcript_st[0]..$transcript_ed[0]";
				}
				#direct
				if($transcript_dir eq "+"){
					$out=$out.$loc."\n";
				}else{
					$out=$out."complement($loc)\n";
				}
				
				#new transcript
				$transcript_no++;
				$exon_no=0;
			}
			#new gene
			if($tp[9] ne $gene_id){
				#newgene
				$chrom=$tp[0];
				$gene_dir=$tp[6];
				$transcript_no=1;
			}
		
			#add information
			$transcript_st[$exon_no]=$tp[3];		
			$transcript_ed[$exon_no]=$tp[4];	
			$exon_no++;				
		}
	}
	#last one
	#output
	my $out= "$transcript_id=$gene_id|$transcript_no|$chrom:";
	my $loc=0;
	if($exon_no>1){
		$loc="join($transcript_st[0]..$transcript_ed[0]";
		for(my $i=2;$i<$exon_no;$i++){
			$loc=$loc.",$transcript_st[$i]..$transcript_ed[$i]";
		}
		$loc=$loc.")";
	}else{
		$loc="$transcript_st[0]..$transcript_ed[0]";
	}
	#direct
	if($transcript_dir eq "+"){
		$out=$out.$loc."\n";
	}else{
		$out=$out."complement($loc)\n";
	}

	print "done\n";


}


sub Gtf_2_gene{
	
	#gtf
	#Chr1	Cufflinks	exon	2895	3255	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr1	Cufflinks	exon	3354	3616	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "2"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr1	Cufflinks	exon	4357	4455	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "3"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr1	Cufflinks	exon	5457	5763	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "4"; gene_name "LOC_Os01g01010"; oId "CUFF.18.2"; nearest_ref "LOC_Os01g01010.1"; class_code "j"; tss_id "TSS1";
	#Chr11	Cufflinks	exon	2087097	2087678	.	-	.	gene_id "XLOC_013958"; transcript_id "TCONS_00022146"; exon_number "1"; gene_name "LOC_Os11g04910"; oId "CUFF.5957.4"; nearest_ref "LOC_Os11g04910.1"; class_code "j"; tss_id "TSS16278";
	#Chr11	Cufflinks	exon	2091028	2091107	.	-	.	gene_id "XLOC_013958"; transcript_id "TCONS_00022146"; exon_number "2"; gene_name "LOC_Os11g04910"; oId "CUFF.5957.4"; nearest_ref "LOC_Os11g04910.1"; class_code "j"; tss_id "TSS16278";
	#Chr11	Cufflinks	exon	2091178	2091321	.	-	.	gene_id "XLOC_013958"; transcript_id "TCONS_00022146"; exon_number "3"; gene_name "LOC_Os11g04910"; oId "CUFF.5957.4"; nearest_ref "LOC_Os11g04910.1"; class_code "j"; tss_id "TSS16278";

	
	my $count=0;
	my $gene_id="NA";
	my $transcript_id="NA";
	my $chrom=0;
	my $gene_st=0;
	my $gene_ed=0;
	my $transcript_st=0;
	my $transcript_ed=0;	
	my $gene_dir="NA";
	my $transcript_dir="NA";
	my $exon_no=0;
	my $transcript_no=0;
	my $gene_name="NA";
	my $nearest_ref="NA";
	
	my $outfg=$ARGV[0].'.gene';
	
	open(OUTG,'>',"$outfg");		
	
	while($lineIN=<IN>){	
		my @tp= &Getheadinfo($lineIN, "=",";",":","\"");	
		my $tp=@tp;
		if($tp[11] eq $transcript_id){
			$transcript_ed=$tp[4];		
			$gene_ed=$tp[4]>$gene_ed?$tp[4]:$gene_ed;
			$exon_no++;
	
		}else{
			if($transcript_id ne "NA"){
				#out put
				print OUT "$transcript_id\t$gene_id\t$chrom\t$transcript_st\t$transcript_ed\t$transcript_dir\t$exon_no\t$gene_name\t$nearest_ref\n";
			}
			#new gene
			if($tp[9] ne $gene_id){
				#out
				if($gene_id ne "NA"){
					print OUTG "$gene_id\t$chrom\t$gene_st\t$gene_ed\t$gene_dir\t$transcript_no\t$gene_name\n";
				}
				#newgene
				$chrom=$tp[0];
				$gene_dir=$tp[6];
				$transcript_no=1;
				$gene_st=$tp[3];
				$gene_ed=$tp[4];
				$gene_name="NA";
				for(my $i=0;$i<$tp;$i++){
					if($tp[$i] eq "gene_id"){
						$gene_id=$tp[$i+1];
					}
					if($tp[$i] eq "gene_name"){
						$gene_name=$tp[$i+1];
					}				
				}					
			}else{
				$transcript_no++;
			}	
			#new transcript	
			$transcript_dir=$tp[6];	
			$nearest_ref="NA";
			for(my $i=0;$i<$tp;$i++){
				if($tp[$i] eq "transcript_id"){
					$transcript_id=$tp[$i+1];
				}
				if($tp[$i] eq "nearest_ref"){
					$nearest_ref=$tp[$i+1];
				}				
			}			
			$exon_no=1;
			$transcript_st=$tp[3];
			$transcript_ed=$tp[4];	
			$gene_ed=$tp[4]>$gene_ed?$tp[4]:$gene_ed;				
		}				
	}
	print OUT "$transcript_id\t$gene_id\t$chrom\t$transcript_st\t$transcript_ed\t$transcript_dir\t$exon_no\t$gene_name\t$nearest_ref\n";
	print OUTG "$gene_id\t$chrom\t$gene_st\t$gene_ed\t$gene_dir\t$transcript_no\t$gene_name\n";
	
	close(OUTG);
}

sub Gff3_2_gtf{
	
	#gff3
	#Sc0002422	KIB	gene	26725	30691	.	-	.	ID=CSA000001;
	#Sc0002422	KIB	mRNA	26725	30691	.	-	.	ID=CSA000001.1;Parent=CSA000001;
	#Sc0002422	KIB	exon	30403	30691	.	-	.	ID=CSA000001.1:exon_1;Parent=CSA000001.1;
	#Sc0002422	KIB	CDS	30403	30691	.	-	0	ID=CSA000001.1:cds_1;Parent=CSA000001.1;

	#gtf
	#Chr1	MSU_osa1r7	exon	2903	3268	.	+	.	transcript_id "LOC_Os01g01010.1"; gene_id "LOC_Os01g01010"; gene_name "LOC_Os01g01010";
	#Chr1	MSU_osa1r7	CDS	3449	3616	.	+	0	transcript_id "LOC_Os01g01010.1"; gene_id "LOC_Os01g01010"; gene_name "LOC_Os01g01010";

	
	my $count=0;
	my $gene="NA";
	
	while($lineIN=<IN>){	
		my @tp= &Getheadinfo($lineIN, "=",";",":");
		if($tp[2] eq "gene"){
			$gene=$tp[9];	
			#print "$tp[9]\n";		
		}elsif($tp[2] eq "mRNA"){
			$mRNA=$tp[9];	
		}elsif($tp[2] eq "exon" || $tp[2] eq "CDS"){
			print OUT "$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]\ttranscript_id \"$mRNA\"; gene_id \"$gene\"; gene_name \"$gene\";\n";
		}				
	}

}


sub Sam_Addloc_RemoveSinglehit{
	#input must sorted
	
	my $count=0;
	
	
	while($lineIN=<IN>){	
		#shoot.13120091	89	Chr1	149928	50	4M1I220N95M	*	0	0	TTAAGAG......
		#MIDN
		my $mn=0;
		my $in=0;
		my $nn=0;
		my $dn=0;
		
		
	}
	
	#output
	my $outt=$ARGV[0].".single";
	open(OUTS,'>',"$outt");	
	for(my $i=0;$i<$count;$i++){
		
	}
	close(OUTS);
}


sub GetNRfirstname{
	while($lineIN=<IN>){	
		if($lineIN=~/>/){
			my $l=index($lineIN,"[",0)-1;
			my $s=substr($lineIN,0,$l);
			print OUT $s."\n";
		}
	}	

}



sub trinotate_annotation_venn{

	#trinotate_annotation_report
	#gene_id				0
	#transcript_id			1	#combin
	#sprot_Top_BLASTX_hit	2	#findbest	by E-value
	#TrEMBL_Top_BLASTX_hit	3	#findbest	by E-value
	#RNAMMER				4	#.
	#prot_id				5	#.
	#prot_coords			6	#.
	#sprot_Top_BLASTP_hit	7	#findbest	by E-value
	#TrEMBL_Top_BLASTP_hit	8	#findbest	by E-value
	#Pfam					9	#pfam combin
	#SignalP				10	#.
	#TmHMM					11	#.
	#eggnog					12	#COG/KOG combin
	#gene_ontology_blast	13	#GO combin
	#gene_ontology_pfam		14	#GO_pfam combin
	#transcript				15	#.
	#peptide				16	#.
	#if have				17	#nr
		
	print OUT "ID\t\tUniPort\tTrEMBL\tCOG_KOG\tGO\tNR\n";
	while($lineIN=<IN>){
		chomp($lineIN);		
		my @tp= split(/\t/,$lineIN);
		my $tp=@tp;
		my $UniPort=0;
		if($tp[2] ne "." || $tp[7] ne "." ){
			$UniPort=1;
		}
		my $TrEMBL=0;
		if($tp[3] ne "." || $tp[8] ne "." ){
			$TrEMBL=1;
		}	
		my $Pfam=0;
		if($tp[9] ne "."){
			$Pfam=1;
		}	
		my $COG_KOG=0;
		if($tp[12] ne "."){
			$COG_KOG=1;
		}	
		my $GO=0;
		if($tp[13] ne "." || $tp[14] ne "." ){
			$GO=1;
		}	
		my $nr=0;
		if($tp[17] ne "." ){
			$nr=1;
		}							
		print OUT "$tp[0]\t\t$UniPort\t$TrEMBL\t$COG_KOG\t$GO\t$nr\n";		
	}	

}


sub PrankforCombinedFile{
	#export PATH=$PATH:/home/zqj/soft/prank/bin
	#prank -d=1.txt -t=Gene_6_align_used_8k.fas.ph -o=2_output -f=nexus -once -showxml -showall +F
	system("export PATH=$PATH:/home/zqj/soft/prank/bin");
	
	#out put combin events file and best.acn.fas file 	

	#my $outd1=$ARGV[0].'.best.anc.fas';
	#open(ANS,'>',"$outd1");	

	my $outd2=$ARGV[0].'.best.indel';
	open(OID,'>',"$outd2");	
	print OID "name\tbranch\tst\ted\tlen\ttype\tterminal\n";
		
	my $name="NA";
	my $tpfile="";
	
	my $treefile="Gene_6_align_used_8k.fas.ph";
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN, "_");
		my $tp=@tp;		
		if($lineIN=~/>/){
			$tp[0]=substr($tp[0],1,length($tp[0]));
			#print "$tp[0]\n";
			if($name ne $tp[0]){
				if($name ne "NA"){
					#out put to tp file and run
					close(TPF);
					my $tprun = "prank -d=".$tpfile." -t=".$treefile." -o=".$name." -f=nexus -once -showxml -showall +F";
					#print "$tprun\n";
					system($tprun); 
					my $tpanc = &read_anc_to_file($name);	
					#print ANS $tpanc."\n";	
					my $tpind = &read_events_to_file($name);
					print OID $tpind;										

				}
				#add new
				$name = $tp[0];
				$tpfile = "$name.fas";
				open(TPF,'>',"$tpfile");			
				print OUT "$name\n";
			}
			print TPF ">$tp[1]\n";
		}elsif(length($lineIN)>0){
			print TPF "$lineIN\n";
		}
	}
	close(TPF);
	my $tprun = "prank -d=".$tpfile." -t=".$treefile." -o=".$name." -f=nexus -once -showxml -showall +F";
	print "$tprun\n"; 	
	system($tprun);
	my $tpanc = &read_anc_to_file($name);	
	#print ANS $tpanc."\n";	
	my $tpind = &read_events_to_file($name);	
	print OID $tpind;		

	#close(ANS);	
	close(OID);	
}

sub read_anc_to_file{
	my $tpfilename=$_[0].".best.anc.fas";
	print "$tpfilename\n";
	my $r="";
	open(TAN,'<',"$tpfilename");
	my $lineTAN=''; 
	#read in and get jap etc
	my $add=1;
	while($lineTAN=<TAN>){
		#print "$lineTAN\n";
		chomp($lineTAN);			
		if($lineTAN=~/>/){
			$add=1;
			if($lineTAN=~/\#/){
				$add=0;
			}else{
				$r=$r.">".$_[0]."_".substr($lineTAN,1,length($lineTAN))."\n";
			}
		}else{
			if($add==1){
				$r=$r.$lineTAN."\n";
			}
		}
	}	
	return $r;
	close(TAN);	
}

sub read_events_to_file{
	my $tpname=$_[0];
	my $tpfilename=$_[0].".best.events";
	#return:  name branch st ed len type terminal
	
	print "$tpfilename\n";
	my $r="";
	
	open(TAN,'<',"$tpfilename");
	my $lineTAN=''; 
	#read in 
	#branch jap
	# 1..29 deletion (terminal)
 	# 1390..1390 insertion
 	# 1497..1497 deletion	
	
	my $head="";
	while($lineTAN=<TAN>){
		#print "$lineTAN\n";
		chomp($lineTAN);
		my @tpe= &Getheadinfo($lineTAN, ".");			
		if($tpe[0] eq "branch"){
			$head="$tpname\t$tpe[1]";
		}elsif($tpe[2] eq "insertion" || $tpe[2] eq "deletion"){
			my $tplen=$tpe[1]-$tpe[0]+1;
			$r=$r.$head."\t$tpe[0]\t$tpe[1]\t$tplen\t$tpe[2]";
			if($tpe[3] eq "(terminal)"){
				$r=$r."\t1\n";
			}else{
				$r=$r."\t0\n";
			}
		}
	}	
	return $r;
	close(TAN);	
}



sub CounttypeAllinfoNo{

#combin the prnak events.detail get by PrankEventsCount
#e.g.
#branch_deletion|r7	21	1	10	10	1	1	1	2	114	8	3	3	32	2	1	11	23	3	2	2	406	7	4	4	7	3	1	1	1	2	1	1	1	2	19	2	19	3	1	9	9	3	9	3	64	4	12	6	9	6	3	13	4	3	4	6	11	7	7	1	4	4	1	41	10	8	1	3	8	10	3	10	11	1	12	6	13	3	4	6	9	3	5	2	1	1	6	52	11	247	3	4	3	5	50	13	3	1	9	1	28	1	8	13	34	138	4	3	10	1	1	5	2048	7	6	5	10	1	5	4	1	2	2	1	2	1	3	11	3	1	1	9	12	6	10	3	4	3	15	6	1	1	8	7	7	1	4	22	3	3	8	1	4	3	9	6	16	5	4	1	46	3	255	7	12	1	7	4	6	11	30	6	13	5	12	3	3	36	4	4	447	248	3	1	1	3	1	1	8	4	8	6	4	108	5	11	6	1	5	64	3	2	4	6	4	150	10	32	2	20	6	21	1	6	4	25	4	1	11	3	6	3	12	9	6	6	8	7	22	5	3	8	1	2	20	1	8	66	178	1	1	2	1	1	1	46	89	5	2	4	1	5	12	6	1	10	14	1	1	1	4	3	9	4	3	2	2	1	2	1	3	2	1	6	2	114	22	1	1	1	5	5	40	11	1	2	1	1	5	1	1	1	5	1	6	259	26	1	1	2	21	1	3	1	1	3	5	4	6	5	9	1	4	1	1	7	4	2	7	1	6	5	6	30	6	2	6	1	3	27	1	3	8	5	13	6	49	13	6	11	4	1	14	5	3	12	5	1	2	1	1	4	2	11	7	161	1	7	1	30	14	2	5	6	15	1	22	1	2	254	8	6	1	7	1	1	1	1	1	16	5	3	6	7	13	3	6	1	32	1	3	1	2	6	2	7	22	1	3576	27	445	6	19	3	4	13	111	8	6	10	7	4	3	6	2	1	3	16	4	15	6	19	19	2	13	1	3	1	1	3	2	8	4	10	4	402	7	3	8	43	2	3	14	454	2	100	5	1	5	2	1	1	1	2	2	2	3	2	15	6	3	1	3	2	2	1	8	1	1	144	16	1	4	29	1038	2	2	1	30	6	4	5	1	2	1	1	1	1	4	2	3	18	9	9	21	41	37	39	7	15	3	3	41	3	1	3	1	3	12	68	2	3	1	2	22	40	3	5	1	5	2	4	2	1	8	5	6	3	9	2	3	1	5	363	4	1	5	1	3	3	3	6	3	4	4	1	1	19	7	7	6	1	12	1	1	3	1	3	5	7	2	1	2	152	1	6	1	16	6	4	1	1	2	3	1	12	3	3	10	1	1	2	4	21	1	1	1	2	29	10	1	2	3	3	2	2	1	5	4	10	32	16	4	3	57	3	1235	2	9	3	1	25	10	1	4	1	12	2	8	10	30	3	18	3	12	3	4	7	1	15	3	102	4	6	8	1	180	3	13	4	3	1127	1	2	1503	42	9	1	2	1	3	27	1	1	1
#branch_insertion|r5	59	175	21	21	1	58	1	1	119	12	257	92	5	51	3	1	786	18	163	1	36	14	205	1	10	25	26	64	9	9	24	28	250	22	10	13
#include mutation

	my $count=1;
	my $typeno=1;
	my $max=0;
	my @all;  #row-subtype col-type
	$all[0][0]="Type";
	while($lineIN=<IN>){
		print ".";
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		#check class
		my $c=-1;
		for(my $i=1;$i<$typeno;$i++){
			if($all[0][$i] eq $tp[0]){
				$c=$i;
				last;
			}
		}
		if($c==-1){
			$c=$typeno;
			for(my $k=1;$k<$count;$k++){
				$all[$k][$c]=0;					
			}			
			$all[0][$c]=$tp[0];
			$typeno++;
			print "+";
		}
		#check sub type	
		#print "$tp ";	
		for(my $i=1;$i<$tp;$i++){
			if($tp[$i]>$max){
				$max=$tp[$i];
			}
			my $r=-1;
			for(my $j=1;$j<$count;$j++){
				if($all[$j][0] eq $tp[$i]){
					$r=$j;
					$all[$r][$c]++;
					last;
				}				
			}
			if($r==-1){
				$r=$count;				
				for(my $k=1;$k<$typeno;$k++){
					$all[$r][$k]=0;					
				}
				$all[$r][0]=$tp[$i];
				$all[$r][$c]++;
				#print "$all[$r][$c]";
				$count++;
				#print ".";
			}			
		}		
	}
	print "\ncol=$count\trow=$typeno\tmax=$max\n";
	#output
	print OUT $all[0][0];
	for(my $k=1;$k<$typeno;$k++){
		print OUT "\t$all[0][$k]";
	}
	print OUT "\n";

	for(my $i=1;$i<=$max;$i++){
		for(my $j=1;$j<=$count;$j++){
			if($all[$j][0] eq $i){
				#print $i;
				print OUT $i;
				for(my $k=1;$k<$typeno;$k++){
					print OUT "\t$all[$j][$k]";
				}
				print OUT "\n";
				last;
			}
		}
	}
	print "\ndone\n";
}


sub RNAhybridOutputCheck{

	#input file like this
	#target: XLOC_008141=1=-|-
	#length: 99
	#miRNA : osa-miR156a
	#length: 20
	#
	#mfe: -20.0 kcal/mol
	#p-value: 0.069938
	#
	#position  1
	#target 5'   C         CAGCU   AGA     A 3'
	#             UUCACUUUC     UUC   UGUCA    
	#             GAGUGAGAG     AAG   ACAGU    
	#miRNA  3' CAC                           5'

	#chang to code
	#0-no match
	#1-match
	#2-target only
	#3-miRNA only
	#output code: target miRNA code miRNA_length target_st target_length match_no mis_match_no targte_bulge_no  miRNA_bulge_no
	print OUT "target\tmiRNA\tcode\t$revcode\tmiRNA_length\ttarget_length\ttarget_st\tmatch_no\tmis_match_no\ttargte_bulge_no\tmiRNA_bulge_no\tOK\n";
	#check rule
	#endogenous microRNA target mimic (eTM).
	#(1) bulges were only permitted at the 5' end ninth to 12th positions of miRNA sequence; 
	#(2) the bulge in eTMs should be composed of only three nucleotides; 
	#(3) perfect nucleotide pairing was required at the 5' end second to eighth positions of miRNA sequence; and
	#(4) except for the central bulge, the total mismatches and G/U pairs within eTM and miRNA pairing regions should be no more than three.

	
	my $count=0;
	my $target="";
	my $miRNA="";
	my $code="";
	my $len=0;
	my $st=0;
	my $tlen=0;
	my $loc=0;
	my $mis=0;
	my $tin=0;
	my $tde=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		#if($tp[0] eq "target:"){
		if($lineIN=~/target:/){
			$count++;
			$target=$tp[1];
			$miRNA="";
			$code="";			#last site may add to 1
			$revcode="";
			$len=0;
			$st=0;
			$tlen=-1;			#last site may add to 1
			$loc=0;	
			$mis=0;
			$tin=0;
			$tde=0;
			#print "$target\n";		
						
		}elsif($lineIN=~/miRNA :/){
			$miRNA=$tp[2];
		}elsif($lineIN=~/length:/){
			$len=$tp[1];
		}elsif($lineIN=~/position/){
			$st=$tp[1];
		}elsif($lineIN=~/target 5/){
			my $tg=substr($lineIN,10,length($lineIN));	
			
			my $mi="";			
			$lineIN=<IN>;
			$lineIN=<IN>;
			$lineIN=<IN>;
			chomp($lineIN);
			if($lineIN=~/miRNA  3/){
				$mi=substr($lineIN,10,length($lineIN));			
			}
			#print "$tg\n$mi\n";		
			#check
			#0-no match 1-match  2-target only  3-miRNA only
			for(my $i=0;$i<length($tg);$i++){
				my $ct=substr($tg,$i,1);	
				my $cm=substr($mi,$i,1);	
				#print "$ct ";
										
				#output code: target miRNA code $revcode miRNA_length target_length target_st match_no mis_match_no targte_bulge_no  miRNA_bulge_no OK?
				if($ct eq "3"){
					print ".";
					#output 	
					#code last 1 must remove
					$code=substr($code,0,(length($code)-1));
					$revcode=substr($revcode,1,(length($revcode)));
					
					my $ok = &checkRNAhybridbyrule($revcode);
					print "$ok\n";
					print OUT "$target\t$miRNA\t$code\t$revcode\t$len\t$tlen\t$st\t$mis\t$tin\t$tde\t$ok\n";
					last;
				}				
				elsif($ct eq " " && $cm eq " "){
					$code=$code."1";
					$revcode="1".$revcode;
					#print "1";
					$tlen++;
				}
				elsif($ct ne " " && $cm ne " "){
					$code=$code."0";
					$revcode="0".$revcode;
					#print "0";
					$tlen++;
					$mis++;
				}
				elsif($ct ne " " && $cm eq " "){
					$code=$code."2";
					$revcode="2".$revcode;
					#print "2";
					$tin++;
					$tlen++;
				}						
				elsif($ct eq " " && $cm ne " "){
					$code=$code."3";
					$revcode="3".$revcode;
					#print "3";
					$tde++;
				}

			}			
		}
	}
	print "\n";
}

sub checkRNAhybridbyrule{
	#check rule
	#endogenous microRNA target mimic (eTM).
	#(1) bulges were only permitted at the 5' end ninth to 12th positions of miRNA sequence; 
	#(2) the bulge in eTMs should be composed of only three nucleotides; 
	#(3) perfect nucleotide pairing was required at the 5' end second to eighth positions of miRNA sequence; and
	#(4) except for the central bulge, the total mismatches and G/U pairs within eTM and miRNA pairing regions should be no more than three.
	
	#target 5'  A      C GC      G  3'
    #    	     GCUCAC C   UCUGU     
    #	    	 CGAGUG G   AGACA     
	#miRNA  3' CA      A AGA     GU 5'

	#revcode 30111113001011111103
	#
	
	#2-8 must 1
	#total 2 must ==3
	my $n2max=3;
	#2 must in 9 to 12
	#3 only in start and end not in middle
	
	print "$_[0]\t";
	my $n2=0;
	my $add=1;
	my $len=length($_[0]);
	for(my $i=0;$i<$len;$i++){
		my $c=substr($_[0],$i,1);
		#print "$c";
		#total 2 must ==3
		if($c eq "2"){			
			$n2++;
			if($n2>$n2max){
				return "no2";
			}
		}
		#2 must in 9 to 12
		if($i<9 || $i >12){
			if($c eq "2"){
				#print "$c";
				return "no1";
			}
		}
		#3 only in start and end not in middle
		if($i>2 && $i <($len-2) && $c eq "3"){
			return "no3";
		}
		#2-8 must 1
		#if($i>2 && $i<8 && $c eq "0"){
		#	return;
		#}		
	}
	return "ok";
		
}
sub PranknexandEventsMutateCount{
	#input is list
	#get length from nex
	#-------------
	##NEXUS
	#begin data;
	#dimensions ntax=6 nchar=717482;
	
	my %tstv = ("A -> G",1, "G -> A",1, "C -> T",1, "T -> C",1, "A -> C",2, "A -> T",2, "C -> A",2, "T -> A",2, "G -> C",3, "G -> T",2, "C -> G",3, "T -> G",2);	
	my $ts=0;			
	my $tv=0;
	my $gc=0;	
	my $mut=0;
	my $len=0;
	my $inno=0;
	my $inlen=0;
	my $deno=0;
	my $delen=0;
	print OUT "ID\tLength\tMutation\tTs\tTv\tInNo\tInLen\tDelNo\tDelLen\n";	
	#readin list
	while($lineIN=<IN>){
		chomp($lineIN);
		
		#nex
		$len=0;
		my $nex="$lineIN.best.nex";
		open(TPIN,'<',"$nex");
		my $lineTPIN;
		$lineTPIN=<TPIN>;
		$lineTPIN=<TPIN>;
		$lineTPIN=<TPIN>;
		chomp($lineTPIN);	
		my @tp1= &Getheadinfo($lineTPIN);
		$len=substr($tp1[2],6,(length($tp1[2])-7));
		close(TPIN);
		
		#events
		# 79 G -> C
		$ts=0;			
		$tv=0;	
		$gc=0;
		$mut=0;	
		$inno=0;
		$inlen=0;
		$deno=0;
		$delen=0;
		my $event="$lineIN.best.events";
		open(TPIN,'<',"$event");
		while($lineTPIN=<TPIN>){
			chomp($lineTPIN);	
			my @tp= &Getheadinfo($lineTPIN);
			if($lineTPIN=~/->/){
				
				my $t=$tstv{"$tp[1] -> $tp[3]"};
				if($t == 1){
					$ts++;
				}elsif($t == 2){
					$tv++;
				}elsif($t == 3){
					$tv++;
					$gc++;
				}else{
					$mut++;
				}
				#print "$t";
			}elsif($tp[1] eq "deletion"){
				my @tploc=&Getheadinfo($tp[0],".");
				#print "$tp[0]=$tploc[1]-$tploc[0]\n";
				my $tplen=$tploc[1]-$tploc[0]+1;
				$deno++;
				$delen+=$tplen;
			}elsif($tp[1] eq "insertion"){
				my @tploc=&Getheadinfo($tp[0],".");
				my $tplen=$tploc[1]-$tploc[0]+1;
				$inno++;
				$inlen+=$tplen;	
			}
		}		
		close(TPIN);
		
			
		$mut=$mut+$ts+$tv;
		print OUT "$lineIN\t$len\t$mut\t$ts\t$tv\t$gc\t$inno\t$inlen\t$deno\t$delen\n";	
	}	

	print "done\n";
}



sub PrankEventsCount{
	#input is a list
	
	#seperate
	#branch seq_name
	# 79 G -> C
 	# 153 A -> G
	# 94..105 deletion
	# 366012..366012 insertion
	#inner brance
	#branch #2#
	
	my $count=-1;
	my @name; 
	my @inno;
	my @deno;
	my @inlen;
	my @delen;
	my @indetail;
	my @dedetail;
	my @ts;			# 
	my @tv;
	
	my %tstv = ("G -> C",1, "C -> G",1, "T -> A",1, "A -> T",1, "A -> G",2, "A -> C",2, "T -> G",2, "T -> C",2, "G -> A",2, "G -> T",2, "C -> A",2, "C -> T",2);	

	my $outn=$ARGV[0].'.detail';	
	open(OUTN,'>',"$outn");		
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($tp[0] eq "branch"){
			$count++;
			$ts[$count]=0;
			$tv[$count]=0;
			$name[$count]=$tp[1];
			$inno[$count]=0;
			$deno[$count]=0;
			$inlen[$count]=0;
			$delen[$count]=0;
			$indetail[$count]="branch_insertion|". $tp[1];
			$dedetail[$count]="branch_deletion|". $tp[1];		
		}elsif($tp[1] eq "deletion"){
			my @tploc=&Getheadinfo($tp[0],".");
			#print "$tp[0]=$tploc[1]-$tploc[0]\n";
			my $tplen=$tploc[1]-$tploc[0]+1;
			$deno[$count]++;
			$delen[$count]+=$tplen;
			$dedetail[$count]=$dedetail[$count]."\t$tplen";
		}elsif($tp[1] eq "insertion"){
			my @tploc=&Getheadinfo($tp[0],".");
			my $tplen=$tploc[1]-$tploc[0]+1;
			$inno[$count]++;
			$inlen[$count]+=$tplen;
			$indetail[$count]=$indetail[$count]."\t$tplen";			
		}elsif($lineIN=~/->/){
			my $t=$tstv{"$tp[0] -> $tp[2]"};
			if($t == 1){
				$ts++;
			}elsif($t == 2){
				$tv++;
			}else{
				print "?";
			}
		}
	}
	$count++;
	print "$count\n";
	print OUT "branch\tinsertion_no\tavg_len\ttotal_len\tdeletion_no\tavg_len\ttotal_len\tMutate_No\tTs\tTv\n";
	for(my $i=0;$i<$count;$i++){
		my $inavg=0;
		if($inno[$i]>0){
			$inavg=$inlen[$i]/$inno[$i];
		}
		my $deavg=0;
		if($deno[$i]>0){
			$deavg=$delen[$i]/$deno[$i];
		}
		my $mutno=$ts[$i]+$tv[$i];
		print OUT "$name[$i]\t$inno[$i]\t$inavg\t$inlen[$i]\t$deno[$i]\t$deavg\t$delen[$i]\t$mutant[$i]$ts[$i]\t$tv[$i]\n";
	}
	print OUT "\n\n###detail####\n";
	for(my $i=0;$i<$count;$i++){
		print OUTN "$indetail[$i]\n$dedetail[$i]\n";
	}	
	
	close(OUTN);
}

sub ReverseFilebyRow{

	my $count=0;
	my @all;
	while($lineIN=<IN>){
		chomp($lineIN);
		$all[$count]=$lineIN;
		$count++;
	}
	print "$count\n";
	for(my $i=$count;$i>0;$i--){
		print OUT "$all[$i-1]\n";
	}

}

sub ReversedOrderForFasta{

	my @name;
	my @seq;
	my 
	$count=-1;
	
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			$count++;
			$name[$count]=$lineIN;			
		}else{
			$seq[$count]=$seq[$count].$lineIN;	
		}
	}
	#out put
	for(my $i=$count;$i>=0;$i--){
		print OUT "$name[$i]\n$seq[$i]\n";
	}
}

sub SeperateTreesReadfromfilelist{
	my $name='NA';
	my $count=0;
	while($lineIN=<IN>){
		$count++;	
		chomp($lineIN);	
		print OUT "$lineIN\t";
		print "$lineIN\n";
		open(INTP,'<',"$lineIN");
		my $lineTP;
		while($lineTP=<INTP>){
			chomp($lineTP);
			my $add=0;
			for(my $i=0; $i<length($lineTP);$i++){
				my $a=substr($lineTP,$i,1);
				if($a eq ":"){
					$add=0;
				}
				if($add==1){
					print OUT $a;
				}
				if($a eq "(" || $a eq ")" || $a eq "," || $a eq ";"){
					print OUT $a;
				}
				if($a eq "_"){
					$add=1;
				}
			}
		}
		close(INTP);
		print OUT "\n";	
    }
    print "$count trees readed!\n";
}



sub Seperategff3locToregionloc{
	#infile 
	#13101.t01015	1	2215	13101.t01015|13101.m01237	-	1	1	1	1	626	627	1926	1927	2215	
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($tp[4] eq "+"){
			for (my $i=0; $i<$tp[5];$i++){
				my $a=1+$i;
				print OUT "$tp[0]\t$tp[8+$i*2]\t$tp[9+$i*2]\t$tp[3].+.UTR5_".$a."\n";					
			}
			for (my $i=0; $i<$tp[6];$i++){
				my $a=1+$i;
				print OUT "$tp[0]\t$tp[8+2*$tp[5]+$i*2]\t$tp[9+2*$tp[5]+$i*2]\t$tp[3].+.exon_".$a."\n";					
			}
			for (my $i=0; $i<$tp[7];$i++){
				my $a=1+$i;
				print OUT "$tp[0]\t$tp[8+2*($tp[5]+$tp[6])+$i*2]\t$tp[9+2*($tp[5]+$tp[6])+$i*2]\t$tp[3].+.UTR3_".$a."\n";					
			}
		}else{
			for (my $i=0; $i<$tp[7];$i++){
				my $a=$tp[7]-$i;
				print OUT "$tp[0]\t$tp[8+$i*2]\t$tp[9+$i*2]\t$tp[3].-.UTR3_".$a."\n";					
			}
			for (my $i=0; $i<$tp[6];$i++){
				my $a=$tp[6]-$i;
				print OUT "$tp[0]\t$tp[8+2*$tp[7]+$i*2]\t$tp[9+2*$tp[7]+$i*2]\t$tp[3].-.exon_".$a."\n";					
			}
			for (my $i=0; $i<$tp[5];$i++){
				my $a=$tp[5]-$i;
				print OUT "$tp[0]\t$tp[8+2*($tp[7]+$tp[6])+$i*2]\t$tp[9+2*($tp[7]+$tp[6])+$i*2]\t$tp[3].-.UTR5_".$a."\n";					
			}		
		}
	}
}

sub getMaskedRegionLocalfromFasta{

	#####chr start length
	
	my $count=0;
	my $name="NA";
	my $st=0;
	my $ed=0;
	my $im=0;   #in NNNN im=1
	my $loc=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp= &Getheadinfo($lineIN);
			print "+";
			if($name ne "NA"){
				if($im == 1){
					print OUT "$name\t$st\t$ed\t$count\n";
					$count++;
					
				}
			}
			$im=0;
			$name=$tp[0];
			$st=0;
			$ed=0;
			$loc=0;
		}else{		
			my @tp=split(//,$lineIN);  
			my $tp=@tp;			
			for(my $i=0;$i<$tp;$i++){
				$loc++;
				if($tp[$i] eq "N" || $tp[$i] eq "n" ){
					if($im == 0){
						$st=$loc;					
						$im=1;
					}
					$ed=$loc;
				}else{
					if($im == 1){
						print OUT "$name\t$st\t$ed\t$count\n";
						#print ".";
						$im=0;
						$st=0;
						$ed=0;	
						$count++;					
					}
				}
			}
		}
	}
	if($im == 1){
		$ed=$loc;
		print OUT "$name\t$st\t$ed\t$count\n";
		$count++;
	}
}

sub RemoveTotalIncludeData{

	#####chr start length 	
	#read in data sort by chr and score big>small £¨maf ¡µtable file etc£©
	
	my $name="NA";
	my $count=0;
	@arr;  #save data  0-stat 1-st 2-ed 3-$lineIN
	my $outd=$ARGV[0].'.delete';
	open(OUTD,'>',"$outd");
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($tp[$_[1]] ne $name){
			if($name ne "NA"){
				#output
				my $info=&checklocoverlap(@arr);
				print "$name total $count line delete $info\n";
			}
			$name=$tp[$_[1]];
			@arr=();
			$count=0;
		}		
		$arr[$count][0]=0;
		$arr[$count][1]=$tp[$_[2]];
		$arr[$count][2]=$tp[$_[2]]+$tp[$_[3]];
		$arr[$count][3]=$lineIN;
		$count++;
	}
	my $info=&checklocoverlap(@arr);
	print "$name total $count line delete $info\n";
	close(OUTD);
}


sub checklocoverlap{
	my $a=@_;
	#check
	my $c=0;
	for(my $i=0;$i<$a;$i++){
		if($_[$i][0] == 0){
			#print "$_[$i][1] $_[$i][2]\n";
			for(my $j=0;$j<$a;$j++){
				if($_[$i][1]>=$_[$j][1] && $_[$i][2] <= $_[$j][2] && $i != $j){
					$_[$i][0]=1;
					last;
				}
				elsif($_[$j][1]>=$_[$i][1] && $_[$j][2] <= $_[$i][2] && $i != $j){
					$_[$j][0]=1;
				}else{
				}

			}
		}
	}
	#output
	for(my $i=0;$i<$a;$i++){
		if($_[$i][0] == 0){
			print OUT "$_[$i][3]\n";
		}else{
			print OUTD "$_[$i][3]\n";
			$c++;
		}
	}
	return $c;
}


sub Methy_gene_Local_count{

	#input 
	#chr	loc	meth un_meth	gene_name	A	1999
	
	#ABCD ·Ö±ðÊÇ up£¨-5k-0£© st(st-half) ed(half-end) down(1-5k)
	
	my %data;  #0-meth; 1-un_meth
	my $win = 10;
	if($_[1]>0){
		$win = $_[1];
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);
		#get loc
		my @tp = split(/\t/,$lineIN);
		my $n=int $tp[6]/$win;
		my $loc=$tp[5];
		for(my $i=5; $i>length($n);$i--){
			$loc=$loc."0";
		}
		$loc=$loc.$n;
		#add
		if($data{$loc}){
			$data{$loc}[0]=$data{$loc}[0]+$tp[2];
			$data{$loc}[1]=$data{$loc}[1]+$tp[3];
		}else{
			$data{$loc}[0]=$tp[2];			
			$data{$loc}[1]=$tp[3];	
		}
	}
	
	#sort
	my @lkeys=keys(%data);
	my @lkeys=sort(@lkeys);
	my $lkeys=@lkeys;
	for(my $i=0;$i<$lkeys;$i++){
		$mp=$data{$lkeys[$i]}[0]/($data{$lkeys[$i]}[0] + $data{$lkeys[$i]}[1]);
		print OUT "$lkeys[$i]\t$data{$lkeys[$i]}[0]\t$mp\n";
	}
	
}



sub MethylationTypeIdentityINGenome{

	#input genome
	my $name="NA";
	my $seq="";
	my $loc=0;
	my $dir=0;   #0- 1+
	my @count;

	my $outn=$ARGV[0].".detail";
	open(OUTN,'>',"$outn");	

	$count[0]=$count[1]=$count[2]=$count[3]=0;
	#print "Chrom\tCpG\tCHG\tCHH\tTotalC\n";
	print OUT "Chrom\tCpG\tCHG\tCHH\tTotalC\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			if($name ne "NA"){
				my @r=&CheckMethylationType($seq,$name);
				#print "$name\t$r[0]\t$r[1]\t$r[2]\t$r[3]\n";
				print OUT "$name\t$r[0]\t$r[1]\t$r[2]\t$r[3]\n";
				$count[0]=$count[0]+$r[0];
				$count[1]=$count[1]+$r[1];
				$count[2]=$count[2]+$r[2];
				$count[3]=$count[3]+$r[3];
			}
			my @tp= &Getheadinfo($lineIN,">");
			$name=$tp[0];
			$seq="";
		}else{
			$seq=$seq.$lineIN;					
		}
	}
	my @r=&CheckMethylationType($seq,$name);
	$count[0]=$count[0]+$r[0];
	$count[1]=$count[1]+$r[1];
	$count[2]=$count[2]+$r[2];
	$count[3]=$count[3]+$r[3];	
	#print "$name\t$r[0]\t$r[1]\t$r[2]\t$r[3]\n";	
	print "Total\t$count[0]\t$count[1]\t$count[2]\t$count[3]\n";
	print OUT "$name\t$r[0]\t$r[1]\t$r[2]\t$r[3]\n";	
	print OUT "Total\t$count[0]\t$count[1]\t$count[2]\t$count[3]\n";	
		
	close(OUTN);
}

sub CheckMethylationType{
	#CHH +H -h
	#CHG +X -x
	#CpG +Z	-z
	my @tp=split(//,uc($_[0]));
	my $tp=@tp;
	my @c;   #0-z 1-x 2-h
	$c[0]=$c[1]=$c[2]=$c[3]=0;
	for(my $i=0;$i<$tp;$i++){
		my $t="N";
		if($tp[$i] eq "C" || $tp[$i] eq "c"){
			 $t="H";
			if($i<$tp-1){
				if($tp[$i+1] eq "G" || $tp[$i+1] eq "g"){
					 $t="Z";
					 $c[0]++;
				}else{
					if($i<$tp-2){
						if($tp[$i+2] eq "G" || $tp[$i+2] eq "g"){
							$t="X";
							$c[1]++;
						}
					}
				}								
			}
		}elsif($tp[$i] eq "G" || $tp[$i] eq "g"){
			$t="h";
			if($i>1){
				if($tp[$i-1] eq "C" || $tp[$i-1] eq "c"){
					 $t="z";
					 $c[0]++;
				}else{
					if($i>2){
						if($tp[$i-2] eq "C" || $tp[$i-2] eq "c"){
							$t="x";
							$c[1]++;
						}
					}
				}								
			}			
		}
		if($t ne "N" ){			
			my $loc=$i+1;
			#print "$_[1]\t$loc\t$t\n";
			print OUTN "$_[1]"."\t$loc\t$t\n";
			$c[3]++;
			
		}
		$c[2]=$c[3]-$c[1]-$c[0];
	}
	return @c;
}


sub FindPairdata{
	#5359    6087    730     1       729     730     96.99   772     94.56   OR_CBa0074D14.f
		my %pair=();
		my %info=();
		my $count=0;
		while($lineIN=<IN>){
			chomp($lineIN);
			my @tp = split(/\t/,$lineIN); 
			my @pairname = split(/\./,$tp[9]); 
			#print "$tp[9] $pairname[0] $pairname[1]\n";
			#
			if($pair{$pairname[0]}){
				if($pair{$pairname[0]} ne $pairname[1]){
					if($pair{$pairname[0]} ne "rl"){
						$info{$pairname[0]}=$info{$pairname[0]}."\t$lineIN";
						$pair{$pairname[0]}="rl";
					}else{
						print "3";
					}
				}else{	
				 print ".";
				}
			}else{
				$pair{$pairname[0]}=$pairname[1];
				$info{$pairname[0]}="$lineIN";
			}
		}	
		#OR_CBa0065I08   3247407 3248167 1       766     761     766     97.13   766     100.00  OR_CBa0065I08.f 29258423        29259257        1       835     835     835     97.84   835     100.00  OR_CBa0065I08.r
		my $outno="$ARGV[0].nosymble";
		open(OUTN,'>',"$outno");
		my @hkeys=keys(%info);
		my $hkeys=@hkeys;
		#
		
		for(my $i=0;$i<$hkeys;$i++){
			my @tp = split(/\t/,$info{$hkeys[$i]}); 
			my $len=$tp[10]-$tp[0];
			if($len<0){
				$len=$len*-1;
			}
			my $dir1="+";
			if(($tp[3]-$tp[2])<0){
				$dir1="-";
			}
			my $dir2="+";
			if(($tp[13]-$tp[12])<0){
				$dir2="-";
			}			
			print OUT "$hkeys[$i]\t$len\t$dir1\t$dir2\t$info{$hkeys[$i]}\n";
			
		}
		close(OUTN);
}



sub GetMafAlignedQualifiedRegion{
## name start alnSize strand seqSize alignment
#a score=40
#s Chr1           24734807 82 + 43270923 TTTTTCTGTGCGGTTGCACTTATGGGCCGCCTGGAAACAAACTAGTCTGGTGTCCAAAAAAATATTTAGTGTAGTAGTGCTA
#s Scaffold10_716      162 82 -      244 TTTTTCTGTGCGGTTCAACTCAAGGATCGCATGGAAATAAAATGGGGTGGTGTTCAGGAAAATCAATTATGTAGTAGTGATA


#step 1: mutation may change the type of CpG CHH CHG£¬so the up and down 2 bp must removed for mutation/indel
#step 2: use slide window remove region < 50

	my $region=50;
	if($_[1]>0){
		$region=$_[1];
	}
	#my $score=1000;
	my $gap=1;
	my $seq1="";
	my $seq2="";
	my $count=0;
	my $start="";
	
	my $outn=$ARGV[0].".detail";
	open(OUTN,'>',"$outn");	
	
	print      "###region length >= $region\n";
	print OUT  "###region length >= $region\n";
	print OUTN "###region length >= $region\n";
	print OUT  "###No\tscore\tname\tstart\tend\talnSize\tstrand\tseqSize\tname\tstart\talnSize\tstrand\tseqSize\tRegionCount\tRegionCount\tsimilarity\n";
	print OUTN "###No\tscore\tname\tstart\tend\talnSize\tstrand\tseqSize\tname\tstart\talnSize\tstrand\tseqSize\tRegionCount\tRegionCount\tsimilarity\tRegionDetail\n";
	my $add=0;
	my $count=0;
	

		
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN,"=");
		if($tp[0] eq "a"){			
			$add++;
			$count++;
			print OUT "$count\t$tp[2]";
			print OUTN "$count\t$tp[2]";
			#if()
		}
		elsif($tp[0] eq "s"){
			if($add==1){
				$seq1=$tp[6];
				$start=$tp[2];
			}
			if($add==2){
				$seq2=$tp[6];
			}
			for(my $i=1;$i<=5;$i++){
				print OUT "\t$tp[$i]";
				print OUTN "\t$tp[$i]";
			}			
			$add++;
		}
		if($add==3){
			#output
			&SeqAlignWindowSearch($seq1,$seq2,$region,$start);
			#print ".";
			$add=0;
			#$count++;
		}
	}
	print "\nTotal $count aligns\n";
	close(OUTN);
		
}
sub SeqAlignWindowSearch{
	#seq
	my @s1=split(//,uc($_[0]));
	my @s2=split(//,uc($_[1]));
	my @st;  
	my @ed;
	my $rc=0;  
	#para
	my $cut=$_[2];
	my $start=$_[3];
	#print "$cut\n";
	#
	my $ts=0;  #total score
	my $r_len=0;
	my $len=length($_[0]);
	#print "$len\n";
	
	#step 1, change seq to region
	
	my $gapst=0;
	my $same=1;
	$st[$rc]=0;
	#print "$ts";
	for(my $i=1;$i<$len;$i++){
		my $tpsc=1;
		if($s1[$i] eq $s2[$i]){
			#score	
			if($same==0){
				$st[$rc]=$i;
				$same=1;	
				
			}
			$gapst=0;	
		}else{
			#step 2, remove region and up/down 2 bp
			if($same==1){
				$ed[$rc]=$i-1;
				#print "$rc $st[$rc] $ed[$rc]\n";
				$rc++;
				
				
			}
			$same=0;
			$tpsc=0;			
			if($s1[$i] eq "-" || $s2[$i] eq "-"){
				if($gapst==0){
					$tpsc=-1;								
				}
				$gapst=1;
			}
		}	
		
		$ts=$ts+$tpsc;
	}	
	$ed[$rc]=$len-1;

	#region count
	my $okc=0;
	my $r="";
	for(my $i=0;$i<$rc;$i++){
		my $tplen=$ed[$i]-$st[$i]+1;
		if($tplen>=$cut){
			$okc++;
			$r_len+=$tplen;
			my $tpst=$st[$i]+$start;
			my $tped=$ed[$i]+$start;
			$r=$r."$tpst..$tped,";
		}
	}
	$ts=$ts*100/$len;
	print OUT "\t$rc\t$okc\t$r_len\t$ts\n";
	print OUTN "\t$rc\t$okc\t$r_len\t$ts\t$r\n";
	#print " $len $rc $okc $r_len $ts $r\n";
}

sub GetMafAlignedTable{

##       batch   0
## name start alnSize strand seqSize alignment

	print OUT "score\tname\tstart\talnSize\tstrand\tseqSize\tname\tstart\alnSize\tstrand\tseqSize\n";
	my $add=0;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN,"=");
		if($tp[0] eq "a"){			
			$add++;
			$count++;
			print OUT "$count\t$tp[2]";
		}
		if($tp[0] eq "s"){
			for(my $i=1;$i<=5;$i++){
				print OUT "\t$tp[$i]";
			}
			$add++;
		}
		if($add==3){
			print OUT "\n";
			$add=0;
			
		}
	}
	print "Total $count aligns\n";
		
}

sub GetUnigeneFromTrinityFasta{
	#>TR5|c0_g1_i3 len=431
	my $name="";
	my $trno="";
	my $seq="";
	my $len=0;
	my $count=0;
	my $add=0;
	while($lineIN=<IN>){
		#chomp($lineIN);
		if($lineIN=~/>/){
			my @tp = split(/ /,$lineIN);
			my @tpn = split(/_/,$tp[0]);
			my $tpname=$tpn[0]."_".$tpn[1];
			my @tpl=split(/=/,$tp[1]);
			my $tplen = $tpl[1];
			#print "$tpl[1]\n";
			if($tpname ne $name){
				if($len > 0){
					print OUT $name."_".$trno."\t$len\n$seq";
				}
				$len=$tplen;
				$name=$tpname;
				$trno=$tpn[2];
				$seq="";
				$add=1;
				$count++;
			}else{
				if($tplen > $len){
					$seq="";
					$add=1;
					$len=$tplen;
					$trno=$tpn[2];
				}else{
					$add=0;
				}
			}
			
		}else{
			if($add ==1){
				$seq=$seq.$lineIN;
			}
		}
	}
	print OUT $name."_".$trno."\t$len\n$seq";
	print "Total $count genes\n";
		
}

sub BismarkExtractorOutCount{

	#HWUSI-EAS611_0006:3:1:1058:15806#0/1 - 6 91793279 z
	#HWUSI-EAS611_0006:3:1:1058:17564#0/1 + 8 122855484 Z

	
	#print OUT "Chrom\tlocal\t+/-\tz\tZ\tx\tX\th\tH\tTotal\n";
	#input
	#0-reads_D 1-direct 2-chr 3-loc 5-M zZ xX hH
	#output
	#0=No 1-chr 2-st 5-all 6-m 7-u 

	# %db  #keys chr-loc value-i in data
	my @data;  #0-chrom 1-loc 2-dir 3-z 4-Z 5-x 6-X 7-h 8-H 9-count
	
	$count=1;
		
	while($lineIN=<IN>){
		chomp($lineIN);
		open(INT,'<',"$lineIN");
		#print "$lineIN\n";
		my $lineINT=''; 
		while($lineINT=<INT>){
			chomp($lineINT);
			my @tp = split(/\t/,$lineINT); 	
			my $tpkey=$tp[2]."|".$tp[3];
			#print "$tp[1]\n";
			if($db{$tpkey}>0){
				#print "$tpkey $db{$tpkey} $data[$db{$tpkey}][2] ne $tp[1]\n";
				if($data[$db{$tpkey}][2] ne $tp[1]){
					$data[$db{$tpkey}][2]="2";
					#print "$tpkey\n";
				}
				$data[$db{$tpkey}][9]++;	
				if($tp[4] eq "z")    { $data[$count][3]++;         }
				elsif($tp[4] eq "Z") { $data[$db{$tpkey}][4]++;   }
				elsif($tp[4] eq "x") { $data[$db{$tpkey}][5]++;   }
				elsif($tp[4] eq "X") { $data[$db{$tpkey}][6]++;   }
				elsif($tp[4] eq "h") { $data[$db{$tpkey}][7]++;   }
				elsif($tp[4] eq "H") { $data[$db{$tpkey}][8]++;	 }
			}else{
				#print "+";
				$data[$count][0]=$tp[2];
				$data[$count][1]=$tp[3];
				$data[$count][2]=$tp[1];
				$data[$count][3]=$data[$count][4]=$data[$count][5]=$data[$count][6]=$data[$count][7]=$data[$count][8]=0;
				$data[$count][9]=1;
				if($tp[4] eq "z")    { $data[$count][3]=1;  }
				elsif($tp[4] eq "Z") { $data[$count][4]=1;  }
				elsif($tp[4] eq "x") { $data[$count][5]=1;  }
				elsif($tp[4] eq "X") { $data[$count][6]=1;  }
				elsif($tp[4] eq "h") { $data[$count][7]=1;  }
				elsif($tp[4] eq "H") { $data[$count][8]=1;  }
				$db{$tpkey}=$count;			
				$count++;
			}
		}

	}
	
	print "Total $count C local Checked\nSort...\n";
	my @mykeys=keys(%db);
	@mykeys = sort(@mykeys);
	#print "@mykeys\n";
	my $mykeys=@mykeys;
	print "write...\n";
	for(my $i=1;$i<$mykeys;$i++){
		#print "$db{$mykeys[$i]}\n";
		for(my $j=0;$j<10;$j++){					
			print OUT "$data[$db{$mykeys[$i]}][$j]\t";
		}
		print OUT "\n";
	}
	
	print "done\n";
}

sub BismarkExtractorOutputSlideCount{

	#HWUSI-EAS611_0006:3:1:1058:15806#0/1 - 6 91793279 z
	#HWUSI-EAS611_0006:3:1:1058:17564#0/1 + 8 122855484 Z
	
	$window = 10000;
	$step   =  1000;
	if($_[1] > 0 && $_[2] > 0){
		$window = $_[1];
		$step   = $_[2];
	}
	if($step > $window){
		print "step must small than window\n";
		return;
	}
	$fold=$window%$step;
	if( $fold != 0){
		print "window must be the fold of step\n";
		return;
	}
	
	print OUT "window size = $window; step size = $step\n";
	#input is sorted by String[2] and int[3] 
	#input
	#0-reads_D 1-direct 2-chr 3-loc 5-M zZ xX hH
	#output
	#0=No 1-chr 2-st 5-all 6-m 7-u 
	my $outn=$lineIN.".w=$window.s=$step.txt";
	open(OUTN,'<',"$outn");	
	
	#chrom
	$chrom="NA";		
	$count=0;
	my @step=();  #0-all 1-m 2-u
		
	my $loc=1;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = split(/ /,$lineIN); 	
		#new chrom
		if($chrom ne $tp[2]){  
			if($chrom ne "NA"){
				&slidewindowsoutput(@step);
			}
			$chrom = $tp[2];
			$step=();
			$count=0;
			$step[$count][0]=$step[$count][1]=$step[$count][2]=0;				
		}
		
		while($loc<=$tp[3]){
			$loc++;
			if($loc==$tp[3]){
				$step[$count][0]++;		
				if($tp[4] eq z || $tp[4] eq x || $tp[4] eq h){  #u			
					$step[$count][2]++;
				}else{                                          #m			
					$step[$count][1]++;
				}
			}
			#check step
			if($loc%$step == 0){
				$count++;
				$step[$count][0]=$step[$count][1]=$step[$count][2]=0;
			}
		}		
	}
	&slidewindowsoutput(@step);
}

sub slidewindowsoutput{
	#@$_ 0-all 1-m 2-u
	#$chrom  $count  $window  $step
	my @c=0;
	$c[0]=$c[1]=$c[2]=0;
	for(my $i=0;$i<$count;$i++){
		my @o;
		$o[0]=$o[1]=$o[2]=0;
		for(my $j=0;$j<$fold;$j++){
			$o[0]+=$step[$i+$j][0];
			$o[1]+=$step[$i+$j][1];
			$o[2]+=$step[$i+$j][2];
			$c[0]+=$step[$i+$j][0];
			$c[1]+=$step[$i+$j][1];
			$c[2]+=$step[$i+$j][2];			
		}
		my $st=$i*$step;
		my $ed=($i+$j+1)*$step-1;
		print OUTN "$i\t$chrom\t$st\t$ed\t$o[0]\t$o[1]\t$o[2]\n"
	}
	print OUT "$chrom\t$c[0]\t$c[1]\t$c[2]\n";
	return;
}


sub ConnectSortedLocal{
	
	my $name="NA";  
	my $id = 0;
	my $st = 1;
	my $ed  =2;
	my $conn=0;
	if($_[1]>=0){
		$conn=$_[1];
	}
	my $data="";
	my $tpst=0;
	my $tped=0;
	my $len=0;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		#my @tp = split(/\t/,$lineIN);	
		if($name ne $tp[$id]){
			$count++;
			if($name ne "NA"){
				my $len=$tped-$tpst+1;
				print OUT "$name\t$tpst\t$tped\t$count\t$len\n";
			}
			$name = $tp[$id];
			$tpst=$tp[$st];
			$tped=$tp[$ed];	
		}else{
			if($tp[$st] <= ($tped+$conn)){
				$tped=$tp[$ed];
			}else{
				$count++;
				my $len=$tped-$tpst+1;
				print OUT "$name\t$tpst\t$tped\t$count\t$len\n";
				$tpst=$tp[$st];
				$tped=$tp[$ed];	
			}
		}	
	}
	my $len=$tped-$tpst+1;
	print OUT "$name\t$tpst\t$tped\t$count\t$len\n";
	print "Connected to $count blocks\n";
}



sub ConnectSortedfasbyName{
	
	my $name="NA";  
	
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp= &Getheadinfo($lineIN, "|");
			if($name ne $tp[0]){
				if($name ne "NA"){
					print OUT "\n";
				}
				$name = $tp[0];
				print OUT "$name\n";
			}
			
		}else{
			print OUT "$lineIN";
		}
	}	
}


sub ConnectfasbyName_sorted{
	
	my $name="";  
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN,$sep);
		if($lineIN=~/>/){
			if($tp[0] eq $name){
				#same
			}else{
				$name=$tp[0];
				print OUT "$name\n";
			}			
		}else{
			print OUT "$lineIN\n";
		}
	}	

}


sub ConnectfasbyName{
	
	my @name;  
	my @seq;
	my $count=0;
	my $wh=0;
	my $sep="|";
	my $loc=0;
	my $insert="";
	if(length($_[1])>0){
		$insert=$_[1];
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN,$sep);
		my $tp=@tp;
		if($lineIN=~/>/){
			$wh=-1;
			for(my $i=0;$i<$count;$i++){
				if($tp[$loc] eq $name[$i]){
					$wh=$i;
					last;
				}
			}
			if($wh==-1){
				$name[$count]=$tp[$loc];
				#print "$tp[$loc]\n";
				$wh=$count;
				$seq[$wh]="";
				$count++;
			}
			
		}else{
			if(length($lineIN)>0){
				$seq[$wh]=$seq[$wh].$insert.$tp[0];
			}
		}
	}	
	for(my $i=0;$i<$count;$i++){
		print OUT "$name[$i]\n$seq[$i]\n";
	}
}




sub gff_2_gtf_for_tophat{
	
	#gff
	#NW_002477246.1	RefSeq	exon	1994433	1994530	.	-	.	ID=id2062;Parent=rna737;Dbxref=Genbank:XM_002382194.1,GeneID:7913739;Note=transcript AFLA_130420A;gbkey=mRNA;partial=true;product=ftsj%2C putative;transcript_id=XM_002382194.1
	#Chr2_A_oryzae_RIB40	AspGD	gene	350	1289	.	+	.	ID=AO090001000001;Name=AO090001000001;Note=Ortholog%28s%29%20have%20role%20in%20austinol%20biosynthetic%20process%2C%20dehydroaustinol%20biosynthetic%20process;orf_classification=Uncharacterized;Alias=AO070180000001
	#Chr2_A_oryzae_RIB40	AspGD	mRNA	350	1289	.	+	.	ID=AO090001000001-T;Parent=AO090001000001;Name=AO090001000001;Note=Ortholog%28s%29%20have%20role%20in%20austinol%20biosynthetic%20process%2C%20dehydroaustinol%20biosynthetic%20process;orf_classification=Uncharacterized;Alias=AO070180000001
	#Chr2_A_oryzae_RIB40	AspGD	exon	350	481	.	+	.	ID=AO090001000001-T-E1;Parent=AO090001000001-T
	#Chr2_A_oryzae_RIB40	AspGD	CDS	350	481	.	+	0	ID=AO090001000001-P;Parent=AO090001000001-T;orf_classification=Uncharacterized;parent_feature_type=ORF


	
	#gtf
	#NW_002477237.1	Cufflinks	exon	905864	905977	.	-	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1";


	my $exon_number=0;	
	my @out;
	
	my $dir="+";
	my $gene_id="";
	my $transcript_id="";
	my $last_transcript_id="NA";
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN,";");
		my $tp=@tp;
		if($tp[2] eq "gene" || $tp[2] eq "pseudogene"){
			if($tp[8]=~/ID=/){			
				$gene_id=substr($tp[8],3,length($tp[8]));
			}else{
				print "err: $lineIN\n";
			}
		}
		elsif($tp[2] eq "mRNA" || $tp[2] eq "rRNA" || $tp[2] eq "tRNA" || $tp[2] eq "ncRNA"){
			if($tp[8]=~/ID=/){			
				$transcript_id=substr($tp[8],3,length($tp[8]));		
			}else{
				print "err: $lineIN\n";
			}	
			
			if($transcript_id ne $last_transcript_id){
				if($last_transcript_id ne "NA"){
					#print OUT
					#print ".";
					if($dir eq "+"){
						for(my $i=0;$i<$exon_number;$i++){
							my $exon=$i+1;
							print OUT "$out[$i]$exon\";\n";
						}
					}else{
						for(my $i=$exon_number;$i>0;$i--){
							my $exon=$exon_number-$i+1;
							print OUT "$out[$i-1]$exon\";\n";
						}						
					}
				}
				#clean
				$exon_number=0;
				@out=();
				$last_transcript_id=$transcript_id;
				$dir=$tp[6];
			}
		}
		elsif($tp[2] eq "exon"){
			my $parent="";
			if($tp[9]=~/Parent=/){			
				$parent=substr($tp[9],7,length($tp[9]));
			}
			if($parent eq $transcript_id){
				#add
				$out[$exon_number]="$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\"; exon_number \"";		
				$exon_number++;
			}else{
				print "err: $parent $transcript_id $lineIN\n";
			}
		}
	}
	
	#last line
	if($dir=="+"){
		for(my $i=0;$i<$exon_number;$i++){
			my $exon=$i+1;
			print OUT "$out[$i]$exon\";\n";
		}
	}else{
		for(my $i=$exon_number;$i>0;$i--){
			my $exon=$exon_number-$i+1;
			print OUT "$out[$i-1]$exon\";\n";
		}						
	}		
}



sub gff_trans_2_pro_fas{
	
#gff_seqname	artemis	CDS	5543137	5545821	0	+	0	ID=ECSP_6001;gene=tagA;locus_tag=ECSP_6001;transl_table=11;product=hypothetical+protein;protein_id=YP_003082135.1;db_xref=GI:254667449%2CGeneID:8206983;translation=MNERWRTPMKLKYLSCTILAPLAIGVFSATAADNNSAIYFNTSQ+PINDLQGSLAAEVKFAQSQILPAHPKEGDSQPHLTSLRKSLLLVRPVKADDKTPVQVE+ARDDNNKILGTLTLYPPSSLPDTIYHLDGVPEGGIDFTPHNGTKKIINTVAEVNKLSD+ASGSSIHSHLTNNALVEIHTANGRWVRDIYLPQGPDLEGKMVRFVSSAGYSSTVFYGD+RKVTLSVGNTLLFKYVNGQWFRSGELENNRITYAQHIWSAELPAHWIVPGLNLVIKQG+NLSGRLNDIKIGAPGELLLHTIDIGMLTTPRDRFDFAKDKEAHREYFQTIPVSRMIVN+NYAPLHLKEVMLPTGELLTDMDPGNGGWHSGTMRQRIGKELVSHGIDNANYGLNSTAG+LGENSHPYVVAQLAAHNSRGNYANGIQVHGGSGGGGIVTLDSTLGNEFSHEVGHNYGL+GHYVDGFKGSVHRSAENNNSTWGWDGDKKRFIPNFYPSQTNEKSCLNNQCQEPFDGHK+FGFDAMAGGSPFSAANRFTMYTPNSSAIIQRFFENKAVFDSRSSTGFSKWNADTQEME+PYEHTIDRAEQITASVNELSESKMAELMAEYAVVKVHMWNGNWTRNIYIPTASADNRG+SILTINHEAGYNSYLFINGDEKVVSQGYKKSFVSDGQFWKERDVVDTREARKPEQFGV+PVTTLVGYYDPEGTLSSYIYPAMYGAYGFTYSDDSQNLSDNDCQLQVDTKEGQLRFRL+ANHRANNTVMNKFHINVPTESQPTQATLVCNNKILDTKSLTPAPEGLTYTVNGQALPA+KENEGCIVSVNSGKRYCLPVGQRSGYSLPDWIVGQEVYVDSGAKAKVLLSDWDNLSYN+RIGEFVGNVNPADMKKVKAWNGQYLDFSKPRSMRVVYK


	my $translation="";
	my $gene="";
	my $product="";
	my $gi="";
	my $gb="";
	my $species="";
	while($lineIN=<IN>){
		my @tp= &Getheadinfo($lineIN,";");
		my $tp=@tp;
		if($tp[2] eq "CDS"){
			for(my $i=8;$i<$tp;$i++){
				if($tp[$i]=~/gene=/){
					$gene=substr($tp[$i],5,length($tp[$i]));
				}	
				elsif($tp[$i]=~/product=/){
					$product=substr($tp[$i],8,length($tp[$i]));
					#$product =~ s/+/_/;
				}	
				#elsif($tp[$i]=~/protein_id=/){
				#	$gb= "gb|".substr($tp[$i],11,length($tp[$i]))."|";
				#}	
				elsif($tp[$i]=~/db_xref=/){
					#db_xref=GI:254667449%2CGeneID:8206983
					$gi= "gi|".substr($tp[$i],11,(index($tp[$i],"%")-11))."|";
					$gb= "gb|".substr($tp[$i],(index($tp[$i],"GeneID:")+7),length($tp[$i]))."|";
				}	
				elsif($tp[$i]=~/translation=/){
					$translation=substr($tp[$i],12,length($tp[$i]));					
					$translation =~ tr/+//d
				}									
			}				
			print OUT ">$gi$gb $gene $product [$species]\n$translation\n";
		}	
	}	
}


	
sub gb_gbk_2_pro_fas{
	
	#only one chromsome, if more ,change this script
#     gene            complement(16476..16610)
#                     /locus_tag="EDL933_p0021"	
#                     /translation="MLILISPAKTLDYQSPLTTTRYTLPELLDNSQQLIHEARKLTPP
#                     WRFVEMGTVNMPPHPFVRPAFDVRSEQAAQVAIARMNRAIDEVLRR"
	my $gene="";
	my $product="";
	my $gb="";
	my $canadd=0;
	my $species="";
	while($lineIN=<IN>){
		
		my @tp= &Getheadinfo($lineIN);
		if($tp[0] eq "SOURCE"){
			$species=substr($lineIN,12,length($lineIN));
		}
		elsif($tp[0] eq "CDS"){
			$loc=$tp[1];
			$gene="";
			$product="";
			$gi="";
			$gb="";
		}
		elsif($lineIN=~/\/protein_id=\"/){
			$gb= "gb|".substr($tp[0],13,(length($tp[0])-1-13))."|";
		}		
		elsif($lineIN=~/\/db_xref=\"GI:/){
			$gi= "gi|".substr($tp[0],13,(length($tp[0])-1-13))."|";
		}			
		elsif($lineIN=~/\/gene=\"/){
			$gene= substr($tp[0],7,(length($tp[0])-1-7));
		}				
		elsif($lineIN=~/\/product=\"/){
			$product= substr($lineIN,31,length($lineIN));
			$product =~ tr/\"//d
		}		
		elsif($lineIN=~/\/translation=/){
			print OUT ">$gi$gb $gene $product [$species]\n";
			$tp[0]=substr($tp[0],14,(length($tp[0])));
			
			if($tp[0]=~/\"/){
				$tp[0]=substr($tp[0],0,(length($tp[0])-1));
				$canadd=0;
				print OUT "$tp[0]\n";  
			}else{
				$canadd=1;
				print OUT "$tp[0]";  
			}		
		}
		elsif($canadd==1){
			if($lineIN=~/\"/){   
				$tp[0]=substr($tp[0],0,(length($tp[0])-1));
				print OUT "$tp[0]\n";  
				$canadd=0;
			}else{
				print OUT "$tp[0]";
			}
		}
	}	
}


sub gb_gbk_2_genelocal{
	
	#only one chromsome, if more ,change this script
#     gene            complement(16476..16610)
#                     /locus_tag="EDL933_p0021"	
#                     /translation="MLILISPAKTLDYQSPLTTTRYTLPELLDNSQQLIHEARKLTPP
#                     WRFVEMGTVNMPPHPFVRPAFDVRSEQAAQVAIARMNRAIDEVLRR"
	my $gene="";
	my $product="";
	my $gb="";
	my $canadd=0;
	my $species="";
	my $loc="";
	while($lineIN=<IN>){
		my @tp= &Getheadinfo($lineIN);
		if($tp[0] eq "SOURCE"){
			$species=substr($lineIN,12,length($lineIN));
		}
		elsif($tp[0] eq "CDS"){
			$loc=$tp[1];
			$product="";
			$gi="";
			$gb="";
			$canadd=1;
		}	
		if($canadd == 1){
			if($tp[0] eq "gene"){
				#print OUT ">$gi$gb\t$loc\t$gene $product [$species]\n";
				print OUT ">$gi\t$loc\t$gene $product [$species]\n";
				$canadd=0;
			}
			#elsif($lineIN=~/\/db_xref=\"GeneID:/){
			#	$gb= "gb|".substr($tp[0],(index($tp[0],":")+1),length($tp[0]))."|";
			#	$gb =~ tr/\"//d;
			#}	
			elsif($lineIN=~/\/db_xref=\"GI:/){
				$gi= "gi|".substr($tp[0],(index($tp[0],":")+1),length($tp[0]))."|";
				$gi =~ tr/\"//d;
			}	
			elsif($lineIN=~/\/gene=\"/){
				$gene= substr($tp[0],7,length($tp[0]));
				$gene =~ tr/\"//d;
			}				
			elsif($lineIN=~/\/product=\"/){
				$product= substr($lineIN,31,length($lineIN));
				$product =~ tr/\"//d;
			}		
		}
	}	
	#print OUT ">$gi$gb\t$loc\t$gene $product [$species]\n";
	print OUT ">$gi\t$loc\t$gene $product [$species]\n";
}


sub ConnectfasbyNNN{
	
	my $n=50;
	if($_[1]>=0){
		$n=$_[1];		
	}
	$lineIN=<IN>;
	print OUT ">$ARGV[0]\n";
	while($lineIN=<IN>){
		#print "$lineIN";
		if($lineIN=~/>/){
			for(my $i=0; $i<$n;$i++){
				print OUT "-";
			}
			#print OUT "\n";
		}else{
			print OUT "A$lineIN";
		}
	}
	
}


sub CombinCorrespondenceElements{
	
	#@db [0]-number [1]string elements
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= split(/\t/,$lineIN);
		my $tp=@tp;
		my $a1=-1;
		my $a2=-1;
		for(my $i=0;$i<$count;$i++){
			my @a=split(/\t/,$db[$i][1]);
			for(my $j=1;$j<$db[$i][0];$j++){
				
				if($a[$j] eq $tp[0]){
					$a1=$i;
				}
				if($a[$j] eq $tp[1]){
					$a2=$i;
				}
			}		
		}
		if($a1>=0 && $a2>=0){
			if($a1 != $a2){
				$db[$a1][0]=$db[$a2][0]+$db[$a1][0];
				$db[$a1][1]=$db[$a1][1]."\t".$db[$a2][1];
				$db[$a2][0]=0;
			}
		}elsif($a1>=0){
			$db[$a1][0]++;
			$db[$a1][1]=$db[$a1][1]."\t".$tp[1];
		}elsif($a2>=0){
			$db[$a2][0]++;
			$db[$a2][1]=$db[$a2][1]."\t".$tp[0];
		}else{
			$db[$count][0]=2;
			$db[$count][1]=$tp[0]."\t".$tp[1];
			$count++;
		}
	}
	
}

sub CompareValuestoCutoff{
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= split(/\t/,$lineIN);
		my $tp=@tp;
		my $n=0;
		for(my $i=0;$i<$tp;$i++){
			if($_[1] eq "gt"){
				if($tp[$i] > $_[2]){
					$n++;
				}
			}elsif($_[1] eq "lt"){
				if($tp[$i] < $_[2]){
					$n++;
				}				
			}elsif($_[1] eq "eq"){
				if($tp[$i] == $_[2]){
					$n++;
				}				
			}else{
			
			}
		}
		print OUT "$n\n";
	}
	print OUT "$seq\n";
}






sub CombinFasbyID{
	
	#>ID_other detail
	my $seq="NA";
	my $name="NA";
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			#new
			my @tp= &Getheadinfo($lineIN, "|");
			#print "$tp[0]\n";
			if($name ne $tp[0]){
			if($name ne "NA"){
				#out put
				print OUT "$seq\n";
			}				
				print OUT "$tp[0]\n";	
				$name=$tp[0];
				$seq="";
			}
		}else{
			$seq=$seq.$lineIN;
		}
	}
	print OUT "$seq\n";
}


sub ListLinetypevalue{
	
	#8_S_35:	0M	2	1M	3	2M	12	3M	37	4M	391	5M	625	
	
	my $len=@_;
	my @key;
	print OUT "ID";
	for(my $i=1;$i<$len;$i++){
		$key[$i]=$_[$i];
		print OUT "\t$_[$i]";
	}
	print OUT "\n";
	while($lineIN=<IN>){
			chomp($lineIN);
			for(my $i=1;$i<$len;$i++){
				$db{$key[$i]}=0;
			}
			my @tp= &Getheadinfo($lineIN);
			my $tp=@tp;
			for(my $i=1;$i<$tp;$i+=2){
				$db{$tp[$i]}=$tp[$i+1];
				#print "$tp[$i] $db{$tp[$i]}\t";
			}			
			print OUT "$tp[0]";
			for(my $i=1;$i<$len;$i++){
				print OUT "\t$db{$key[$i]}";
				#print ".";
			}
			print OUT "\n";
	}

}


sub CombinAlternativeSplicingCDSloc{
	
#	chr1	9475309	9487684	+	ENSSSCG00000029876	3	ENSSSCT00000031745	9475309	9487684	5	9475629	9475651	9476006	9476208	9480034	9480150	9482724	9482903	9484697	9484839
#	chr1	9475309	9487684	+	ENSSSCG00000029876	3	ENSSSCT00000034561	9475621	9485188	4	9475629	9475651	9476006	9476208	9482724	9482903	9484697	9484839
#	chr1	9475309	9487684	+	ENSSSCG00000029876	3	ENSSSCT00000036564	9475621	9485216	4	9475629	9475651	9476006	9476208	9480034	9480150	9484697	9484839

	my @cb;
	my $cbn=0;
	my $name="NA";
	my @oc;
	my $ocn=0;
	my $st=0;
	my $ed=0;


	while($lineIN=<IN>){
			chomp($lineIN);
			my @tp= &Getheadinfo($lineIN);
			my $tp=@tp;
			if($name ne $tp[4]){
				print "+";
				if($name ne "NA"){
					#output
					my $addnew=0;   #0-inter  1-exon
					print "$cbn";
					for(my $i=$st;$i<=$ed;$i++){
						my $add=0;  #0-no 1-yes						
						for(my $j=0;$j<$cbn;$j++){
							#print "\t$cb[$j][0]";
							for(my $k=1;$k<(2*$cb[$j][0]);$k+=2){
								if($i>=$cb[$j][$k] && $i<=$cb[$j][$k+1]){
									$add=1;
									#print ".";
									last;
								}
							}
							if($add==1){
								last;
							}
						}
						if($addnew==0){
							if($add==1){
								$oc[$ocn][0]=$i;
								$oc[$ocn][1]=$i;
								$addnew=1;
							}
						}else{
							if($add==1){
								$oc[$ocn][1]=$i;
							}else{
								$addnew=0;
								$ocn++;
							}							
						}
					}
					if($addnew==1){
						
						$ocn++;
						
					}
					#print "\t$ocn\t$oc[0][0]\n";
					print OUT "$ocn";
					for(my $i=0;$i<$ocn;$i++){
						print OUT "\t$oc[$i][0]\t$oc[$i][1]";
					}
					print OUT "\n";
				}	
				#new
				print OUT "$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\tCombin\t$tp[1]\t$tp[2]\t";								
				$name= $tp[4];
				$st=$tp[1];
				$ed=$tp[2];
				@cb=();
				@oc=();
				$ocn=0;
				$cbn=0;										
			}
			$cb[$cbn][0]=$tp[9];
			my $c=1;
			for(my $i=10;$i<$tp;$i++){
				$cb[$cbn][$c]=$tp[$i];
				$c++;					
			}				
			$cbn++;	
	}
	
	#last output
	my $addnew=0;   #0-inter  1-exon
	print "\t$cbn\t$cb[0][0]";
	for(my $i=$st;$i<=$ed;$i++){
		my $add=0;  #0-no 1-yes						
		for(my $j=0;$j<$cbn;$j++){
			#print "\t$cb[$j][0]";
			for(my $k=1;$k<(2*$cb[$j][0]);$k+=2){
				if($i>=$cb[$j][$k] && $i<=$cb[$j][$k+1]){
					$add=1;
					#print ".";
					last;
				}
			}
			if($add==1){
				last;
			}
		}
		if($addnew==0){
			if($add==1){
				$oc[$ocn][0]=$i;
				$oc[$ocn][1]=$i;
				$addnew=1;
			}
		}else{
			if($add==1){
				$oc[$ocn][1]=$i;
			}else{
				$addnew=0;
				$ocn++;
			}							
		}
	}
	if($addnew==1){
		$ocn++;
		
	}
	print "\t$ocn\t$oc[0][0]\n";
	print OUT "$ocn";
	for(my $i=0;$i<$ocn;$i++){
		print OUT "\t$oc[$i][0]\t$oc[$i][1]";
	}
	print OUT "\n";

	print "\n";	
}


sub ChangeGffouttoLoc{
	#
	#1	1979288	1987260	+	ENSSSCG00000027112	1	ENSSSCT00000004438	1979288	1987260	3	1979288	1979325	1980125	1980399	1987208	1987260

	#CDS
	#JZDT01000187.1  841     1875    gene0   rna0    +       0       4       0       900     1085    1146    1508    1568    1612    1666    1875
	#NW_002477246.1  100227  101173  gene42  rna42   +       0       3       0       100227  100436  100534  100894  101016  101173
	
	#gene
	#000000F_arrow_pilon     1941843 1943757 +      gene-19.4    4       1941843 1941903 1942253 1942410 1942428 1943083 1943556 1943757
	

	#complement(join(2354845..2355351,2355357..2355806))
	my $name=substr($ARGV[0],0,index($ARGV[0],".",));
		
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;		
		if($_[1] eq "CDS"){			
			my $o="$name=$tp[3]\_$tp[4]|$tp[7]|$tp[0]:";	
			if($tp[5] eq "-"){
				$o=$o."complement(";
			}		
			if($tp[7]==1){
				$o=$o."$tp[9]..$tp[10]";

			}else{				
				$o=$o."join($tp[9]..$tp[10]";	
				for(my $i=11;$i<$tp;$i+=2){		
					$o=$o.",$tp[$i]..$tp[$i+1]";
				}
				$o=$o.")";
			}
			if($tp[5] eq "-"){
				$o=$o.")";
			}			
			print OUT "$o\n";			
		}elsif($_[1] eq "gene"){
			my $c=0;
			for(my $i=6;$i<$tp;$i+=2){
				$c++;
				if($tp[3] eq "+"){
					print OUT ">$tp[4]|$c/$tp[5]|$tp[0]:$tp[$i]..$tp[$i+1]\n";
				}else{
					print OUT ">$tp[4]|$c/$tp[5]|$tp[0]:complement($tp[$i]..$tp[$i+1])\n";
				}
			}	
		}else{
			my $c=0;
			for(my $i=9;$i<$tp;$i+=2){
				$c++;
				if($tp[3] eq "+"){
					print OUT ">$tp[4]\_$tp[6]|$c/$tp[9]|$tp[0]:$tp[$i]..$tp[$i+1]\n";
				}else{
					print OUT ">$tp[4]\_$tp[6]|$c/$tp[9]|$tp[0]:complement($tp[$i]..$tp[$i+1])\n";
				}
			}
		}			
	}
}


sub GetGeneExpfromSoftFile{

	#read gene ID
	my $max=10;   #max sample
	my @gene;  		#0-ID 1-ORF 2-info
	my %geneID;
	my $count=0;  #gene count;
	my $add=0;		#0-no add  1-add gene   2-add sample
	#get sample name
	my $sc=0;     #sample number
	my @sampleName;     #sample name
	#get expression data
	while($lineIN=<IN>){
		chomp($lineIN);

		if($add==0){
			if(index($lineIN,"platform_table_begin",0)>=0){
				$lineIN=<IN>;
				$add=1;
			}			
		}elsif($add==1){
			if(index($lineIN,"platform_table_end",0)>=0){
				$add=2;
			}else{
				my @tp= &Getheadinfo($lineIN);
				$gene[$count][0]=$tp[0];
				$gene[$count][1]=$tp[1];	
				$gene[$count][2]=$lineIN;
				for(my $i=3;$i<$max;$i++){
					$gene[$count][$i]=0;
				}	
				$geneID{$tp[0]}=$count;		
				$count++;
			}										
		}elsif($add==2){
			if(index($lineIN,"Sample_title",0)>=0){
				my @tp= &Getheadinfo($lineIN);
				$sampleName[$sc]=$tp[2];
			}
			if(index($lineIN,"sample_table_begin",0)>=0){
				$lineIN=<IN>;
				$add=3;
			}							
		}elsif($add==3){	
			if(index($lineIN,"sample_table_end",0)>=0){
				$add=2; 
				$sc++;
			}else{	
				my @tp= &Getheadinfo($lineIN);
				my $loc=$geneID{$tp[0]};
				$gene[$loc][$sc+3]=$tp[1];	
			}
		}

	}
	
	#output
	print OUT "Gene_ID\tORF_ID";  
	for(my $j=0;$j<$sc;$j++){
		print OUT "\t$sampleName[$j]";
	}
	print OUT "\t Detail\n";
	for(my $i=0;$i<$count;$i++){
		print OUT "$gene[$i][0]\t$gene[$i][1]";    
		for(my $j=3;$j<(3+$sc);$j++){
			print OUT "\t$gene[$i][$j]";
		}
		print OUT "\t$gene[$i][2]\n";
		#print OUT "\n";
	}
	
}


sub getLineFromSymble{

	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/$_[1]/){
			print OUT $lineIN."\n";
		}
	}
}

sub getLineWithoutSymble{

	while($lineIN=<IN>){
		chomp($lineIN);
		if(index($lineIN,$_[1],)<0){
			print OUT $lineIN."\n";
		}
	}
}

sub getLinfromName{
	print "$_[1] $_[2]\n";
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		#my @tp = split(/\t/,$lineIN);
		#print "$tp[$_[2]]\t$_[1]\n";
		if($tp[$_[2]] eq $_[1]){
			print OUT $lineIN."\n";
		}
	}
}

sub getLinWithoutName{

	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($tp[$_[2]] ne $_[1]){
			print OUT $lineIN."\n";
		}		
	}
}

sub getLinefromNumber{
	my $end=$_[1]+$_[2];
	my $count=0;
	my $canadd=0;
	while($lineIN=<IN>){
		$count++;	
		if($count==$_[1]){
			$canadd=1;
		}	
		if($count==$end){
			last;
		}		
		if($canadd==1){
			print OUT "$lineIN";
		}				
	}
}

sub SeperateFilebyLine{
	
	my $c=1;
	my $outt=$ARGV[0].".$c.out";
	open(OUTT,'>',"$outt");
	print OUT "$outt\n";
	print "$outt\n";
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		$count++;	
		print OUTT "$lineIN\n";
		if($count==$_[1]){
			close(OUTT);
			$count=0;
			$c++;
			$outt=$ARGV[0].".$c.out";
			print OUT "$outt\n";
			print "$outt\n";
			open(OUTT,'>',"$outt");
		}		
				
	}
	close(OUTT);
}

sub FindLonglinebyCutoff{
		my $count=0;
		while($lineIN=<IN>){
		chomp($lineIN);
		$count++;	
		my $len=length($lineIN);		
		if($len > $_[1]){
			print OUT "$count $len $lineIN\n";
			print "$count $len $lineIN\n";
		}
	}
}




sub getlocMapEverageDepth{
	#depth cannot change line
	my $m1=0; #media
	my $m2=0; #modal
	my $m3=0; #perent $m3n average
	my $m3n=50;	 #percent of length	
	my $m4=0; #perent $m4n average
	my $m4n=80;  #percent of length	
	my $n=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		$n++;
		if($n==10000000){
			$n=0;
			print ".";
		}
		if($lineIN=~/>/){
			print OUT "$lineIN";
		}else{
			my @tp = split(/ /,$lineIN); 
			my $tp=@tp;			
			my $t1=0;
			my %t2;
			for(my $i=0;$i<$tp;$i++){
				$t1+=$tp[$i];
				if($t2{$tp[$i]}){
					$t2{$tp[$i]}=$t2{$tp[$i]}+1;
				}else{
					$t2{$tp[$i]}=1;
				}
			}
			$m1=0;
			$m2=0;	
			$m3=0;
			$m4=0;	
			my $m3len=$tp*$m3n/100;
			$m3len=&Floatcut($m3len,0);	
			my $m4len=$tp*$m4n/100;
			$m4len=&Floatcut($m4len,0);
			#print "$m4len $tp\n";									
			if($tp>0){								
				$m1=$t1/$tp;
				my @k=keys(%t2);
				my $k=@k;	
				my @value;
				for(my $i=0;$i<$k;$i++){
					$value[$i]="$t2{$k[$i]} $k[$i]";
					my $len=length($t2{$k[$i]});
					for($j=$len;$j<=20;$j++){
						$value[$i]="0".$value[$i];
					}	
					#print "$value[$i]\n";								
				}								
				@value=sort(@value);
				my @value2=reverse (@value);
				my $value=@value;
				
				my $n=0;
				my $add=1;
				for(my $i=0;$i<=$value;$i++){
					my @tpm= split(/ /,$value2[$i]); 
					if($i==0){
						$m2=$tpm[1];
					}
					if($n>$m3len && $add==1){
						$add=0;
							$m3=$m4/$n;
					}
					if($n>$m4len){
						$m4=$m4/$n;
						last;
					}	
					$m4=$m4+($tpm[1]*$tpm[0]);
					$n=$n+$tpm[0];
				}
				#print "$value2[0] $value2[1] $value2[2] $value2[3]\n";
			}
			$m1=&Floatcut($m1,2);
			$m3=&Floatcut($m3,2);
			$m4=&Floatcut($m4,2);
			print OUT "\t$tp\t$m1\t$m4\t$m3\t$m2\n";
		}		
	}
	print "\n";
}

sub Floatcut{
	#print "$_[0]\t";
	my $o=$_[0];
	my $l=index($_[0],".",0);
	my $leave=0;
	if($_[1]>0){
		$leave=$leave+1+$_[1];
	}
	if($l>0){
		my $a=substr($_[0],$l+$leave+1,1);
		$o=substr($_[0],0,$l+$leave);
		if($a<5){		
		}else{
			$o++;
		}
	}	
	return $o;	
}



sub BatchCombinFilebyRowName{

	#read in list
	#CK_24=CK_48.tsv	CK_24_R1.reads.count	CK_24_R2.reads.count	CK_24_R3.reads.count	CK_48_R1.reads.count	CK_48_R2.reads.count	CK_48_R3.reads.count
	#CK_24=CK_72.tsv	CK_24_R1.reads.count	CK_24_R2.reads.count	CK_24_R3.reads.count	CK_72_R1.reads.count	CK_72_R2.reads.count	CK_72_R3.reads.count
	
	#each file must have same elements
	
	my $incount=0;

	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		
		my $outn=$tp[0].".combin";		
		open(OUTN,'>',"$outn");	
		print OUTN "ID\t$tp[0]";
		$incount++;
		print OUT "file $incount: $outn\n";		
		#data in %db
		&Filetohash($tp[0]);  #read first file to hash %db
		
		for(my $i=1;$i<$tp;$i++){
			print OUTN "\t$tp[$i]";
			&Addinfo($tp[$i]);
		}
		print OUTN "\n";		
		
		#output 

		my @outkey=sort(keys(%db));
		my $outkey=@outkey;
		for(my $i=0;$i<$outkey;$i++){
			print OUTN "$outkey[$i]";
			my @tpout= &Getheadinfo($db{$outkey[$i]});
			my $tpout=@tpout;
			for(my $j=0;$j<$tpout;$j++){
				print OUTN "\t$tpout[$j]";
			}
			print OUTN "\n";
		}
		close(OUTN);	
		
	}
}

sub Addinfo{   #0-add file2 1-row number
	
	open(INB,'<',"$_[0]");	
	my $count=0;
	while($lineINB=<INB>){
		chomp($lineINB);
		my @tpinb = &Getheadinfo($lineINB); 
		my $tpinb=@tpinb;
		my $r="";
		if($db{$tpinb[0]}){
        	$r = "\t".$tpinb[1];
        	#for ($i=2; $i<$tpinb; $i++)
        	#{
        	#	$r = $r."\t".$tpinb[$i];
       		#}
       		#print "$db{$tpinb[0]}\t";
       		$db{$tpinb[0]}=$db{$tpinb[0]}.$r;
       		#print "$db{$tpinb[0]}\n";
		}else{
			print OUT "err:$tpinb[0]\n";
		}
		$count++;
	}	
	print OUT "$_[0]\t$count\n";	
	close(INB);
}

sub Filetohash{	  #0-file name
	#first row to key
	#other row to value
	open(INTP,'<',"$_[0]");	
	my @tpinfo;
	my $count=0;
	while($lineINTP=<INTP>){
		chomp($lineINTP);
		@tpinfo = ('');
		$tpvalue = '';
        @tpinfo = &Getheadinfo($lineINTP);
        $tpinfo = @tpinfo;
        #print $tpinfo[0]." ".$tpinfo[1]." ".$tpinfo."\n";
        $tpvalue = $tpinfo[1];
        for ($i = 2; $i < $tpinfo; $i++)
        {
        	$tpvalue = $tpvalue."\t".$tpinfo[$i];
       	}
       	if (!$tpvalue){
       		$db{$tpinfo[0]} = "x";	
       	}else{
        	$db{$tpinfo[0]} = $tpvalue;
        }
        $count++;
	}
	print OUT "$_[0]\t$count\n";
	close(INTP);
}


sub BatchAverageValueCount{

#Compound	1-1	1-2	2-1	2-2	3-1	3-2	4-1
#M1	1.032278745	0.386359919	0.13262574	0.699764301	0.28544503	0.196760573	0.501724272
#M2	0.782141996	0.634079615	0.566509371	0.599049721	0.957193281	0.381254052	0.560427154

	my $outn=$ARGV[0].'.data';
	open(OUTN,'>',"$outn");
	

	#head
	$lineIN=<IN>;
	chomp($lineIN);	
	#@tp = &Getheadinfo($lineIN);
	my @head=split(/\t/,$lineIN);	
	my $head=@head;
	my $count=0;
	
	#read in file
	#@data=();  #contail value
	while($lineIN=<IN>){
		chomp($lineIN);	
		#my @tp = &Getheadinfo($lineIN);
		my @tp=split(/\t/,$lineIN);
		$data[$count]=[@tp];
		$count++;
		#print $data[0][0]."\6";
	}	
	print "line $count\n";
	
	#out
	my @avgv=();
	my @avgname=();
	my $avgno=0;
	
	#temp
	my @add=();
	my $tpno=1;
	my $name="NA";
	my $last=$head-1;
	for(my $i=1;$i<$head;$i++){
		my @tpname=split(/-/,$head[$i]);
		if($tpname[0] ne $name || $i==$last){
			if($name ne "NA"){
				#add				
				$avgname[$avgno]=$name;
				#count each row
				for(my $j=0;$j<$count;$j++){
					my $totalv=0;
					print OUTN "$avgname[$avgno]\t$data[$j][0]";
					for(my $k=0;$k<$tpno;$k++){
						$totalv+=$data[$j][$add[$k]];
						print OUTN "\t$data[$j][$add[$k]]";
					}	
					print OUTN "\n";
					#avg
					$avgv[$j][$avgno]=$totalv/$tpno;
					
				}
				$avgno++;
				print "$name $tpno\n";
			}		
			
			#clean
			@add=();
			$add[0]=$i;
			$tpno=1;
			$name=$tpname[0];

		}else{
			$add[$tpno]=$i;
			$tpno++;
		}
	}
	#out put head
	print OUT "id";
	for(my $i=0;$i<$avgno;$i++){
		print OUT "\t$avgname[$i]";
	}
	print OUT "\n";
	#out put avg value
	for(my $i=0;$i<$count;$i++){
		print OUT "$data[$i][0]";
		for(my $j=0;$j<$avgno;$j++){
			print OUT "\t$avgv[$i][$j]";
		}
		print OUT "\n";
	}	

	close(OUTN);
	print "avg no $avgno\ndone\n";

}



sub BatchAlntoFasbyList{

	#read in list

	while($lineIN=<IN>){
		chomp($lineIN);	
		my $outn=$lineIN.".aln.fas";
		my $inn=$lineIN.".aln";
		print "$outn ";
		open(OUTN,'>',"$outn");	
		open(INN,'<',"$inn");
		$lineINTP=<INN>;
		$lineINTP=<INN>;
		$lineINTP=<INN>;
		my %seq;
		my $count=0;
		my $add=1;
		while($lineINTP=<INN>){
			#print $lineINTP."\t";
			chomp($lineINTP);	
			my @tp= &Getheadinfo($lineINTP);	
			$tp=@tp;	
			if($add==1){
				if((index($lineINTP,'*',0)>=0) || $tp == 0){
					$add=0;					
				}else{
					$seq{$tp[0]}=$tp[1];
					$count++;
				}
			}else{
				if($seq{$tp[0]}){
					$seq{$tp[0]}=$seq{$tp[0]}.$tp[1];
				}
			}
		}
		@o=keys(%seq);
		@o=sort(@o);
		#print OUT "$count\n";
		for(my $i=0;$i<$count;$i++){
			print OUT ">$o[$i]\n$seq{$o[$i]}\n";
			print OUTN ">$o[$i]\n$seq{$o[$i]}\n";
		}		
		close(OUTN);	
		close(INN);
		print OUT "\n";
		
	}
}

sub BatchClustalWbyList{

	#read in list

	while($lineIN=<IN>){
		chomp($lineIN);	
		
		my $aln=$lineIN; #input
		my $clustalw = "ClustalW2"; #clustalw local
		my $tempoutput=$lineIN.'.aln'; #output file
		my $system_check=system("$clustalw -KTUPLE=2 -INFILE=$aln -OUTFILE=$tempoutput"); 		
		
	}
}



sub CombinFilesbyList{

	#read in list

	while($lineIN=<IN>){
		chomp($lineIN);	
		open(OUTN,'<',"$lineIN");			
		while($lineINTP=<OUTN>){
			print OUT "lineINTP";
		}
		print "$lineIN\n";
		close(OUTN);

	}
}


sub ChangeSoapsnpConsensus_toFas{
	##Chr10   3367    A       A       1       A       0       0       0       T       0       0       0       0       1.00000 255.000 0
	my $n=3;
	if($_[1]>0){
		$n=$_[1];
	}
	my $seq="";
	my $name="";
	my $l=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);	
		$tp=@tp;
		if($tp[0] ne $name){
			if($name ne ""){
				#output
				$l=length($seq);
				for(my $i=0;$i<1;$i+=100){
					my $o=substr($seq,$i,100);
					print OUT "$l\n$o\n";
					print "$l\n";
				}
			}
			$name=$tp[0];
			print OUT ">$name\t";
			print "$name\t";
			$seq="";
		}
		$seq=$seq.$tp[3];	
	}
	#last one
	$l=length($seq);
	for(my $i=0;$i<1;$i+=100){
		my $o=substr($seq,$i,100);
		print OUT "$l\n$o\n";
	}
	print "$l\n";
}


sub Change_Soapsnp_Consensus_toFas{
	##Chr10   3367    A       A       1       A       0       0       0       T       0       0       0       0       1.00000 255.000 0
	my $n=3;
	if($_[1]>0){
		$n=$_[1];
	}
	my $seq="";
	my $name="";
	my $l=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);	
		$tp=@tp;
		if($tp[0] ne $name){
			$name=$tp[0];
			print OUT ">$name\n";
			print "$name\t";
		}
		print OUT "$tp[3]";	
	}
	print "\n";
}


sub CountSoapCoverageWindowsAvgDepth{
	
	my $len=0;
	my $dep=0;
	my $cut=10000;
	my $count=0;
	my $cn=1;
	if($_[1]>0){
		$cut=$_[1];
	}
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($lineIN=~/>/){
			print OUT "$lineIN\n";
		}else{
			my @tp = split (' ', $lineIN);
			my $tp=@tp;
			for(my $i=0;$i<$tp;$i++){
				$count++;
				$len=$len+$tp[$i];
				if($count==$cut){
					$dep=$len/$count;
					my $o=$cn*$count;
					$cn++;
					$len=0;
					$count=0;
					print OUT "$o\t$dep\n";
				}
			}
		}
	}
	$dep=$len/$count;
	my $o=$cn*$count;
	$cn++;
	$len=0;
	$count=0;
	print OUT "$o\t$dep\n";	
}

sub count_col_average_value{
	#
	my $count=0;
	my $value=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		#my @tp= &Getheadinfo($lineIN);	
		my @tp= split(/\t/,$lineIN);
		$count++;
		$value+=$tp[$_[1]];
	}
	my $avg=$value/$count;
	print OUT "$avg\n";
}


sub CountSoapCoverageAvgDepth{
	#scaffold385: 143786/170737	Percentage:84	Depth:7.9
	my $len=0;
	my $dep=0;
	my $cut=10000;
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN,"/",":");	
		if($lineIN=~/Depth/){
			if($tp[1]>$cut && $tp[6] ne "nan"){
				#print "$tp[1]\t$tp[6]\n";
				#print ".";
				$len+=$tp[1];
				$dep=$dep+($tp[1]*$tp[6]);
			}
		}
	}
	my $o=0;
	if($len>0){
		$o=$dep/$len;
	}
	print OUT "Len\t$len\nDepth\t$o\n";
	print "Len\t$len\tDepth\t$o\n";
}


sub CountlineElementsNo{
	#

	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);	
		$tp=@tp;
		print OUT "$tp\n";
	}
	print "\n";
}

sub CombinBlastclustalforBothLTR{
	#combin different class for LTR5 and LTR3
	#readin and loop combin
	#output remove 5 and 3
	#each line had a value for type number 
	
	#name length
	my $l=8;
	if($_[1]>0){
		$l=$_[1];
	}
	print "name length = $l\n";
	my $outfname=$ARGV[0].".combin";
	open(OUTF,'>',"$outfname");	
	
	#read in
	my %hash;	
	my $count=1;
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);			
		$tp=@tp;
		my $value=$count;

		#test if same as before
		my @cover=();
		my $coverno=0;
		my $tpvalue=0;
		my @which=();		
		for(my $i=0; $i<$tp; $i++){
			if($hash{$tp[$i]}>0){
				if($tpvalue != $hash{$tp[$i]}){
					$tpvalue=$hash{$tp[$i]};
					$cover[$coverno]=$hash{$tp[$i]};
					$which[$coverno]=$tp[$i];
					$coverno++;
				}
			}
		}
		
		#add new count
		if($coverno ==0){
			print "+";
			for(my $i=0; $i<$tp; $i++){
				$hash{$tp[$i]}=$count;				
			}	
			$count++;		
		}
		
		#if same as one class
		elsif($coverno == 1){
			print ".";
			for(my $i=0; $i<$tp; $i++){
				$hash{$tp[$i]}=$tpvalue;				
			}			
		}
		
		#if cover muti class, combin to first one
		else{
			print "-";
			print OUTF "$coverno $which[0]";
			
			for(my $i=0; $i<$tp; $i++){
				$hash{$tp[$i]}=$cover[0];				
			}
			for(my $i=1;$i<$coverno;$i++){
				print OUTF " $which[$i]";
				#find all valus same
				my @hkeys=keys(%hash);
				my $hkeys=@hkeys;
				for(my $j=0;$j<$hkeys;$j++){
					if($hash{$hkeys[$j]}==$i){
						$hash{$hkeys[$j]}=$cover[0];
					}
				}
			}	
			print OUTF "\n";					
		}
	}	
	
	#class
	my %n;
	my %namelist;
	
	my @hkeys=keys(%hash);
	my $hkeys=@hkeys;
	print "\ntotal  $hkeys seq ";
	for(my $i=0;$i<$hkeys;$i++){
		$tpv=$hash{$hkeys[$i]};
		#print "$tpv ";
		if($n{$tpv}>0){
			$n{$tpv}++;
			$namelist{$tpv}.=" $hkeys[$i]";
		}else{
			
			$n{$tpv}=1;
			$namelist{$tpv}="$hkeys[$i]";
		}		
	}
	
	my @nkeys=keys(%n);
	my $nkeys=@nkeys;	
	print "classified to  $nkeys families\n";
	for(my $i=0;$i<$nkeys;$i++){
		print OUT "$n{$nkeys[$i]}\t$namelist{$nkeys[$i]}\n"; 
	}
	close(OUTF);
}



sub CombinCDHitClstrToClass{
	# * is used one
	#>Cluster 2
	#0	4769nt, >J7L_2109_3... at +/98.70%
	#1	4775nt, >J7L_2109_5... *
	my $name="NA";
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);	
		$tp=@tp;
		if($lineIN=~/>Cluster/){	
			if($name ne "NA"){
				print OUT "\n$count\t$name";
			}
			$name="";
			$count=0;
		}else{
			if($lineIN=~/>/){
				my $tpname=substr($tp[2],1,(index($tp[2],0,'...')-2));
				my $tpname=$tp[2];
				if($tp[3] eq "*"){
					$name=$tpname."\t".$name;
				}else{
					$name=$name."\t".$tpname;
				}
				$count++;
			}
		}
	}
	print OUT "$count\t$name";
	print "\n";
}

sub count01Matrix{
	#input value matrix
	#first row line is title
	# type	b	c	d
	# a1	0	1	0
	# a2	1	1	0
	my @row;
	my $crow=0;
	my @line;
	my $cline=0;
	my @matrix;
	
	#first row
	$lineIN=<IN>;
	my @tp= &Getheadinfo($lineIN);	
	my $tp=@tp;	
	for(my $i=1;$i<$tp;$i++){
		$row[$i-1]-$tp[$i];
	}
	$crow=$tp-1;
	print OUT "$lineIN\n";
	print "$lineIN\nType\t";
	#read in lines
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);	
		$tp=@tp;
		
		my $tpl=$cline;
		for(my $j=0;$j<$cline;$j++){
			if($tp[0] eq $line[$j]){
				$tpl=$j;
				#print "1";
				last;
			}	
		}
		if($tpl==$cline){
			#print "2";
			$line[$tpl]=$tp[0];
			print "$tp[0]\t";
			for(my $i=1;$i<$tp;$i++){
				$matrix[$tpl][$i-1]=$tp[$i];
			}	
			$cline++;
		}else{
			for(my $i=1;$i<$tp;$i++){
				$matrix[$tpl][$i-1]=$matrix[$tpl][$i-1]+$tp[$i];
			}			
		}
	}
	for (my $i=0; $i<$cline; $i++){
		print OUT "$line[$i]";
		for (my $j=0; $j<$crow; $j++){
			print OUT "\t$matrix[$i][$j]";
		}		
		print OUT "\n";
	}
	print "\n";
}

sub ChechPatternDiffer{
	#input pattern1 pattern2 
	# use pattern check pattern2
	#out pattren1 pattern2 pattern_differ pattern2_1_differ pattern2_change 
	while($lineIN=<IN>){
		chomp($lineIN);	
		print OUT "$lineIN\t";
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		my $len=length($tp[0]);
		my $dif="";
		my $dif2-"";
		my $cha="";
		for(my $i=0; $i<$len; $i++){
			my $p1=substr($tp[0],$i,1);
			my $p2=substr($tp[1],$i,1);
			if($p1 eq $p2){
				$dif.="0";
				$dif2.="0";
				$cha.="$p2";
			}else{
				$dif.="1";
				if($p2 eq "1"){
					$dif2.="2";
					$cha.="$p2";					
				}else{
					$dif2.="3";
					$cha.="$p1";		
				}
			}
		}
		print OUT "$dif\t$dif2\t$cha\n";
	}
}


sub MatrixToList{
	#input value matrix
	my $count=1;
	
	$lineIN=<IN>;
	chomp($lineIN);	
	my @row= &Getheadinfo($lineIN);	
	my $rowno=@row;	
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		
		for(my $i=1;$i<$tp;$i++){
			print OUT "$tp[0]\t$row[$i]\t$tp[$i]\n";			
		}	
	}
	print "down\n";
}


sub MatrixToFrame{
	#input value matrix
	my $count=1;

	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		for(my $i=0;$i<$tp;$i++){
			my $a1=$i+1;
			print OUT "$a1	$count	$tp[$i]\n";
		}	
		$count++;
	}
}


sub TableChangeTo01ByMutiIf{

	#input value matrix
	$lineIN=<IN>;
	chomp($lineIN);	
	my @tp= &Getheadinfo($lineIN);	
	my $tp=@tp;
	print OUT $tp[0];
	for(my $i=1;$i<$tp;$i+=4){
		print OUT "\t$tp[$i]";
		#print "\t$tp[$i]";
	}
	print OUT "\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		print OUT $tp[0];
		for(my $i=1;$i<$tp;$i+=4){
			if($tp[$i+3] eq "NA"){
				print OUT "\t0";
			}
			elsif($tp[$i+3]<=0.05){
				if($tp[$i+1]>=1){
					print OUT "\t1";
				}
				elsif($tp[$i+1]<=-1){
					print OUT "\t2";
				}else{
					print OUT "\t0";
				}			
			}else{
				print OUT "\t0";
			}

		}	
		print OUT "\n";
	}
}


sub TableChangeTo01ByCutoff{
	my $cut=$_[1];
	#input value matrix
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		for(my $i=0;$i<$tp;$i++){
			if($tp[$i]>$cut){
				print OUT "1\t";
			}else{
				print OUT "0\t";
			}
		}	
		print OUT "\n";
	}
}

sub CounttypePercentRegionValue{
	#input must be sorted
	my @cut=(0,5,25,50,75,95,100);
	my $cut=@cut;
	my @value;
	my $type="NA";
	my $count=0;
	print OUT "Type	No";
	print "Type	No";
	for(my $i=0;$i<$cut;$i++){
		print OUT "\tP".$cut[$i];
		print " P".$cut[$i];
	}
	print OUT "\n";
	print "\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if($tp[0] ne $type){
			if($type ne "NA"){
				#count
				#@value=sort{$a <=> $b}@value;
				#print "@value\n";
				print OUT "$type\t$count\t$value[0]";
				for(my $i=1;$i<($cut-1);$i++){
					
					my $o=$count*$cut[$i]/100;	
					my $c=$o;
					my $l=index($o,'.',0);
					if($l>=0){
						if($l==0){
							$c=0;
						}else{
							$c=substr($o,0,$l);
						}
						my $c5=substr($o,($l+1),1);
						if($c5>=5){
							$c++;
						}
					}
					if($c<=0){
						$c=1;
					}
					#print "$c\t";
					print OUT "\t$value[$c-1]";
				}
				print OUT "\t$value[$count-1]\n";
				print ".";
			}
			#add new
			$type = $tp[0];
			$count=0;
			@value=();
		}
		$value[$count]=$tp[1];
		$count++;
	}
	
	#last
	@value=sort{$a <=> $b}@value;
	print OUT "$type\t$count\t$value[0]";
	for(my $i=1;$i<($cut-1);$i++){
		my $o=$count*$cut[$i]/100;
		my $c=$o;
		my $l=index($o,'.',0);
		if($l>0){
			$c=substr($o,0,$l);
			my $c5=substr($o,($l+1),1);
			if($c5>=5){
				$c++;
			}
		}
		print OUT "\t$value[$c-1]";
	}
	print OUT "\t$value[$count-1]";
	print "\n";
}

sub ChechSameFamilybyRow{
	my $r1=$_[1];
	my $r2=$_[2];
	my $cut=$_[3];
	my $dif="";
	my $outf=$ARGV[0].'.out.differ';
	open(OUTF,'>',"$outf");
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);		
		my @tp1=&Getheadinfo($tp[$r1],$cut);
		my @tp2=&Getheadinfo($tp[$r2],$cut);
		if($tp1[0] eq $tp2[0]){
			print OUT $lineIN."\n";
		}else{
			print OUTF $lineIN."\n";
		}
	}
}


sub getBlastRelationship{
	# input gene must sorted
	#TE1 TE2 
	
	my @rel;	 #0-TE1 1-TE2 2-number
	my $rcounr=0;
	my $ccount=0;
	my %TEclass; #save class for each
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);		
		#my $tp=@tp;
		my $canadd=1;
		for(my $i=0;$i<$rcount;$i++){
			#print "$tp[0] eq $rel[$i][0] && $tp[1] eq $rel[$i][1]\n";
			if($tp[0] eq $rel[$i][0] && $tp[1] eq $rel[$i][1]){
				#add number for already have
				$rel[$i][2]++;
				$canadd=-1;
				#print "-";
				last;
			}
			elsif($tp[1] eq $rel[$i][1] && $tp[0] eq $rel[$i][0]){
				#old
				$canadd=-1;
				#print "-";
				last;
			}
		}
		#add brand new
		if($canadd==1){
			#print "\nnew ";
			$rel[$rcount][0]=$tp[0];
			$rel[$rcount][1]=$tp[1];
			$rel[$rcount][2]=1;
			#check class
			my $c0=$TEclass{$rel[$rcount][0]};
			my $c1=$TEclass{$rel[$rcount][1]};
			#print "add $rcount $c0 $c1\n";
			if($c0 ne ""){
				if($c1 ne ""){
					if ($c0 != $c1){
						#combin
						my @TE=keys(%TEclass);
						$TE=@TE;
						for(my $j=0;$j<$TE;$j++){
							if($TEclass{$TE[$j]} == $c1){
								$TEclass{$TE[$j]} == $c0;
							}
						}
						print "c";
					}else{
						#same
						print "t";
					}
				}else{
					$TEclass{$rel[$rcount][1]}=$c0;
					print ".";
				}
			}else{
				if($c1 ne ""){
					$TEclass{$rel[$rcount][0]}=$c1;
					print ".";
				}else{
					#add new
					$TEclass{$rel[$rcount][0]}=$ccount;
					$TEclass{$rel[$rcount][1]}=$ccount;
					$ccount++;
					print "+";
				}
			}
			$rcount++;	
		}
	}
	#output hash
	my @TE=keys(%TEclass);
	my $TE=@TE;
	my @out;   #0-number 1-string
	for(my $i=0;$i<$ccount;$i++){
		$out[$i][0]=0;
	}
	for(my $j=0;$j<$TE;$j++){
		$out[$TEclass{$TE[$j]}][0]++;
		$out[$TEclass{$TE[$j]}][1]=$out[$TEclass{$TE[$j]}][1]."\t".$TE[$j];
	}
	for(my $i=0;$i<$ccount;$i++){
		print OUT "$i\t$out[$i][0]\t$out[$i][1]\n";
	}
	print OUT "\n-----------------\n";
	
	for(my $i=0;$i<$rcount;$i++){
		my $c=$TEclass{$rel[$i][0]};
		print OUT "$rel[$i][0]\t$rel[$i][1]\t$rel[$i][2]\t$c\n";
	}	
	print "\n$rcount $ccount\n";
}


sub RemoveNestingCodeinGene{
	# input gene must sorted
	# gene_ID code1 TE_info other_info
	my $coutoff=1000;
	if($_[1]>0){
		$cutoff=$_[1];
	}
	my $st=5;
	my $ed=6;
	
	my $name="NA";
	my @info;  #2-dimensional array, each :0-ID 1-code 2-st 3-ed 4-info
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		my $addtp="$tp[2]";
		for(my $i=3;$i<$tp;$i++){
			$addtp=$addtp."\t".$tp[$i];
		}	
		#print "$addtp\n";
		my @tptp= &Getheadinfo($tp[2],"|",".");
		if($name ne $tp[0]){
			if($name ne "NA"){
				#out
				print "+";
				for(my $i=0;$i<$count;$i++){
					if($info[$i][0] ne "d"){
						print OUT "$info[$i][0]\t$info[$i][1]\t$info[$i][2]\t$info[$i][3]\t$info[$i][4]\n";
					}					
				}
			}
			#add new
			$name=$tp[0];
			@info=();
			#6 7
			#

			$info[0]=[$tp[0],$tp[1],$tptp[$st],$tptp[$ed],$addtp];
			$count=1;
		}else{
			#print "z $tp[1] $info[$i][1]";
			#check code
			my $add=-1;   #-1 add tail, other add loc
			
			for(my $i=0;$i<$count;$i++){
				#print "bbb $i $info[$i][0] $info[$i][1]\n";
				if($tp[1] eq $info[$i][1]){ 
					#same code check
					#1-2 before corse 1 2-1 before corse 2 3-2 in 1 4-1 in 2 5-same loc 6-no corse 1 before 2 7-no corse 2 before 1
					my @o=&Comparelocal($name,$tptp[$st],$tptp[$ed],$name,$info[$i][2],$info[$i][3]);
					#print "count $count $i  $tptp[$st] $tptp[$ed] $info[$i][2] $info[$i][3] @o\n";
					if($o[0]==1){
						$tptp[$st]=$info[$i][2];
						$info[$i][0]="d";
						$add=$i;
					}
					elsif($o[0]==2){
						$tptp[$ed]=$info[$i][3];
						$info[$i][0]="d";
						$add=$i;						
					}
					elsif($o[0]==3 || $o[0]==5){
						$info[$i][0]="d";
						$add=$i;						
					}
					elsif($o[0]==4){
						$tptp[$st]=$info[$i][2];
						$tptp[$ed]=$info[$i][3];
						$info[$i][0]="d";
						$add=$i;					
					}		
					elsif($o[0]==6 && $o[1]<=$cutoff){
						$tptp[$ed]=$info[$i][3];
						$info[$i][0]="d";
						$add=$i;						
					}
					elsif($o[0]==7 && $o[1]<=$cutoff){
						$tptp[$st]=$info[$i][2];
						$info[$i][0]="d";
						$add=$i;					
					}			
				}
			}
			#print "add=$add\n";
			#add info	
			if($add==-1){			
				$info[$count]=[$tp[0],$tp[1],$tptp[$st],$tptp[$ed],$addtp];
				print ".";
				$count++;				
			}else{
				#combin to last one
				#print "$add\n";
				my $infoadd="$info[$add][4]\t$addtp";
				$info[$add]=[$tp[0],$tp[1],$tptp[$st],$tptp[$ed],$infoadd];
				print "c";
			}			
		}
	}
	for(my $i=0;$i<$count;$i++){
		if($info[$i][0] ne "d"){
			print OUT "$info[$i][0]\t$info[$i][1]\t$info[$i][2]\t$info[$i][3]\t$info[$i][4]\n";
		}					
	}	
	print "\n";	
}



sub GetGeneCodebyTECode{
	# input gene must sorted
	# gene_ID code1 code2 code3 code4 code5 code6
	# gene_ID code1 code2 code3 code4 code5 code6
	
	my $name="NA";
	my @code=(0,0,0,0,0,0);
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
		if($name ne $tp[0]){
			if($name ne "NA"){
				#out
				my $total=$code[0]+$code[1]+$code[2]+$code[3]+$code[4]+$code[5];
				print OUT "$name\t$total\t@code\n";

			}
				#add new
			$name=$tp[0];
			@code=(0,0,0,0,0,0);
		}
		#check code
		for(my $i=1;$i<$tp;$i++){
			if($tp[$i]>$code[$i-1]){
				$code[$i-1]=$tp[$i];
			}
		}			
	}
	my $total=$code[0]+$code[1]+$code[2]+$code[3]+$code[4]+$code[5];
	print OUT "$name\t$total\t@code\n";
}

sub ChangeStringtoCode{
	my @cut=(0.7, 0.6, 0.5);  #<100 >1000
	
	#input
	#jap-r5(resource)
	#jap-r5-r7-w10-w8-w9
	#cut length use the smallest resource 
	my @resname=('jap','w9','w8','r5','w10','r7');
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @code=(0,0,0,0,0,0);
		#get res smallest lenth
		my @res=&Getheadinfo($lineIN,'-');
		my $res=@res;
		#get code
		#print "$res @res\n";
		for(my $i=0;$i<$res;$i++){
			for(my $j=0;$j<6;$j++){
				if($res[$i] eq $resname[$j]){
					$code[$j]=1;
					last;
				}
			}

		}
		#print "$len\n";
		print OUT "$lineIN\t$code[0]$code[1]$code[2]$code[3]$code[4]$code[5]\n";
		#print "$lineIN\t@code\n";
	}

}


sub Cut_to_Max_value
{
	my $loc=$_[1];
	my $max=$_[2];  
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp=&Getheadinfo($lineIN);
		if($tp[$loc] > $max){
			$tp[$loc] = $max;
		}
		print OUT "@tp\n";
	}
}


sub CutLenToCode{
	my @cut=(0.7, 0.6, 0.5);  #<100 >1000
	
	#input
	#len l1 l2 l3 l4 l5 l6 jap-r5(resource)
	#jap-r5-r7-w10-w8-w9
	#cut length use the smallest resource 
	my @resname=('jap','w9','w8','r5','w10','r7');
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
		
		#get res smallest lenth
		my @res=&Getheadinfo($tp[$tp-1],'-');
		my $res=@res;
		#print "@res\n";
		#get res smallest lenth
		my $len=99999;
		for(my $i=0;$i<$res;$i++){
			for(my $j=0;$j<6;$j++){
				if($res[$i] eq $resname[$j]){
					my $tplen=$tp[$j];
					if($len>$tplen){
						$len=$tplen;
					}
					last;
				}
			}

		}
		#print "$len\n";
		if($tp[0]<100){
			$len=$len*$cut[0];
		}
		elsif($tp[0]<1000){
			$len=$len*$cut[1];
		}else{
			$len=$len*$cut[2];
		}
		#print "$len\n";
		my $code="";
		for(my $i=0;$i<$tp-1;$i++){
			if($tp[$i]>=$len){
				$code=$code."1";
			}else{
				$code=$code."0";
			}
		}
		print OUT "$lineIN\t$code\n";
		#print "$lineIN\t$code\n";
	}

}

sub CutSortedDataToMeanValue{
	my $cut=$_[1];
	my $row=$_[2];
	#input muste be sorted
	
	my $count=0;
	my $total=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		print OUT "$lineIN\n";
		my @tp= &Getheadinfo($lineIN);		
		my $tpc=($total+$tp[$row])/($count+1);
		if($tpc<=$cut){
			print "$count\t$total\t+1\t$tpc\n";
			last;
		}	
		$count++;
		$total=$total+$tp[$row];
		#print "$total\n";	
	}
}



sub StasticGetRandlines{
	
	my $n=$_[1];
	my $count=0;
	my @a;
	while($lineIN=<IN>){
		chomp($lineIN);	
		$a[$count]=$lineIN;
		$count++;
	}
	for(my $i=0;$i<$n;$i++){
		$int=int(rand($count));
		print OUT "$a[$int]\n";
	}
}

sub StasticEachRowDivideMean{
	#Each row Divide by Mean
	#all must be value, no title

	#read in data, line to array
	my $count1=0;
	my $count2=0;
	my @array;
	my @tp;
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);		
		$array[$count1]=[@tp];
		$count1++;
	}
	$count2=@tp;
	print "$count1 X $count2\n";
		
	#count
	@out;
	for(my $j=0;$j<$count2;$j++){
		my $mid=0;
		for(my $i=0;$i<$count1;$i++){
			$mid+=$array[$i][$j];
		}
		$mid=$mid/($count1-1);
		print "min=$mid\n";
		for(my $i=0;$i<$count1;$i++){
			$out[$i][$j]=$array[$i][$j]/$mid;
			#print "$array[$i][$j] $out[$i][$j]\n";
		}
	}
	
	#output
	for(my $i=0;$i<$count1;$i++){
		for(my $j=0;$j<$count2-1;$j++){
			print OUT $out[$i][$j]."\t";
		}
		print OUT $out[$i][$count2-1]."\n";
	}
			
}



sub StasticCountAlternativeSplicingTEloc{
	
	#infile: 		0-Gene 1-mRNA No 2-mRNA 3-TEloc
	#outfile: 	0-Gene 1-mRNA No 2-mRNA.1 TEloc 3-mRNA.2 TEloc 3-mRNA.3 TEloc
	#outfile-d:	0-list 1-Gene 2-mRNA.1	3-TEloc 
	print OUT "Gene	AltNo	StateNo	nRMA.1	nRMA.2	nRMA.3	nRMA.4\n";
	my $name="NA";
	my @mRNA;
	my $count=0;
	my $geneno=1;
	my $ln="NA";
	my $outf1=$ARGV[0].'.detail';	
	open(OUTD,'>',"$outf1");		
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($tp[0] ne $name){
			if($name ne "NA"){
				#print $count;
				my $o="";
				my $a=$mRNA[0];
				my $oa=0;
				for(my $i=0;$i<$count;$i++){
					$o= "$o\t$mRNA[$i]";
					my $ii=$i+1;
					print OUTD "$geneno\t$name\t$name.$ii\t$mRNA[$i]\n";
					if($a ne $mRNA[$i]){
						$oa=1;
					}
				}
				print OUT "\t$oa$o\n";
				$geneno++;
			}			
			$name=$tp[0];
			$count=$tp[1];
			print OUT "$name\t$count";
			@mRNA=();
			for(my $i=0;$i<$count;$i++){
				$mRNA[$i]="0";
			}
			
		}

		my $n=substr($tp[2],(index($tp[2],".",)+1),length($tp[2]))-1;
		#print "$n ";
		my $te=0;
		if($lineIN=~/exon|UTR/i){	
			$te=2;
		}else{
			$te=1;
		}
		$mRNA[$n]=$mRNA[$n]>$te?$mRNA[$n]:$te;

	}
	my $o="";
	my $a=$mRNA[0];
	my $oa=0;
	for(my $i=0;$i<$count;$i++){
		$o= "$o\t$mRNA[$i]";
		if($a ne $mRNA[$i]){
			$oa=1;
		}
	}
	print OUT "\t$oa$o\n";
	close(OUTD)

}

sub StasticCountAverageBy0or1{
	#first row V=value S=state, only count S=1 mean
	#out put x=v y=s 
	#output last line is total mean
	my @av;			#data
	my @as;
	my $cv=0;		#count
	my $cs=0;
	my $lv=0;		#number
	my $ls=0;
	my @title;
	my @state;

	#first line
	$lineIN=<IN>;
	chomp($lineIN);	
	@title= &Getheadinfo($lineIN);
	#second line
	$lineIN=<IN>;
	chomp($lineIN);	
	my @state= &Getheadinfo($lineIN);
	my $state=@state;
	for(my $i=0;$i<$state;$i++){
		if($state[$i] eq "V" || $state[$i] eq "v"){			
			$av[$lv][0]=$title[$i];
			$lv++;
		}elsif($state[$i] eq "S" || $state[$i] eq "s"){
			$as[$ls][0]=$title[$i];
			$ls++;
		}
	}
	#read in data
	#print "@state\n";
	my $count=1;
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		$cv=0;
		$cs=0;
		for(my $i=0;$i<$tp;$i++){
			if($state[$i] eq "V" || $state[$i] eq "v"){
				$av[$cv][$count]=$tp[$i];
				$cv++;
			}elsif($state[$i] eq "S" || $state[$i] eq "s"){
				$as[$cs][$count]=$tp[$i];
				$cs++;
			}
		}
		$count++;
	}
	print "value = $lv & state = $ls\n";
	#add total mean
	$as[$ls][0]="mean";
	for (my $i=1;$i<$count;$i++){
		$as[$ls][$i]=1;
		#print "$as[14][$i]";
	}
	$ls++;
	
	#count	
	print OUT "Title";
	for (my $j=0;$j<$ls;$j++){
		print OUT "\t$as[$j][0]";
	}
	print OUT "\n";
	for (my $i=0;$i<$lv;$i++){
		print OUT "$av[$i][0]";
		for (my $j=0;$j<$ls;$j++){			
			my $o=0;
			my $n=0;
			for(my $k=1;$k<$count;$k++){
				if($as[$j][$k]==1){
					if($av[$i][$k] ne "x"){
						$o+=$av[$i][$k];
						#print "$av[$i][$k] ";
						$n++;
					}
				}
			}
			if($n>0){
				$o=$o/$n;
			}else{
				$o="NA";
			}
			print OUT "\t$o";
		}
		print OUT "\n";
	}
}


sub sortRow{
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		@tp=sort(@tp);
		print OUT "@tp\n";
	}

}

sub RemoveSequencWithSpecialCharacters{
	
	my $name="NA";
	my $seq="";
	
	my $canadd=1;
	
	while($lineIN=<IN>){
		chomp($lineIN);			
		if($lineIN=~/>/){			
			if($name ne "NA" && $canadd==1){
				print OUT "$name$seq";
				
			}
			$name=$lineIN."\n";;
			$seq="";
			$canadd=1;
		}else{
			#test
			if($lineIN=~/B|[D-F]|[H-S]|[U-Z]/i){
			#if($lineIN=~/B/i){
				$canadd=0;
				print "$lineIN\n";
			}
			$seq=$seq.$lineIN."\n";;
		}
	}
	if($canadd==1){
		print OUT "$name$seq";				
	}
}

sub RemoveSequenceBy_NNN_MoreThanCutoff{
	my $cut=11;
	if ($_[1]>0){
		$cut=$_[1];
	}
	
	my $s1="";
	my $s2="";
	
	my $name="NA";
	my $seq="";
	
	for(my $i=0;$i<$cut;$i++){
		$s1=$s1.'N';
		$s2=$s2.'n';
	}
	print "Search $s1\t$s2\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);			
		if($lineIN=~/>/){			
			if($name ne "NA"){
				#test NNN
				#print "$seq\n";
				my $a=index($seq,$s1,1);
				my $b=index($seq,$s2,1);
				#print "$a $b\t";
				if($a>0 || $b>0){
					print ".";
				}else{
					print OUT "$name\n$seq\n";		
				}
				
			}
			$name=$lineIN;
			#print "$name\t";
			$seq="";
		}else{
			$seq=$seq.$lineIN;
		}
	}
	my $a=index($seq,$s1,);
	my $b=index($seq,$s2,);
	#print "$a $b\t";
	if($a>0 || $b>0){
		print ".";
	}else{
		print OUT "$name\n$seq\n";		
	}
}


sub Find_NNN_MoreThanCutoff{
	my $cut=11;
	if ($_[1]>0){
		$cut=$_[1];
	}
	
	my $s1="";
	my $s2="";
	
	my $name="NA";
	my $seq="";
	
	for(my $i=0;$i<$cut;$i++){
		$s1=$s1.'N';
		$s2=$s2.'n';
	}
	print "Search $s1\t$s2\n";
	my $outd=$ARGV[0].'.delete';
	open(OUTD,'>',"$outd");	
	while($lineIN=<IN>){
		chomp($lineIN);			
		if($lineIN=~/>/){			
			if($name ne "NA"){
				#test NNN
				#print "$seq\n";
				my $a=index($seq,$s1,);
				my $b=index($seq,$s2,);
				#print "$a $b\t";
				if($a>=0 || $b>=0){
					#print "$a $b\t";
				#if(index($seq,0,$s1)>0 || index($seq,0,$s2)>0){
					print OUT "$name\t$seq\n";
					#print "$name\t$seq\n";				
				}else{
					print OUTD "$name\n$seq\n";
				}
			}
			$name=$lineIN;
			#print "$name\t";
			$seq="";
		}else{
			$seq=$seq.$lineIN;
		}
	}
	my $a=index($seq,$s1,);
	my $b=index($seq,$s2,);
	#print "$a $b\t";
	if($a>0 || $b>0){
		#print "$a $b\t";
	#if(index($seq,0,$s1)>0 || index($seq,0,$s2)>0){
		print OUT "$name\t$seq\n";
		#print "$name\t$seq\n";				
	}else{
		print OUTD "$name\n$seq\n";
	}
	close(OUTD);
}


sub AddNewClass{
	#when blastclustaladd new data, use old class name
			
	my $newsb=$_[1];
	my $outf=$ARGV[0].'.out.error';
	open(OUTF,'>',"$outf");

	
	while($lineIN=<IN>){
		chomp($lineIN);			
		my @tp= &Getheadinfo($lineIN);
				
		my $tp=@tp;
		my @type;
		my $typecount=0;
		my @newarray;
		my $count=0;
		for(my $i=$st; $i<$tp; $i++){
			if(index($tp[$i],$newsb)>=0){
				$newarray[$count]=$tp[$i];
				$count++;
			}else{
				my @tpsub= &Getheadinfo($tp[$i],'_');
				my $add=1;
				for(my $i=0;$i<$typecount;$i++){
					if($type[$i] eq $tpsub[0]){
						$add = -1;
						last;
					}
				}
				
				if($add == 1){
					$type[$typecount]=$tpsub[0];
					$typecount++;
				}
			}
		}

		#print ".";
		if($typecount>0){
			for (my $i=0;$i<$count;$i++){
				print OUT "$newarray[$i]\t@type\n";
			}
		}else{
			for (my $i=0;$i<$count;$i++){
				print OUT "$newarray[$i]\t0\n";
			}
		}
	}
	print "\n";
	close(OUTF);
}



sub TableCombinTo01Matrix{
	my $cutoff=0.5;
	if($_[1]>0){
		$cutoff=$_[1];
	}
	#0-name 1-length 2-3-4-5-6-7-table value
	while($lineIN=<IN>){
		chomp($lineIN);	
		print OUT "$lineIN\t";
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my $cut=$cutoff*$tp[1];
		for(my $i=2; $i<$tp; $i++){
			if($tp[$i] eq "N"){
				print OUT "N";
			}else{
				my $o=$tp[$i]>$cut?1:0;
				print OUT "$o";
			}
		}
		print OUT "\n";
	}
	print "\n";
}


sub countAccumulationValue{
	my $v=0;
	if($_[1]>0){
		$v=$_[1];
	}
	my $a=0;
	while($lineIN=<IN>){
		chomp($lineIN);			
		my @tp= &Getheadinfo($lineIN);
		$a+=$tp[$v];
		print OUT "$lineIN\t$a\n";
	}
}


sub SequenceRemoveGapLocCorrespondence{
	
	my $g='-';	
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($lineIN=~/>/){
			print OUT "$lineIN\n";
		}else{
			for(my $i=1; $i<length($lineIN); $i++){
				if(substr($lineIN,$i,1) ne $g){
					print OUT "$i,";
				}					
			}
			print OUT "\n";
		}
	}
}



sub TransposonsNestCombin{
	
	#$LTR	z_other	RM287-int	13101.t00053	2496	2764
	##DNA	TcMar-Stowaway	rnd-2_family-570	13101.t00053	2535	2764
	#$LTR	z_other	RM275-int	13101.t00053	2764	2905
	#sINE	z_other	rnd-3_family-1321	13101.t00053	2767	2904 165
	#0-type	1-subtype 2-family 3-chrom 4-st 5-ed 6-len
	#if 3 noly contain -int, 
	#
	
	my @array;  #6-len
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);
		$tp=@tp;
		my $o=&FindOverlap($array[4], $array[5], $tp[4], $tp[5]);
		if($o>0 && $array[3] eq $tp[3]){
			my $score_a=0;
			my $score_t=0;
			if(index($array[0],"LTR",0)>=0){
				if (index($array[2],"-LTR",0)>0){
					$score_a+=1;
				}else{
					$score_a-=1;										
				}
			}elsif(index($array[0],"Unknown",0)>=0){
				$score_a-=2;
			}
			else{
				$score_a+=2;
			}
			if(index($tp[0],"LTR",0)>=0){
				if (index($tp[2],"-LTR",0)>0){
					$score_t+=1;
					#print "1";
				}else{
					$score_t-=1;
					#print "2";
				}
			}elsif(index($tp[0],"Unknown",0)>=0){
				$score_t-=2;
				#print "3";
			}else{
				$score_t+=2;
				#print "4";
			}
			if($array[6] > $tp[6]){
				$score_a+=1;
			}else{
				$score_t+=1;
			}
			if($array[6] > ($tp[6]*2)){
				$score_a+=1;
			}
			if($tp[6] > ($array[6]*2)){
				$score_t+=1;
			}
			if($score_a == $score_t){
				if($array[6] > $tp[6]){
					$score_a++;
				}else{
					$score_t++;
				}
			}
			if($score_a > $score_t){
				print OUT "-@tp\t$score_t\n";
				$array[0]="+".$array[0];
				@array=(@array,$score_a);
			}else{
				print OUT "-@array\t$score_a\n";
				@array=(@tp,$score_t);
				$array[0]="+".$array[0];				
			}
		}else{
			print OUT "@array\n";
			@array=@tp;
		}
	}
}


sub RemoveSameNameFas{
	
	my %name;
	my $add=1;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		if($lineIN=~/>/){
			if($name{$tp[0]} ==1){
				$add=0;
				print $lineIN."\n";
			}else{
				$name{$tp[0]} =1;
				$add=1;
			}
		}
		if($add==1){
			print OUT $lineIN."\n";
		}
		
	}
}
sub RemoveLineWtihAllSameRow{
	
	$lineIN=<IN>;
	chomp($lineIN);	
	print OUT "$lineIN\n";
	my $no=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		for(my $i=2;$i<$tp;$i++){
			if($tp[$i-1] ne $tp[$i]){
				print OUT "$lineIN\n";
				last;
			}		
		}
	}
	print "\n";
}


sub RemoveLinesBytotalValue{
	$lineIN=<IN>;
	chomp($lineIN);	
	print OUT "$lineIN\n";
	my $cut =$_[1];
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN,",");	
		my $tp=@tp;
		my $tpvalue=0;
		for(my $i=1;$i<$tp;$i++){
			$tpvalue+=$tp[$i];
		}
		if($tpvalue >=$cut ){
			print OUT "$lineIN\n";
		}
	}
	print "\n";
}


sub RemoveLinesAll0{
	$lineIN=<IN>;
	chomp($lineIN);	
	print OUT "$lineIN\n";
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN,",");	
		my $tp=@tp;
		for(my $i=1;$i<$tp;$i++){
			if($tp[$i] != 0){
				print OUT "$lineIN\n";
				last;				
			}	
		}
	}
	print "\n";
}

sub RemvoeChar{
	while($lineIN=<IN>){
		chomp($lineIN);	
		for(my $i=0;$i<length($lineIN);$i++){
			my $char=substr($lineIN,$i,1);
			if($char ne $_[1]){
				print OUT $char;
				#print $char;
			}
		}
		print OUT "\n";
	}
}


sub FindlinesOverlap{
	
	#loc line
	my $stl = 2;
	#back search number
	my $bac = 10;
	
	my $tpout=$outf.".overlap";
	open(OUTF,'>',"$tpout");
	$lineIN=<IN>;
	chomp($lineIN);	
	my @tp= &Getheadinfo($lineIN,':',".");	
	my $chrom=$tp[$_[1]];
	my $st=$tp[$_[1]+1];
	my $ed=$tp[$_[1]+2];
	my $lastline=$lineIN;
	print OUT "$lineIN\n";
	while($lineIN=<IN>){
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN,':',".");
		my $o=0;		
		if ($chrom eq $tp[$_[1]]){
			$o=&FindOverlap($st, $ed, $tp[$_[1]+1], $tp[$_[1]+2]);
			if($o==1){
				print OUTF "$lastline\n";
				print OUT "$lineIN\n";
				print ".";
			}
		}
		if($out == 0){
			print OUT $lastline."\n";
		}	
		$lastline=$lineIN;	
		$chrom=$tp[$_[1]];
		$st=$tp[$_[1]+1];
		$ed=$tp[$_[1]+2];		
	}
	print OUT "$lineIN\n";
	print "\n";
	close(OUTF);
}



sub Find6alignhit{
	#overlap big than 40% ok
	#overlap <40% list out
	#type	loc	start	end	species	number	
	#out put to:
	#type	loc	start	end	count	JAP-number	NIV-number	GLA-number	BAR-number	GLU-number	MER-number

	my $tpout=$outf.".check";
	open(OUTF,'>',"$tpout");
	my $len=0;
	my @array;	#type	loc	start	end	count	JAP-number	NIV-number	GLA-number	BAR-number	GLU-number	MER-number	JAP-code	NIV-code	GLA-code	BAR-code	GLU-code	MER-code
	my $checkout;
	my $checkadd=-1;   #
	my %sp = ("JAP",0,"NIV",1,"GLA",2,"BAR",3,"GLU",4,"MER",5);
	my $count=1;
	my $count1=1;
	my $max=0;
	
	$lineIN=<IN>;
	my @tp= &Getheadinfo($lineIN);	
	my $spno=5+$sp{$tp[4]};		#which species
	@array=($tp[0],$tp[1],$tp[2],$tp[3],1,0,0,0,0,0,0,0,0,0,0,0,0);
	$array[$spno]=$tp[5];
	$array[$spno+6]=1;
	$max=$tp[5];
	
	while($lineIN=<IN>){
		$count++;
		chomp($lineIN);	
		@tp= &Getheadinfo($lineIN);	
		$spno=5+$sp{$tp[4]};
		#print "$array[2],$array[3],$tp[2],$tp[3]";
		my $o=&OverlapLength($array[2],$array[3],$tp[2],$tp[3]);
		#print " $o ";
		
		if($tp[0] ne $array[0] || $tp[1] ne $array[1]){
			#print "1\n";
			if($checkadd==1){
				$checkadd=-1;			
			}
			#output
			$checkout=$lineIN."\n";
			print OUT "$count1\t$max\t@array\n";
			$count1++;
							
			#add new
			@array=($tp[0],$tp[1],$tp[2],$tp[3],1,0,0,0,0,0,0,0,0,0,0,0,0);
			$array[$spno]=$tp[5];
			$array[$spno+6]=1;	
			$max=$tp[5];					
		}else{				
			if($o==-1){	#output
				#print "2\n";
				if($checkadd==1){
					$checkadd=-1;			
				}
				#output
				$checkout=$lineIN."\n";
				print OUT "$count1\t$max\t@array\n";
				$count1++;
				#add new
				@array=($tp[0],$tp[1],$tp[2],$tp[3],1,0,0,0,0,0,0,0,0,0,0,0,0);
				$array[$spno]=$tp[5];
				$array[$spno+6]=1;
				$max=$tp[5];
				
			}else{				
				#add info to array
				if($array[$spno+6]==1){
					#print "3\n";
					print OUTF "$count\t$count1\t$lineIN\n";
					print OUT "@";
					print ".";
				}else{
					#print "4\n";
				}
				
				$array[2]=$array[2]>$tp[2]?$tp[2]:$array[2];
				$array[3]=$array[3]<$tp[3]?$tp[3]:$array[3];
				$array[$spno]=$tp[5];
				$array[$spno+6]=1;
				$array[4]++;
				$max=$array[5]>$tp[5]?$array[5]:$tp[5];
			}
		}
	}
	print "\n";
	print OUT "$count1\t$max\t@array\n";
	close(OUTF);
}



sub Find6alignbigestRegion{
	########### must have no gap local, else, wrong
	#Find6alignbigestRegion
	#LTR	z_other	RM012-LTR	13112.t04165	fas_r7	4122	4523	4122	4523
	#0-type	1-subtype	2-family 3-chrom 4-species 5-st 6-ed 7-nogap_st 8-nogap_ed
	#loc use largest one, do not combin
	my $tpout=$outf.".overlap";
	open(OUTF,'>',"$tpout");
	my $len=0;	
	my $count=0;
	
	
	$lineIN=<IN>;
	chomp($lineIN);	
	my @array=&Getheadinfo($lineIN);
	$array[9]=$array[8]-$array[7]+1;
	my $add=$array[4];
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		$tp[9]=$tp[8]-$tp[7]+1;
		if($tp[3] ne $array[3]){
			if($count==0){
				$count++;
			}else{
				print OUT "@array\t$add\n";				
			}
			@array=@tp;
			$add=$tp[4];
		}else{		
			my $overlap=&OverlapLength($array[5],$array[6],$tp[5],$tp[6]);
			if($array[9]>$tp[9]){
				$len =$tp[9]*0.5;
			}else{
				$len =$array[9]*0.5;
			}
			
			#print OUT "$len\t$overlap\n";
			if($overlap >=0){
				if($overlap >= $len){
					print ".";
        			if($tp[0] eq $array[0] && $tp[1] eq $array[1]){        				
        				$array[6]=$array[9]>$tp[9]?$array[6]:$tp[6];
						$array[5]=$array[9]>$tp[9]?$array[5]:$tp[5];
						$array[9]=$array[9]>$tp[9]?$array[9]:$tp[9];
						$add=$add."-".$tp[4];
					}else{
						print ".";
						#DNA-TE first, LTR delete
						if($array[0] eq "LTR"){
							#drop @array
							print OUTF "\$"."@array\t$add\n";
							$add=$tp[4];
							@array=@tp;		
						}elsif($tp[0] eq "LTR"){
							#drop @tp
							print OUTF "\$"."@array\t$add\n";
							$add=$tp[4];
							@array=@tp;
						}else{
							#keep the big one
							print OUTF "\$"."@array\t$add\n";		
							if($array[9]<$tp[9]){
								$add=$tp[4];
								@array=@tp;	
							}
										
						}			
					}
				}else{
					
					if($tp[0] eq $array[0] && $tp[1] eq $array[1]){
						print "C";
        				#$array[6]=$array[6]>$tp[6]?$array[6]:$tp[6];
        				$array[6]=$array[9]>$tp[9]?$array[6]:$tp[6];
						$array[5]=$array[9]>$tp[9]?$array[5]:$tp[5];
						$array[9]=$array[9]>$tp[9]?$array[9]:$tp[9];
						$add=$add."-".$tp[4];
					}else{
						print ".";
						print OUT "@array\t$add\n";
						$add=$tp[4];
						@array=@tp;
					}
				}
			}else{
					print OUT "@array\t$add\n";
					$add=$tp[4];
					@array=@tp;				
			}	
		}
	}
	print OUT "@array\t$add\n";
	print "\n";
	close(OUTF);
}


sub SearchAlignRegionGetFas{
	#input
	#No	File	Start	End	Len	Score	ScoreGap	ScoreN
	#1	100.out	3604	3655	51	9.0980	9.0980	10
	my $sp=6;
	my $count=0;
	$lineIN=<IN>;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		my $filen=substr($tp[1],(index($tp[1],'/',)+1),(index($tp[1],'.',)-index($tp[1],'/',)-1));
		my $filename="..//alignments6genomes//".$filen."//mavid.mfa";
		print "$filename\n";
		my @seq=&FasttoArray($filename);
		for(my $i=0;$i<$sp;$i++){
			my $outseq=substr($seq[$i][1],$tp[2]*10,$tp[4]*10);
			my $outname=$tp[0]."_".$filen."_".$seq[$i][0];
			print OUT ">$outname\n";
			print OUT "$outseq\n";
			print ".";
		}
		print OUT "\n";
		print "\n";
	}
	print "\n";
}

sub SearchAlignPercentfromPhyFile{
	my $w=10;  
	#my $cutoff=5;
	#system("mkdir out");
	#read in list
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tpf= &Getheadinfo($lineIN,"/");	
		my $tpf=@tpf;

		#print "$lineIN\t$tpf\t$tpout\n";
		open(INF,'<',"$lineIN");
		#for each file
		my @name;		#seq name
		my @score;		#each score
		my @scoreN;		
		my @scoreGap;	
		my @aln;
		my $countlane=0;
		#my $countConn=0;
		
		my $lineINF=<INF>;
		chomp($lineINF);	
		my @tp= &Getheadinfo($lineINF);	
		my $tp=@tp;
		my $len=$tp[1]; #align length
		my $n=$tp[0];   #align seq number	
		my $lane=$len/(50);
		print OUT "$lineINF\t$len\t$n\t$lane\n";
		print "$len $n $lane\n";
		
		my $times_Gap;
		my $times_N;
		my $times_all;
		while($countlane<$lane){
			for(my $i=0;$i<$n;$i++){
				$lineINF=<INF>;
				chomp($lineINF);	
				@tp= &Getheadinfo($lineINF);	
				$tp=@tp;
				if($countlane==0){
					$name[$i]=$tp[0];
					for(my $j=1;$j<$tp;$j++){
						#$aln[$i]=$tp[$i];
						$times_Gap = (@temp) = ($tp[$j] =~ /(-)/g);					
						$times_N = (@temp) = ($tp[$j] =~ /(N)/g);
						$times_all=$times_Gap+$times_N;
						if($times_N==10){
							$times_N=a;
						}
						if($times_Gap==10){
							$times_Gap=a;
						}
						if($times_all==10){
							$times_all=a;
						}
						$score[$i]=$score[$i].$times_all;
						$scoreN[$i]=$scoreN[$i].$times_N;
						$scoreGap[$i]=$scoreGap[$i].$times_Gap;
					}
				}else{
					for(my $j=0;$j<$tp;$j++){
						#$aln[$i]=$tp[$i];
						my $times_Gap = (@temp) = ($tp[$j] =~ /(-)/g);					
						my $times_N = (@temp) = ($tp[$j] =~ /(N)/g);
						$times_all=$times_Gap+$times_N;
						if($times_N==10){
							$times_N=a;
						}
						if($times_Gap==10){
							$times_Gap=a;
						}
						if($times_all==10){
							$times_all=a;
						}
						$score[$i]=$score[$i].$times_all;
						$scoreN[$i]=$scoreN[$i].$times_N;
						$scoreGap[$i]=$scoreGap[$i].$times_Gap;
					}
				}
			}
			$countlane++;	
			$lineINF=<INF>;	
			print ".";
		}
		close(INF);
		
		my $tpout="out//$tpf[$tpf-2].out";
		open(OUTF,'>',"$tpout");
		print OUTF "$n\t$len\n";
		for(my $i=0;$i<$n;$i++){
			print OUTF "$name[$i]\t$score[$i]\n";
		}
		print OUTF "\n";		
		for(my $i=0;$i<$n;$i++){
			print OUTF "$name[$i]-N\t$scoreN[$i]\n";
		}
		print OUTF "\n";	
		for(my $i=0;$i<$n;$i++){
			print OUTF "$name[$i]-Gap\t$scoreGap[$i]\n";
		}	
		close(OUTF);
		
		print "\n";
	}
	print "\n";
}


sub SearchAlignPefectRegion{
	my $sp=6; #six species
	#cutoff
	my $st_cut=6;				#st score
	my $st_no_cut=3;			#continue loc > $st_cut can start a region
	my $mid_tmp_score_cut=0; 	#mid gap casued low score, nearest 10
	my $total_score_cut=7;		#ok total score
	my $mid_N_cut=1;
			
	#start consider both N & -
	#mid consider only N
	#list: all-N-Gap
	print OUT "No\tFile\tStart\tEnd\tLen\tScore\tScoreGap\tScoreN\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tpf= &Getheadinfo($lineIN,"/");	
		print "$lineIN\n";
		
		#for each file
		my @name;		#seq name
		my @score;		#each score
		my @scoreN;		
		my @scoreGap;
		my $nocount=0;
			
		open(INF,'<',"$lineIN");
		my $len=0;
		my $n=0;
		my $lineINF=<INF>;
		chomp($lineINF);	
		my @tp= &Getheadinfo($lineINF);	
		my $tp=@tp;
		$len=substr($tp[1],0,(length($tp[1])-1))+1; #align length
		$n=$tp[0];   #align seq number	
		print "$n\t$len\n";
		if($n==$sp){		
			for(my $i=0;$i<$n;$i++){
				$lineINF=<INF>;
				chomp($lineINF);	
				@tp= &Getheadinfo($lineINF);	
				$name[$i]=$tp[0];
				$score[$i]=$tp[1];
			}
			$lineINF=<INF>;
			for(my $i=0;$i<$n;$i++){
				$lineINF=<INF>;
				chomp($lineINF);	
				@tp= &Getheadinfo($lineINF);	
				$scoreN[$i]=$tp[1];
			}
			$lineINF=<INF>;
			for(my $i=0;$i<$n;$i++){
				$lineINF=<INF>;
				chomp($lineINF);	
				@tp= &Getheadinfo($lineINF);	
				$scoreGap[$i]=$tp[1];
			}
			close(INF);
			#print "score\t@score\n";	
			#print "N\t@scoreN\n";	
			#print "Gap\t@scoreGap\n";
					
			#analysis
			$len=length($score[0]);
			my @region;   #selected region
			my $count=0;
			
			#place each score
			my @tpregionAll;
			my @tpregionGap;
			my @tpregionN;
			
			#place total score (6 species) in upter
			my @tpScoreAll;
			my @tpScoreGap;
			my @tpScoreN;
			
			#place add new composite score (add one)
			my @new_ScoreN=0;
			my @new_ScoreAll=0;
			my @new_ScoreGap=0;			
	    	
			#place smallest score
			#my @mid_10Score;
			my $min_midScore=0;
			my $min_Score=0;
			my $mid_0_no=0;
				
			my $GapCount=0;
			my $stcount=0;
			my $stloc=0;
			my $canadd=0;  #0-no start  1-start 
			
			my $stop_reson;
						
			for(my $i=0;$i<$len;$i++){	
				#print "\n$i ";
				my $st=10;
				my $gap=10;	
				my $midN=10;
				my $minno=10;
				$min_Score=$len*10;	
				$min_midScore=100;
				
				my $clean=0;
				
				#find smallest score
				for(my $j=0;$j<$n;$j++){
					$tpregionAll[$j][$stcount]=substr($score[$j],$i,1);
					$tpregionAll[$j][$stcount]=$tpregionAll[$j][$stcount] eq "a" ? 0 : (10-$tpregionAll[$j][$stcount]);
					$tpregionGap[$j][$stcount]=substr($scoreGap[$j],$i,1);
					$tpregionGap[$j][$stcount]=$tpregionGap[$j][$stcount] eq "a" ? 0 : (10-$tpregionGap[$j][$stcount]);
					$tpregionN[$j][$stcount]=substr($scoreN[$j],$i,1);
					$tpregionN[$j][$stcount]=$tpregionN[$j][$stcount] eq "a" ? 0 : (10-$tpregionN[$j][$stcount]);
					
					$new_ScoreN[$j]+=$tpregionN[$j][$stcount];
					$new_ScoreAll[$j]+=$tpregionAll[$j][$stcount];
					$new_ScoreGap[$j]+=$tpregionGap[$j][$stcount];	
					
					$st=$st < $tpregionAll[$j][$stcount] ? $st : $tpregionAll[$j][$stcount];
					$gap=$gap < $tpregionGap[$j][$stcount] ? $gap : $tpregionGap[$j][$stcount];
					$midN=$midN < $tpregionN[$j][$stcount] ? $midN : $tpregionN[$j][$stcount];
					#print "$stcount:$new_ScoreN[$j] $st $gap $midN ";
																
					$min_Score = $min_Score < $new_ScoreAll[$j] ? $min_Score : $new_ScoreAll[$j];	
									
					my $mid10=0;
					for(my $k=0;$k<10;$k++){
						if(($stcount-$k)>=0){
							$mid10+=$tpregionAll[$j][$stcount-$k];
							$minno=$k;	
						}else{
							last;
						}
					}			
					$min_midScore = $min_midScore < $mid10 ? $min_midScore : $mid10;		
				}
				if($stcount>=0){
					$min_Score = $min_Score/($stcount+1);
				}	
				if($minno>=0){
					$min_midScore = $min_midScore/($minno+1);
				}	
				#print "$min_Score=$min_midScore ";
						
				#print "$canadd ";
				if($canadd == 0){  #no							
					if($st >= $st_cut){
						#print "|";
						$canadd=1;
						$stloc=$i;
						$stcount++;	
						@tpScoreAll=@new_ScoreAll;
						@tpScoreGap=@new_ScoreGap;
						@tpScoreN=@new_ScoreN;														
					}else{
						$clean=1;	#clean
					}
				}
				elsif($canadd == 1){  #test start	
					#print ".";				
					if($st >= $st_cut){
						$stcount++;	
						@tpScoreAll=@new_ScoreAll;
						@tpScoreGap=@new_ScoreGap;
						@tpScoreN=@new_ScoreN;						
						#print "&";	
						if($stcount>=$st_no_cut){
							$canadd=2;					
						}					
					}else{
						$canadd=0;					
						$clean=1;	#clean
					}
				}
				elsif($canadd == 2){  #start					
					my $add2=0;
					if($midN < $mid_N_cut){  #end
						$add2=1;
						$stop_reson="N";				
					}else{ 
						if($min_midScore<$mid_tmp_score_cut || $min_Score<$total_score_cut){ 	#end
							$add2=1;
        					$stop_reson="G";
						}else{
							if($st==0){								
								$mid_0_no++;
							}else{
								$mid_0_no=0;
							}
							#print "$mid_0_no";
							$stcount++;
							#add new score
							@tpScoreAll=@new_ScoreAll;
							@tpScoreGap=@new_ScoreGap;
							@tpScoreN=@new_ScoreN;												
						}	
					}
					if($add2==1){
						#print "+";
						$nocount++;
						my $tplen=$i-$stloc;
						my $outScoreAll=$tpScoreAll[0];
						my $outScoreGap=$tpScoreGap[0];
						my $outScoreN=$tpScoreN[0];
						#print "=$new_ScoreN[0]=$tpScoreAll[0] ";					
						
						for(my $j=0;$j<$n;$j++){
							$outScoreAll=$outScoreAll < $tpScoreAll[$j] ? $outScoreAll : $tpScoreAll[$j];
							$outScoreGap=$outScoreGap < $tpScoreGap[$j] ? $outScoreGap : $tpScoreGap[$j];
							$outScoreN=$outScoreN < $tpScoreN[$j] ? $outScoreN : $tpScoreN[$j];
						}
						
						#remove all score 0 in tail
						my $i_to=$i-$mid_0_no;	
						$tplen_to=$tplen-$mid_0_no;
						
						#print "$outScoreN $tplen ";
						$outScoreAll=substr($outScoreAll/$tplen_to,0,6);
						$outScoreGap=substr($outScoreGap/$tplen_to,0,6);
						$outScoreN=substr($outScoreN/$tplen,0,6);
						#print "$outScoreN ";
											
						$region[$count]=($stloc,$i_to,$tplen_to,$outScoreAll,$outScoreGap,$outScoreN);
						print OUT "$nocount\t$lineIN\t$stloc\t$i_to\t$tplen_to\t$outScoreAll\t$outScoreGap\t$outScoreN\t$stop_reson\n";
						print "$nocount\t$lineIN\t$stloc\t$i_to\t$tplen_to\t$outScoreAll\t$outScoreGap\t$outScoreN\t$stop_reson\n";
						$count++;
						$canadd=0;
						$clean=1;	#clean
					}					
				}
        	
				#clean
				if($clean==1){
					#print "() ";
					#place each score
					@tpregionAll=();
					@tpregionGap=();
					@tpregionN=();
					
					#place total score (6 species) in upter
					@tpScoreAll=();
					@tpScoreGap=();
					@tpScoreN=();
					
					#place add new composite score (add one)
					@new_ScoreN=();
					@new_ScoreAll=();
					@new_ScoreGap=();			
	    			
					#place smallest score
					#@mid_10Score=();
					$min_midScore=0;
					$min_Score=0;
					$mid_0_no=0;
					
					$GapCount=0;
					$stcount=0;
					$stloc=0;
				}
			}
			print "\n";
		}
	}
}


sub Remove_0_row_percentage{

	my $cut=0.99;
	if($_[1] > 0){
		$cut=$_[1];
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		my $n=$tp*$cut;
		my $count=0;
		for(my $i=1;$i<$tp;$i++){
			if($tp[$i] == 0){
				$count++;
			}
		}
		if($count > $n){
			print "Delete: $tp[0] $count > $n\n";
		}else{
			print OUT "$lineIN\n";
		}
	}
	print "\n";
}


sub RemoveLinefromNumber{
	my $l=@_;
	my %d;
	for(my $i=1;$i<$l;$i++){
		$d{$_[$i]}=1;
		#print "$_[$i]";
	}

	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		my $st=1;
		for(my $i=0;$i<$tp;$i++){

			if($d{$i} ==1){
				#delete
			}else{
				if($st==1){
					$st=0;
				}else{
					print OUT "\t";
				}
				print OUT "$tp[$i]";
			}
		}
		print OUT "\n";
	}
	print "done\n";
}

sub RemoveAllSameLineExceptRows{
	my $l=@_;
	my @d;
	for(my $i=2;$i<$l;$i++){
		$d[$i-2]=$_[$i];
	}
	my $d=@d;
	my @lp;
	my $st=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		my $tp=@tp;
		if($st == 0){			
			@lp=@tp;
			$lineIN=<IN>;
			chomp($lineIN);	
			@tp= &Getheadinfo($lineIN);	
			$st=1;		
		}	
		for(my $i=2;$i<$l;$i++){
			my $cp=0;
			for(my $j=2;$j<$d;$j++){
				if($i==$d[$j]){
					$cp=1;
					last;
				}
			}
			if($cp==0 && $tp[$i] ne $lp[$i]){
				print OUT "@lp\t\n";
				@lp=@tp;
			}
		}
	}
	print "\n";
}

sub FindFirstSymble{
	my $loc=-1;
	my $n=@_;
	for(my $i=1;$i<$n;$i++){
		my $l=index($_[0],$_[$i]);
		if($l>0){
			if($loc==-1){
				$loc=$l;
			}elsif($l<$loc){
				$loc=$l;
			}
		}
	}
	return $loc;
}




sub ChechIndelType{
	#non 0  delete 1  insertion 2 both 3
	#>13102.t04614	DNA#MITE/ORSgTEMT01601924_230_T05261	all=231	r5=3.46=0	r7=196.53=1.73	w8=3.46=0	w9=0=0	w10=3.46=0.43
	#input must be i1 d1 i2 d2 i3 d3 ...
	my $cutoff=60;
	if($_[1]>0){
		$cutoff=$_[1];
	}
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp;
		if($lineIN=~/>/){
			my @atp=&Getheadinfo($lineIN, "=");
			my $atp=@atp;
			my $c=0;
			for(my $i=5; $i<$atp;$i+=3){
				$tp[$c]=$atp[$i];
				$c++;
				$tp[$c]=$atp[$i+1];
				$c++;
			}
		}else{
			@tp= &Getheadinfo($lineIN);
		}
		my $tp=@tp;
		#print "@tp\n";
		my @de;
		my @in;
		my @indel;
		for(my $i=0; $i<$tp;$i+=2){
			#in
			if($tp[$i]>$cutoff){
				$in[$i/2]=1;
				$indel[$i/2]=2;
			}else{
				$in[$i/2]=0;
				$indel[$i/2]=0;			
			}
			#de
			if($tp[$i+1]>$cutoff){
				$de[$i/2]=1;
				$indel[$i/2]+=1;
			}else{
				$de[$i/2]=0;
				$indel[$i/2]+=0;					
			}			
		}
		print OUT "@in\t@de\t@indel\n";
	}
}


sub AlignedFasGetPercent{
	my $name="NA";
	my @type;
	my @seq;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);			
		if($lineIN=~/>/){
			my @tp= &Getheadinfo($lineIN);
			my @tpname=&Getheadinfo($tp[0], "_");
			my $tpsy=$tpname[0]."\t".$tp[1];
			if($tpsy ne $name){				
				if($name ne "NA"){
					my @a=split(/-/,$seq[0]);
					my $a=join("",@a);
					my $l=length($a);
					print OUT "$name\t$type[0]=$l";
					for(my $i=1;$i<$count;$i++){
						my @r=&Aligncompair($seq[0],$seq[$i]);
						print OUT "\t$type[$i]=$r[0]=$r[1]";
					}
					print OUT "\n";
				}
				$name = $tpname[0]."\t".$tp[1];
				$count=0;
				@type=();
				@seq=();
				print "\n$name ";
			}
			$type[$count]=$tpname[1];
			#print "$type[$count] ";
		}else{
			$seq[$count]=$lineIN;
			$count++;
		}
	}
	my @a=split(/-/,$seq[0]);
	my $a=join("",@a);
	my $l=length($a);
	print OUT "$name\t$type[0]=$l";
	for(my $i=1;$i<$count;$i++){
		my @r=&Aligncompair($seq[0],$seq[$i]);
		print OUT "\t$type[$i]=$r[0]=$r[1]";
	}
	print OUT "\n";
}

sub Aligncompair{  
	my $de=0;
	my $in=0;
	my $l=0;
	for(my $i=0;$i<length($_[0]);$i++){
		$a1=substr($_[0],$i,1);
		$a2=substr($_[1],$i,1);
		if($a1 eq "-" && $a2 ne "-"){
			$in++;
		}
		if($a1 ne "-" && $a2 eq "-"){
			$de++;
		}
		if($a1 ne "-"){
			$l++;
		}
	}
	my @r;
	#print "$l $in $de\n";
	if($l>0){
		$r[0]=$in*100/$l;
		$r[1]=$de*100/$l;
		my $d0=index($r[0],".");
		my $d1=index($r[1],".");
		if($d0>0){
			$r[0]=substr($r[0],0,($d0+3));	
		}
		if($d1>0){
			$r[1]=substr($r[1],0,($d1+3));	
		}	
	}
	return @r;

}

sub AlignedFasFilter{
	my $cut=0.5;
	if($_[1]>0){
		$cut=$_[1];
	}
	my $name="NA";
	my $info="";
	my $len=0;
	my $canadd==1;
	my $new=1;
	my $count=0;
	my $dinfo;
	my $outd=$ARGV[0].'.delete';
	open(OUTD,'>',"$outd");
	while($lineIN=<IN>){
		chomp($lineIN);			
		if($lineIN=~/>/){
			my @tp= &Getheadinfo($lineIN, "_");
			if($tp[0] ne $name){				
				if($canadd==1 && $name ne "NA"){
					print OUT "$info";
				}else{
					print OUTD "delete $name $dinfo\n";
					print "delete $name $dinfo\n";
					$count++;
				}				
				$name = $tp[0];
				$info="";
				$canadd=1;
				$len=0;
				$new=1;
				#print " $name";
			}else{
				$new=0;
			}
		}else{
			my @seq=split(/-/,$lineIN);
			my $seq=join("",@seq);
			if($new==1){
				$len=length($seq)*$cut;
				$new=0;
			}else{
				my $tplen=length($seq);
				if($len>$tplen){
					$canadd=0;
					$dinfo="$len $tplen";
				}else{
					#print "+";
				}	
			}		
		}
		$info=$info.$lineIN."\n";	
	}	
	if($canadd==1){
		print OUT "$info";
	}else{
		print OUTD "delete $name $dinfo\n";
		print "delete $name $dinfo\n";
	}
	print "$count seq deleted\n";
	close(OUTD);
}

sub RemoveInfoAfterSymble{
	my $symble="_";
	if($_[1]){
		$symble=$_[1]; 
	}
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		for(my $i=0; $i<$tp;$i++){
			if($tp[$i]=~/$symble/){
				$tp[$i]=substr($tp[$i],0,index($tp[$i],$symble,));	
			}
			print OUT $tp[$i];
			
			if($i== $tp-1){
				print OUT "\n";
			}else{
				print OUT "\t";
			}
		}
	}
}



sub AddLogforPl{
	my @key=("system");
	my $key=@key;
	while($lineIN=<IN>){
		chomp($lineIN);	
		for(my $i=0;$i<$key;$i++){
			if($lineIN=~/$key[$i]/){
				my $o=substr($lineIN,(index($lineIN,"\"",)+1),length($lineIN));
				my $a=index($o,")",)-1;
				$o=substr($o,0,$a);
				print OUT "print \"do-$o\\n\";\n";
				last;
			}
		}
		print OUT $lineIN."\n";
	}
}



sub SeperateGOannoResult{
	#input	
	#Level	GO_ID	GO_Name	GO_Type	Parents_(ACC)	Parents_(Name)	Nodescore	#Seqs	Sequence_Names

	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= split(/\t/,$lineIN);
		my $info="$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]\n";
		
		my @seq= split(/,/,$tp[8]);
		my $seq=@seq;
		for(my $i=0;$i<$seq;$i++){
			print OUT "$seq[$i]\t$info";
		}
	}
}


sub CombinTypetoLines{
	my $type="NA";
	my $symble=" ";
	if($_[1]){
		$symble=$_[1]; 
	}
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= split(/\t/,$lineIN);
		if($tp[0] ne $type){
			if($type ne "NA"){
				print OUT "\n";
			}
			$type=$tp[0];
			print OUT "$type\t";
		}
		print OUT "$tp[1]$symble";
	}
}


sub SignlocalLargeThanMedian{
	my $l=@_;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my $c=0;
		for(my $i=0;$i<$tp;$i++){
			if($i != $_[1]){
				if($_[2] eq "gt"){
					if($tp[$i] > $_[3]*$tp[$_[1]]){
						print OUT "1";
						$c++;
					}else{
						print OUT "0";
					}
				}
				if($_[2] eq "lt"){
					if($tp[$i] < $_[3]*$tp[$_[1]]){
						print OUT "1";
						$c++;
					}else{
						print OUT "0";
					}
				}
			}
		}
		print OUT "\t$c\n";
	}
}


sub RemoveAllSameRowbyType{
	open(OUTD,'>',"$outf.delete");
	my $l=@_;
	my @info;  
	my $type="NA";
	my $count=0;
	my $len=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($type ne $tp[$_[1]]){
			if($type ne "NA"){
				#check
				my $add=0;
				for(my $j=1;$j<$len;$j++){
					for(my $i=2;$i<$l;$i++){
						if ($info[$j][$_[$i]+1] ne $info[0][$_[$i]+1]){
							$add=1;
							#print "$info[$j][$_[$i]] $info[0][$_[$i]]\n";
							last;	
						}
					}
					if($add==1){
						last;
					}
				}
				if($add==1){
					for(my $j=1;$j<$len;$j++){
						for(my $i=1;$i<=$info[$j][0];$i++){
							print OUT "$info[$j][$i]\t";
						}
						print OUT "\n";
					}
				}else{
					$count++;
					print "Delete $type\n";
					print OUTD "$type\n";
				}								
			}
			$type=$tp[$_[1]];
			$len=0;
			@info=();
		}
		@tp=($tp,@tp);
		$info[$len]=[@tp];
		$len++;
	}
	close(OUTD);
}


sub RemoveElementSameRow{
	open(OUTD,'>',"$outf.delete");
	my $l=@_;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my $d=1;
		for(my $i=2;$i<$l;$i++){
			if ($tp[$_[1]] ne $tp[$_[$i]]){
				$d=0;
				print OUT $lineIN."\n";
				last;	
			}
		}
		if($d==1){
			print OUTD $lineIN."\n";
		}
	}
	close(OUTD);
}


sub MatrixRowLineChange{
	my @a;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		for(my $i=0;$i<$tp;$i++){
			$a[$i][$count]=$tp[$i];
		}
		$count++;
	}
	$a=@a;
	print "$a $count\n";
	for(my $i=0;$i<$a;$i++){
		for(my $j=0;$j<$count;$j++){
			print OUT $a[$i][$j]."\t";
		}
		print OUT $a[$i][$count]."\n";
	}
}



sub Seperate_or_Type{
	my $loc=$_[1];
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my @type=&Getheadinfo($tp[$loc],"|");
		my $type=@type;
		#print "$tp $type\n";
		#output
		for(my $i=0;$i<$type;$i++){		
			my $out=$type[$i];
			#for(my $j=0;$j<$tp;$j++){
			#	if($j != $loc){
			#		$out=$out.$tp[$j]."\t";
			#	}else{
			#		$out=$out.$type[$i]."\t";
			#	}
			#}
			print OUT "$lineIN\t$out\n";
		}
	}
}



sub SeperateLinebyType{
	my $loc=0;
	if($_[1]>0){
		$loc=$_[1];
	}
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		my $type=$tp[$loc];
		for(my $i=0;$i<$tp;$i++){
			if($i != $loc){
				print OUT "$type\t$tp[$i]\n";
			}
		}
	}
}


sub MummerResultCombin{
	#input need sort??
	#[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [FRM]   [TAGS]	[querry]	[target]
	#133680	141590	1	7911	7911	7911	99.96	266237	16395165	2.97	0.05	1	-1	scaffold3492	Oglum_chr3s
	
	#0-name-chrom-direct 1-q-start 2-qend 3-t-start 4-t-end 5-q-length 6-number 7-q-map-length 8-q-percent 
	my @scaf;
	my $remove="";
	my $coun=0; 
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $a="$tp[13]\t$tp[14]\t$tp[12]";
		my $canadd=1;
		if($tp[12] eq "-1"){
			my $t=$tp[2];
			$tp[2]=$tp[3];
			$tp[3]=$t;	
		}
		for(my $i=0;$i<$count;$i++){
			if($scaf[$i][0] eq $a){
				$canadd=0;
				#if near
				if($tp[12] eq "1"){
					my $ta=$tp[1]-$scaf[$i][2];  
					my $qu=$tp[3]-$scaf[$i][4]; 
					if($ta>0 && $qu>0 && $ta<$scaf[$i][7] && $qu<$scaf[$i][7]){
						$scaf[$i][2]=$tp[1];
						$scaf[$i][4]=$tp[3];
						$scaf[$i][7]=$scaf[$i][4]-$scaf[$i][3];
						$scaf[$i][8]=$scaf[$i][7]*100/$tp[7];
						$scaf[$i][6]++;
					}else{
						print ".";
						if($scaf[$count][5]>$tp[7]){
							$scaf[$i][1]=$tp[0];
							$scaf[$i][2]=$tp[1];
							$scaf[$i][3]=$tp[2];
							$scaf[$i][4]=$tp[3];
							$scaf[$i][5]=$tp[7];
							$scaf[$i][7]=$tp[4];
							$scaf[$i][8]=$tp[9];						
						}else{
							$remove=$remove.$lineIN."\n";
						}
					}	
				}else{
					my $ta=$tp[0]-$scaf[$i][1];
					my $qu=$tp[3]-$scaf[$i][4];
					if($ta<0 && $qu>0 && $ta<$scaf[$i][7] && $qu<$scaf[$i][7]){
						$scaf[$i][1]=$tp[0];
						$scaf[$i][4]=$tp[3];
						$scaf[$i][7]=$scaf[$i][4]-$scaf[$i][3];
						$scaf[$i][8]=$scaf[$i][7]*100/$tp[7];
						$scaf[$i][6]++;
					}else{
						print ".";
						if($scaf[$count][5]>$tp[7]){
							$scaf[$i][1]=$tp[0];
							$scaf[$i][2]=$tp[1];
							$scaf[$i][3]=$tp[2];
							$scaf[$i][4]=$tp[3];
							$scaf[$i][5]=$tp[7];
							$scaf[$i][7]=$tp[4];
							$scaf[$i][8]=$tp[9];						
						}else{
							$remove=$remove.$lineIN."\n";
						}
					}					
				}
			}	
		}
		#add new
		if($canadd==1){
			$scaf[$count][0]=$a;
			$scaf[$count][1]=$tp[0];
			$scaf[$count][2]=$tp[1];
			$scaf[$count][3]=$tp[2];
			$scaf[$count][4]=$tp[3];
			$scaf[$count][5]=$tp[7];
			$scaf[$count][7]=$tp[4];
			$scaf[$count][8]=$tp[9];
			$scaf[$count][6]=1;
			$count++;
		}
	}
	for(my $i=0;$i<$count;$i++){
		print OUT "$i\t$scaf[$i][0]";
		for(my $j=1;$j<8;$j++){
			print OUT "\t".$scaf[$i][$j];
		}
		my $l=length($scaf[$i][8]);
		if(index($scaf[$i][8],".",)>0){
			$l=index($scaf[$i][8],".",)+3;
		}
		my $o=substr($scaf[$i][8],0,$l);
		print OUT "\t$o\n";
	}
	print OUT "-----------------------------------------\n$remove";
	print "\n";
}


sub LocalCombinbyDirectandDistancd{
	#input
	#DNA	En-Spm	rnd-2_family-263 chr01_13101	11009664	11010094	430	460	871	872	49.3119266055046	+
	#0-type	1-subtype  2-name	3-chrom		4-st		5-ed	    6-hitlen	7-qst	8-qed   9-qlen	    10-percent	     11-direct
	
	my @last;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp=&Getheadinfo($lineIN);
		my $add=1;
		if($tp[2] eq $last[2] && $tp[3] eq $last[3] && $tp[11] eq $last[11]){
			my $dChrom=$tp[4]-$last[5];
			if($tp[11] eq "+"){
				my $dTE=$tp[7]-$last[8];				
				if($dTE>-20 && $dChrom<$dTE+100){
					$add=0;
					#combin
					$last[5]=$tp[5];
					$last[8]=$tp[8];
					$last[6]=$tp[5]-$tp[4];
					$last[10]=$last[6]/$last[9];
					print "+";
				}
			}else{
				my $dTE=$last[7]-$tp[8];
				if($dTE>-20 && $dChrom<$dTE+100){
					$add=0;
					$last[5]=$tp[5];
					$last[7]=$tp[7];
					$last[6]=$tp[5]-$tp[4];
					$last[10]=$last[6]/$last[9];
					print "-";
				}
				
			}
			
			my $dis=$tp[4]-$last[5];
			
		}
		if($add==1){
			if($count>0){
				for(my $i=0;$i<11;$i++){
					print OUT $last[$i]."\t";
				}
				print OUT $last[11]."\n";
			}
			@last=@tp;
		}
		$count++;		
	}
	if($count>0){
		for(my $i=0;$i<11;$i++){
			print OUT $last[$i]."\t";
		}
		print OUT $last[11]."\n";
	}
	print "\n";
}

sub ChangeLTRinfoTolocinfo{
	#in
	#RL346	Chr3	+	7834284	7834488	7838849	7839051	4361
	
	my $outf3=$ARGV[0].'.LTR3';
	open(OUT3,'>',"$outf3");
	my $outf5=$ARGV[0].'.LTR5';
	open(OUT5,'>',"$outf5");
	my $count=0;
	while($lineIN=<IN>){
		$count++;
		chomp($lineIN);	
		my @tp=&Getheadinfo($lineIN);
		print OUT5 ">$tp[0]-5-$count-$tp[7]|$tp[1]:";
		print OUT3 ">$tp[0]-3-$count-$tp[7]|$tp[1]:";
		if($tp[2] eq "+"){
			print OUT5 "$tp[3]..$tp[4]\n";
			print OUT3 "$tp[5]..$tp[6]\n";
		}else{
			print OUT5 "complement($tp[3]..$tp[4])\n";
			print OUT3 "complement($tp[5]..$tp[6])\n";
		}
	}
	close(OUT5);
	close(OUT3);
}


sub ReadinfilebyList{
	#read in list
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		#get aln
		open(TPIN,'<',"$lineIN");	
		my $tplineIN;
		while($tplineIN=<TPIN>){
			chomp($tplineIN);	
			print OUT "$lineIN\t$tplineIN\n";
		}
		$count++;
	}
	print "$count file readed\n"
}




sub ReadCLUSTALfilebyList{
	#read in list
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		#get aln
		my $fas=&ReadAln($lineIN);
		print OUT $fas;
		$count++;
	}
	print "$count aln readed\n"
}


sub RemoveExtremeValuebyType{
	#input mast be: Name Type Value
	#use two times mean +- sd
	#read in
	my $rmtimes=1;
	if($_[1]>0){
		$rmtimes=$_[1];
	}
	open(OUTD,'>',"$outf.delete");
	print "Repeat $rmtimes times\n";
	my @d;
	my @t1;
	my $count=0;
	my $name="NA";
	while($lineIN=<IN>){
		my @tp=&Getheadinfo($lineIN);	
		if($name ne $tp[1]){
			if($name ne "NA"){
				for(my $i=0;$i<$rmtimes-1;$i++){
					@t1=&RemovebyMeanandSD(@t1);					
				}
				my $t1=@t1;
				my $t1d=\@t1;
				my @r=&GetMeanandStandardDeviation($t1d, $t1);
				for(my $i=0;$i<$t1;$i++){
					if($d[$i][2]<($r[0] + $r[1]) && $d[$i][2]>($r[0] - $r[1])){
						print OUT "$d[$i][0]\t$d[$i][1]\t$d[$i][2]\n"
					}else{
						print OUTD "$d[$i][0]\t$d[$i][1]\t$d[$i][2]\n";
					}
				}					
			}
			$name=$tp[1];
			$count=0;
			@d=();
			@t1=();
		}
		chomp($lineIN);	
		my @tp=&Getheadinfo($lineIN);
		$t1[$count]=@tp[2];
		$d[$count]=[@tp];
		$count++;
	}
	for(my $i=0;$i<$rmtimes;$i++){
		@t1=&RemovebyMeanandSD(@t1);
	}
	my $t1=@t1;
	my $t1d=\@t1;
	my @r=&GetMeanandStandardDeviation($t1d, $t1);
	for(my $i=0;$i<$t1;$i++){
		if($d[$i][2]<($r[0] + $r[1]) && $d[$i][2]>($r[0] - $r[1])){
			print OUT "$d[$i][0]\t$d[$i][1]\t$d[$i][2]\n"
		}else{
			print OUTD "$d[$i][0]\t$d[$i][1]\t$d[$i][2]\n";
		}
	}
	close(OUTD);
}

sub RemovebyMeanandSD{
	#input array(name type value), out put delated array
		my $data = \@_;
		my $l=@_;
		my @r=&GetMeanandStandardDeviation($data, $l);
    	#first	
		my @o;
		for(my $i=0;$i<$l;$i++){
			if($_[$i]<($r[0] + $r[1]) && $_[$i]>($r[0] - $r[1])){
				@o=(@o, $_[$i]);
			}
		}	
		return @o;
}
	
sub GetMeanandStandardDeviation{	
	#input $data out put 0-mena 1-standard deviation
    my $mean     = gsl_stats_mean($_[0], 1, $_[1]);
	my $stats_sd_m=gsl_stats_sd_m($_[0], 1, $_[1], $mean); 
	
    print qq{ 
    Load Dataset Number : $_[1]
    Sample mean           $mean 
    standard deviation	  $stats_sd_m
    };
    
    my @r;
    $r[0]=$mean;
    $r[1]=$stats_sd_m;
	return @r;
}



sub SelectbyLargestSubtypeValue{
	#must sorted
	#subtype combined
	#input must be: type subtype value
	print "Cutoff = $_[1] %\n";
	my $type="NA";
	my @subtype;
	my $count=0;
	my $un;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp=&Getheadinfo($lineIN);	
		if($tp[0] ne $type){
			if($type ne "NA"){
				$un=1;
				if($count == 1){
					print OUT "$subtype[0][0]\t$subtype[0][1]\t$subtype[0][2]\n";
					#print ".";
					$un=0;
				}
				elsif($count >1){
					for(my $i=0;$i<$count;$i++){
						my $canadd=1;
						for(my $j=0;$j<$count;$j++){
							if($i != $j){

								if($subtype[$i][2]-$subtype[$j][2] < $_[1]){
									$canadd=0;
									last;
								}
							}
						}
						if($canadd==1){
							print OUT "$subtype[$i][0]\t$subtype[$i][1]\t$subtype[$i][2]\n";
							$un=0;
							#print "0";
							last;
						}
					}
				}
				if($un==1){
					#print "1";
					print OUT "$subtype[0][0]\tUnmap\t$count\n";
				}
			}
			@subtype=();
			$count=0;
			$type=$tp[0];
		}
		$subtype[$count]=[@tp];
		#print $subtype[$count][0];		
		$count++;
	}
	$un=1;
	if($count == 1){
		print OUT "$subtype[0][0]\t$subtype[0][1]\t$subtype[0][2]\n";
		$un=0;
	}
	elsif($count >1){
		for(my $i=0;$i<1;$i++){
			my $canadd=1;
			for(my $j=0;$j<1;$j++){
				if($i != $j){
					if($subtype[$i][2]-$subtype[$j][2] < $_[1]){
						$canadd=0;
						last;
					}
				}
			}
			if($canadd==1){
				print OUT "$subtype[$i][0]\t$subtype[$i][1]\t$subtype[$i][2]\n";
				$un=0;
				last;
			}
		}
	}
	if($un==1){
		print OUT "$subtype[0][0]\tUnmap\t$count\n";
	}					
}


sub SelectLargestValueLine{
	my @type;
	my @lines;
	my @value;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp=&Getheadinfo($lineIN);	
		my $canadd=1;
		$count=@value;
		for(my $i=0;$i<$count;$i++){
			if($tp[$_[1]] eq $type[$i]){
				if($tp[$_[2]] > $value[$i]){
					 splice(@type,$i,1);
					 splice(@lines,$i,1);
					 splice(@value,$i,1);				 
				}else{
					$canadd=0;	
				}
				last;
			}
		}
		if($canadd==1){
			$count=@value;
			$type[$count]=$tp[$_[1]];
			$lines[$count]=$lineIN;
			$value[$count]=$tp[$_[2]];
		}		
	}
	$count=@value;
	for(my $i=0;$i<$count;$i++){
		print OUT "$lines[$i]\n";
	}						
}



sub CombinlSameType{
	my $row=0;
	if($_[1]>0){
		$row=$_[1];
	}
	my $count=0;
	my $subtype="";
	my $type='NA';
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		my $tp=@tp;
		my $tpsub="";
		#print "$tp\n";
		for (my $i=0;$i<$tp;$i++){
			if($i != $row){
				$tpsub=$tpsub."$tp[$i]|";
		#		print "$tpsub\n";
			}
		}
		if ($type ne $tp[$row]){
			if ( $type ne 'NA'){
				#output
				print OUT "$type\t$count=$subtype\n";
			}
			#add new
			$subtype="";
			$type = $tp[$row];
			$count=0;
		}
		$subtype=$subtype."=".$tpsub;
		$count++;			
	}
	print OUT "$type\t$count=$subtype\n";
}


sub CombinltoLargeSubType{
	#only two subtype
	my $count=1;
	my @subtype;
	my @subtypevalue;
	my $type='start';
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		if ($type eq $tp[$_[1]]){
			my $canadd=0;
			for (my $i=0;$i<$count;$i++){
				if ($tp[$_[2]] eq $subtype[$i]){
					$subtypevalue[$i]+=$tp[$_[3]];	
					$canadd=1;	
				}
			}
			if ($canadd==0){
				$subtype[$count]=$tp[$_[2]];
				$subtypevalue[$count]=$tp[$_[3]];
				$count++;
			}
		}else{
			if ( $type ne 'start'){
				if($count>=2){
					my $total=$subtypevalue[0]+$subtypevalue[1];
					if($subtypevalue[0] >= $subtypevalue[1] ){						
						print OUT "$type\t$subtype[0]\t$total\n";
					}else{
						print OUT "$type\t$subtype[1]\t$total\n";
					}
				}else{
					print OUT "$type\t$subtype[0]\t$subtypevalue[0]\n";	
				}	
			}
			$type = $tp[$_[1]];
			$count=1;
			$subtype[0]=$tp[$_[2]];
			$subtypevalue[0]=$tp[$_[3]];
		}			
	}
	if($count>=2){
		my $total=$subtypevalue[0]+$subtypevalue[1];
		if($subtypevalue[0] >= $subtypevalue[1] ){						
			print OUT "$type\t$subtype[0]\t$total\n";
		}else{
			print OUT "$type\t$subtype[1]\t$total\n";
		}
	}else{
		print OUT "$type\t$subtype[0]\t$subtypevalue[0]\n";	
	}						
}


sub MisaResult_samedirect_Seperate{
	#j10004n_all	1	c	(GCC)4gggccgagctagttgaaggaggagagatggcggaaattggggaggaaagatgcaagttgctca(GCT)6gg(GA)15	125	94	218	260	AG
	#seperate to
	#j10004n_all	c	GCC	4	125	94	218	260	AG
	#j10004n_all	c	GCT	6	125	94	218	260	AG
	#j10004n_all	c	GA	15	125	94	218	260	AG

	my $outc=$outf.'.cx';
	open(OUTC,'>',"$outc");
	print "complex overlap saved to: $outc\n";
	
	my $loc=3;
	my $namecount=0;
	my $name=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($name ne $tp[0]){
			$namecount=1;
			$name=$tp[0];
		}
		my $out=$tp[1]."\t".$tp[2]."\t";
		my $tail;
		for(my $i=4;$i<$tp;$i++){
			$tail=$tail."\t".$tp[$i];
		}
        
		my $len=length($tp[3]);
		my @type;
		my @copy;
		my $count=-1;
		my $add=0;
		my @locadd;
		for(my $i=0;$i<$len;$i++){
			my $a=substr($tp[3],$i,1);
			if($a eq "a" || $a eq "t" || $a eq "g" || $a eq "c" || $a eq "n"){
				$add=0;
				$locadd[$count]++;
			}
			elsif($a eq "("){
				$add=1;
				$count++;
				$locadd[$count]=0;
				#print "1";
			}
			elsif($a eq ")"){
				$add =2;	
			}
			elsif($add == 2){
				$copy[$count].=$a;
			}
			elsif($add==1){
				$type[$count].=$a;
			}
        
		}
		my $oa=0;
		my $st=$tp[5];
		my $ed=$tp[6];
		#print OUT "$count ";
		for(my $i=0;$i<=$count;$i++){
			my $o=&SSR_SameDirect_ConsensusSequence($type[$i]);
			my $ol=length($o);
			if($i>0){
				$oa=$oa+($copy[$i-1]*length($type[$i-1]))+$locadd[$i-1];
				$st+=$oa;
			}
			$ed=$st-1+$ol*$copy[$i];
 			if($lineIN=~/\*/){
				print OUTC "$o|$ol|$tp[0]\t$namecount\t$out\t$type[$i]\t$copy[$i]\t$st\t$ed\n";
			}else{       
				print OUT "$o|$ol|$tp[0]\t$namecount\t$out\t$type[$i]\t$copy[$i]\t$st\t$ed\n";
				$namecount++;
			}
		}
		
	}
	close(OUTC);
}


sub MisaResultSeperate{
	#j10004n_all	1	c	(GCC)4gggccgagctagttgaaggaggagagatggcggaaattggggaggaaagatgcaagttgctca(GCT)6gg(GA)15	125	94	218	260	AG
	#seperate to
	#j10004n_all	c	GCC	4	125	94	218	260	AG
	#j10004n_all	c	GCT	6	125	94	218	260	AG
	#j10004n_all	c	GA	15	125	94	218	260	AG

	my $loc=3;
	my $namecount=0;
	my $name=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($name ne $tp[0]){
			$namecount=1;
			$name=$tp[0];
		}
		my $out=$tp[1]."\t".$tp[2]."\t";
		my $tail;
		for(my $i=4;$i<$tp;$i++){
			$tail=$tail."\t".$tp[$i];
		}

		my $len=length($tp[3]);
		my @type;
		my @copy;
		my $count=-1;
		my $add=0;
		my @locadd;
		for(my $i=0;$i<$len;$i++){
			my $a=substr($tp[3],$i,1);
			if($a eq "a" || $a eq "t" || $a eq "g" || $a eq "c" || $a eq "n"){
				$add=0;
				$locadd[$count]++;
			}
			elsif($a eq "("){
				$add=1;
				$count++;
				$locadd[$count]=0;
				#print "1";
			}
			elsif($a eq ")"){
				$add =2;	
			}
			elsif($add == 2){
				$copy[$count].=$a;
			}
			elsif($add==1){
				$type[$count].=$a;
			}

		}
		my $oa=0;
		my $st=$tp[5];
		my $ed=$tp[6];
		#print OUT "$count ";
		for(my $i=0;$i<=$count;$i++){
			my $o=&SSRConsensusSequence($type[$i]);
			my $ol=length($o);
			if($lineIN=~/\*/){
				print OUT "\#";
			}
			if($i>0){
				$oa=$oa+($copy[$i-1]*length($type[$i-1]))+$locadd[$i-1];
				$st+=$oa;
			}
			$ed=$st-1+$ol*$copy[$i];


			print OUT "$o|$ol|$tp[0]\t$namecount\t$out\t$type[$i]\t$copy[$i]\t$st\t$ed\n";
			$namecount++;
		}
	}
}


sub SSRCountConsensusNumber{
	my $l=$_[1];
	#use @dba save array
	#use %db as hash
	$dblevel=$_[1];
	print "input:	$dblevel\n";
	&digui($l-1,4);
	$dba=@dba;

	print "Total:	$dba\n";
	#@dbcount;  #for digui level
	for (my $i=0;$i<$dba;$i++){
		$dba[$i] =~s/1/A/g;
		$dba[$i] =~s/2/T/g;
		$dba[$i] =~s/3/G/g;
		$dba[$i] =~s/4/C/g;
		my $a=&SSRConsensusSequence($dba[$i]);
		if($db{$a}){
			$db{$a}++;
		}else{
			$db{$a}=1;
		}
	}
	my @db = sort{$a <=> $b} keys(%db);
	my $db=@db;
	my $n=$db-2;
	print OUT "$l\t$n\n";
	print OUT "4 need mise 2, 6 need mise 2 and 3, etc\n";
	print "Types:\t$n\n";
	for(my $i=0;$i<$db;$i++){
		print OUT "$db[$i]\t$db{$db[$i]}\n";
	}
}

sub digui{
	my $k=$_[0];
	my $n=$_[1];
    for(my $i=1; $i<=$n; $i++){
        $dbcount[$k] = $i;
        if ($k > 0){
            digui($k - 1, $n);
        }else{
            my $str="";
           
            for (my $j=0; $j < $dblevel; $j++){
                $str=$str.$dbcount[$j];
            }
            @dba=(@dba,$str);
        }
    }
}

sub SSRFindConsensusSequence{
	my $l=0;
	if ($_[1]>0){
		$l=$_[1];
	}
	#print "1\n";
	while($lineIN=<IN>){
		chomp($lineIN);
		my $o1=&SSRSameDirectConsensusSequence($lineIN);
		my $o2=&SSRConsensusSequence($lineIN);
		print OUT "$lineIN\t$o1\t$o2\n";	
		
    }
}


sub SSRConsensusSequence{
	my $in=$_[0].$_[0];
	#print "$in\n";
	my @arr;
	my $len=length($_[0]);
	my $count=0;
	#+
	for(my $i=0;$i<$len;$i++){
		$arr[$count]=substr($in,$i,$len);
		$count++;
	}
	#-
	my $cin = &Getcomplementseq($in);
	for(my $i=0;$i<$len;$i++){
		$arr[$count]=substr($cin,$i,$len);
		$count++;
	}
	#sort

	@arr=sort(@arr);
	return $arr[0];
}

sub SSRSameDirectConsensusSequence{
	my $in=$_[0].$_[0];
	#print "$in\n";
	my @arr;
	my $len=length($_[0]);
	my $count=0;
	for(my $i=0;$i<$len;$i++){
		$arr[$count]=substr($in,$i,$len);
		$count++;
	}
	@arr=sort(@arr);
	return $arr[0];
}

sub SSR_SameDirect_ConsensusSequence{
	my $in=$_[0].$_[0];
	#print "$in\n";
	my @arr;
	my $len=length($_[0]);
	my $count=0;
	#+
	for(my $i=0;$i<$len;$i++){
		$arr[$count]=substr($in,$i,$len);
		$count++;
	}
	@arr=sort(@arr);
	return $arr[0];
}


sub ListandFindMaxSubtypeNumber{
	#ListsubtypeNumber 0 1 + -
	my $type;
	my %subtype;
	my %subtypename;
	my @number;	

	#@_ 1-type row 2-subtype row 3-4-5-6-sub type list
	$n=@_;
	#print OUT "Type\ttotal\ttypeno\t";
	print OUT "Type\t";
	for(my $i=3;$i<$n;$i++){
		$subtype{$_[$i]}=$i-2;
		$subtypename{$i-3}=$_[$i];
		print OUT "$_[$i] ";
	}
	print OUT "MixSubtype\ttMixSubtypeNo.\n";
	my $typeno=$n-3;
	@number=&arrayvalue($typeno,0);
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if ($type ne $tp[$_[1]]){
			if(length($type)>0){
				my $maxtype="";
				my $maxno=0;
				my $totalno=0;
				my $subtypeno=0;
				$l=@number;	
				for(my $i=0;$i<$l;$i++){
					if($number[$i]>0){
						$totalno+=$number[$i];
						$subtypeno++;						
					}
					if($number[$i]>$maxno){
						$maxno=$number[$i];
						$maxtype=$subtypename{$i};
					}elsif($number[$i]==$maxno){
						$maxtype=$maxtype." ".$subtypename{$i};
					}
				}			
				#print OUT "$type\t$totalno\t$subtypeno\t@number\t$maxtype\t$maxno\n";
				print OUT "$type\t@number\n";
			}
			$type = $tp[$_[1]];
			@number=&arrayvalue($typeno,0);		
		}
		my $loc=$subtype{$tp[$_[2]]};
		if($loc<=0){
			print "$tp[$_[2]] ";
		}else{
			$number[$loc-1]++;
		}		
	}

	my $totalno=0;
	my $subtypeno=0;
	my $maxtype="";
	my $maxno=0;	
	$l=@number;	
	for(my $i=0;$i<$l;$i++){
		if($number[$i]>0){
			$totalno+=$number[$i];
			$subtypeno++;
		}
		if($number[$i]>$maxno){
			$maxno=$number[$i];
			$maxtype=$subtypename{$i};
		}elsif($number[$i]==$maxno){
			$maxtype=$maxtype." ".$subtypename{$i};
		}		
	}	
	#print OUT "$type\t$totalno\t$subtypeno\t@number\t$maxtype\t$maxno\n";
	print OUT "$type\t@number\n";
	
}

sub ListsubtypeNumber{
	#ListsubtypeNumber 0 1 + -
	my $type;
	my %subtype;
	my @number;


	#@_ 1-type row 2-subtype row 3-4-5-6-sub type list
	$n=@_;
	print OUT "Type\ttotal\ttypeno\t";
	for(my $i=3;$i<$n;$i++){
		$subtype{$_[$i]}=$i-2;
		print OUT "$_[$i] ";
	}
	my $typeno=$n-3;
	@number=&arrayvalue($typeno,0);
	print OUT "\n";
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if ($type ne $tp[$_[1]]){
			if(length($type)>0){
				my $totalno=0;
				my $subtypeno=0;
				$l=@number;	
				for(my $i=0;$i<$l;$i++){
					if($number[$i]>0){
						$totalno+=$number[$i];
						$subtypeno++;
					}
				}			
				print OUT "$type\t$totalno\t$subtypeno\t@number\n";
			}
			$type = $tp[$_[1]];
			@number=&arrayvalue($typeno,0);		
		}
		my $loc=$subtype{$tp[$_[2]]};
		if($loc<=0){
			print "$tp[$_[2]] ";
		}else{
			$number[$loc-1]++;
		}		
	}
	my $totalno=0;
	my $subtypeno=0;
	$l=@number;	
	for(my $i=0;$i<$l;$i++){
		if($number[$i]>0){
			$totalno+=$number[$i];
			$subtypeno++;
		}
	}	
	print OUT "$type\t$totalno\t$subtypeno\t@number\n";
	
}



sub Listsubtypevalue{
	#Listsubtypevalue 0 1 7 all w9 r5 w8 w10 r7
	my $type;
	my %subtype;
	my @value;

	#@_ 1-type row 2-subtype row 3-subtype value row 4-5-6-sub type list
	$_=@_;
	print OUT "Type\t";
	for(my $i=4;$i<$_;$i++){
		$subtype{$_[$i]}=$i-4;
		print OUT "$_[$i] ";
	}
	my $no=$_-4;
	@value=&arrayvalue($no,0);
	print OUT "\n";
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if ($type ne $tp[$_[1]]){
			if(length($type)>0){
				print OUT "$type\t@value\n";
			}
			$type = $tp[$_[1]];
			@value=&arrayvalue($no,0);		
		}
		my $loc=$subtype{$tp[$_[2]]};
		if($loc<0){
			print "$tp[$_[2]] ";
		}else{
			$value[$loc]+=$tp[$_[3]];
		}		
	}
	print OUT "$type\t@value\n";
	
}

sub arrayvalue{
	my @r;
	for(my $i=0;$i<$_[0];$i++){
		$r[$i]=$_[1];
	}
	return @r;
}

sub SelectMutiRowbyValue{
	my @info;
	$_=@_;
	my $count=0;
	for(my $i=1;$i<$_;$i+=3){
		$info[$count][0]=$_[$i];
		$info[$count][2]=$_[$i+2];
		$info[$count][1]=$_[$i+3];
		$count++;
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $add=1;			
		for(my $i=0;$i<$count;$i++){
			if ($info[$i][2] eq "gt" ){
				if ($tp[$info[$i][0]] < $info[$i][1]){
					$add=0;
				}
			}	
			elsif ($info[$i][2] eq "=" ){
				if ($tp[$info[$i][0]] != $info[$i][1]){
					$add=0;
				}
			}
			elsif ($info[$i][2] eq "lt" ){
				if ($tp[$info[$i][0]] > $info[$i][1]){
					$add=0;
				}
			}
			elsif ($info[$i][2] eq "eq" ){
				if ($tp[$info[$i][0]] ne $info[$i][1]){
					$add=0;
				}
			}
			elsif ($info[$i][2] eq "ne" ){
				if ($tp[$info[$i][0]] eq $info[$i][1]){
					$add=0;
				}
			}
		}
		if($add==1){
			print OUT $lineIN."\n";
		}
	}

}


sub CutSequence_by_length{
	my $len=1;
	if($_[1]>0){
		$len=$_[1];
	}
	my $name;
	my $seq="";
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			if(length($seq)>0){
				$seq=substr($seq,0,$len);
				print OUT "$name\n$seq\n";
			}
			$name=$lineIN;
			$seq="";
		}else{
			$seq=$seq.$lineIN;
		}

    }
   	$seq=substr($seq,0,$len);
	print OUT "$name\n$seq\n";
}


sub CutfastafromEnd{
	my $len=3;
	my $name;
	my $seq="";
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			if(length($seq)>0){
				$seq=substr($seq,0,(length($seq)-3));
				print OUT "$name\n$seq\n";
			}
			$name=$lineIN;
			$seq="";
		}else{
			$seq=$seq.$lineIN;
		}

    }
   	$seq=substr($seq,0,(length($seq)-3));
	print OUT "$name\n$seq\n";
}

sub CutfastqfromEnd{
	my $len=0;
	my $name;
	while($lineIN=<IN>){
		$name=$lineIN;
		$lineIN=<IN>;
		chomp($lineIN);
		$len=length($lineIN);
		if($len>=$_[1]){
			print OUT $name;
			my $seq=substr($lineIN,0,($len-$_[1]));
			print OUT "$seq\n"; 
			$lineIN=<IN>;
			print OUT $lineIN;
			$lineIN=<IN>;
			$seq=substr($lineIN,0,($len-$_[1]));
			print OUT "$seq\n"; 			
		}else{
			$lineIN=<IN>;
			$lineIN=<IN>;
		}
    }
}

sub CutfastqfromStart{
	my $len=0;
	my $cut=25;
	my $name;
	
	my $outd=$ARGV[0].'.delete';
	open(OUTD,'>',"$outd");
	
	while($lineIN=<IN>){
		$name=$lineIN;
		$lineIN=<IN>;
		chomp($lineIN);
		$len=length($lineIN);
		if($len>=($_[1]+$cut)){
			print OUT $name;
			my $seq=substr($lineIN,($_[1]+1),$len);
			print OUT "$seq\n"; 
			$lineIN=<IN>;
			print OUT $lineIN;
			$lineIN=<IN>;
			chomp($lineIN);
			$seq=substr($lineIN,($_[1]+1),$len);
			print OUT "$seq\n"; 			
		}else{
			print OUTD "$name";
			$lineIN=<IN>;
			$lineIN=<IN>;
		}
   }
   close(OUTD);
}


sub CutfastqtoLength{      
	# $_[1] min $_[2] max
	my $name;
	
	my $outd=$ARGV[0].'.delete';
	open(OUTD,'>',"$outd");
	
	while($lineIN=<IN>){
		$name=$lineIN;
		$lineIN=<IN>;
		chomp($lineIN);
		$len=length($lineIN);
		if($len>=$_[1]){
			print OUT $name;   
			chomp($lineIN);
			my $seq=substr($lineIN,0,$_[2]);
			print OUT "$seq\n"; 
			$lineIN=<IN>;
			print OUT $lineIN;
			$lineIN=<IN>;
			chomp($lineIN);
			$seq=substr($lineIN,0,$_[2]);
			print OUT "$seq\n"; 			
		}else{
			print OUTD "$name";
			$lineIN=<IN>;
			$lineIN=<IN>;
		}
	}
	close(OUTD);
}


sub countfastqtotallen_GC{
	my $count=0;
	my $len=0;
	my $gc=0;
	my $at=0;
	my $nn=0;
	while($lineIN=<IN>){
		$count++;
		$lineIN=<IN>;
		chomp($lineIN);
		$len+=length($lineIN);
		for(my $i=0;$i<length($lineIN);$i++){
			my $aa=substr($lineIN,$i,1);
			if($aa eq 'G' || $aa eq 'g' || $aa eq 'C' || $aa eq 'c'){
				$gc++;
			}elsif($aa eq 'A' || $aa eq 'a' || $aa eq 'T' || $aa eq 't'){
				$at++;	
			}else{
				$nn++;
			}	
		}	
		$lineIN=<IN>;
		$lineIN=<IN>;
    }
    print OUT "$ARGV[0]\tTotal_Length $len\t";
    print OUT "Reads_number $count\t";
	print OUT "GC $gc\tAT $at N $nn\n";
}

sub ChangeFastaN{
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			print OUT $lineIN."\n";		
		}else{
			$lineIN =~ s/N/\*/;
			print OUT $lineIN."\n";	
		}
    }
}


sub RegulateFastaFile{
	my $cut=100;
	if($_[1]>0){
		$cut=$_[1];
	}
	my $seq="";
	my $name="";
	while($lineIN=<IN>){
		chomp($lineIN);

		if($lineIN=~/>/){
			#print "1";
			my $len=length($seq);
			if($len<1){
				print $name."\n";	
			}
			if($len>0){
				for(my $i=0;$i<$len;$i+=$cut){
					my $o=substr($seq,$i,$cut);
					$o =~ tr/[a-z]/[A-Z]/;
					print OUT "$o\n";
				}
			}
			$seq="";
			print OUT $lineIN."\n";	
			$name=$lineIN;
				
		}else{
			$seq=$seq.$lineIN;
		}
  }
  my $len=length($seq);
	for(my $i=0;$i<$len;$i+=$cut){
		my $o=substr($seq,$i,$cut);
		$o =~ tr/[a-z]/[A-Z]/;
		print OUT "$o\n";
	}
}

sub countFastqLength{
	my $count=0;
	my $name="NA";
	if($_[1] eq "d"){
		my $outt=$ARGV[0].".detail";
		print "Detail save to $outt";
		open(OUTT,'>',"$outt");	
		while($lineIN=<IN>){
			chomp($lineIN);
			$name = $lineIN;
			$lineIN=<IN>;
			chomp($lineIN);
			my $len=length($lineIN);
			$count+=$len;
			print OUTT "$name\t$len\n";
			$lineIN=<IN>;
			$lineIN=<IN>;	
    	}
    	close(OUTT);	
	
	}else{
		while($lineIN=<IN>){
			$lineIN=<IN>;
			chomp($lineIN);
			$count+=length($lineIN);
			$lineIN=<IN>;
			$lineIN=<IN>;	
    	}
    	
	}
	print OUT $count."\n";
    print "\nTotal Length $count\n";	
}

sub GiveEachLineMemberASymble{
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		for(my $i=1;$i<$tp;$i++){
			print OUT $tp[$i]."\t".$tp[0]."\n";
		}
		$count++;
    }
}

sub FreecBlastOverlapClear{
	#RL003_092	Chr7	24214	25295	1082	0	 0 0 0 0 0 0 0 0 0 0 0
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		print OUT "$lineIN\t";
		
		my $st=substr($tp[2],(length($tp[2])-2),2);
		print $st."..";
		my $count=0;
		my @no;
		if($st<20){
			print OUT "$tp[6] ";	
			$no[$count]=$tp[6];
			$count++;	
		}		
		
		for (my $i=7;$i<$tp-1;$i++){
			print OUT "$tp[$i] ";
			$no[$count]=$tp[$i];
			$count++;	
		}
		my $ed=substr($tp[3],(length($tp[3])-2),2);		
		if($ed>80){
			print OUT "$tp[$tp-1]\t";
			$no[$count]=$tp[$tp-1];
			$count++;		
		}else{
			print OUT "\t";	
		}
		print $ed."\t$count\t";
		@no=sort{$a <=> $b}(@no);
		
		#remove max and min get average number, if ><50% delete it, get avg again
		#print "@no\t";
		my $sc=0;
		my $sa;
		if($count>2){
			for(my $i=1;$i<$count-1;$i++){
				$sc+=$no[$i];
				#print $no[$i]." ";
			}
			$sa=$sc/($count-2);
		}elsif($count==2){
			$sa=($no[0]+$no[1])/2;
		}else{
			$sa=$no[0];	
		}
		
		my $sa2=$sa/2;
		my $sac=0;	
		my $scount=0;
		for(my $i=0;$i<$count;$i++){
			if($sa2<$no[$i] && ($sa2+$sa)>$no[$i]){
				$sac+=$no[$i];
				$scount++;
			}
		}	
		
		my $ot;
		if($scount>0){
			$ot=$sac/$scount;
			if(index($ot,".")>0){
				$ot=substr($ot,0,(index($ot,".")+3));	
			}
		}else{
			$ot="=";
		}
		
		print OUT "$ot\n";
		print "$ot\t$tp[5]\n";
    }
}

sub BlastnChangeTargettoLocal{
	#RLOth_323	chr10_13110	100.00	2655	0	0	1	2655	554894	557548	0.0	5263	4437 percent direct
	#add loss
	my $p=0;
	if($_[1]>0){
		$p=$_[1];
	}
	my $count=1;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		print OUT ">$tp[0]|$count|$tp[1]:";
		if ($tp[8]<$tp[9]){
			my $st=$tp[8]-$tp[6]+1-$p;
			my $ed=$tp[9]+$tp[12]-$tp[7]+$p;
			print OUT "$st..$ed\n";
		}else{
			my $st=$tp[9]-$tp[6]+1-$p;
			my $ed=$tp[8]+$tp[12]-$tp[7]+$p;
			print OUT "complement($st..$ed)\n";
		}
		
    }
	print "\n";
}


sub BlastnChangeQuerrytoLocal{
	#RUF_057n	J7L_1643	93.09	246	17	0	6876	7121	246	1	3.00E-97	353	459	0.535947712	-1
	#add loss
	my $p=0;
	if($_[1]>0){
		$p=$_[1];
	}
	my $count=1;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		print OUT ">$tp[0]|$count|$tp[0]:";
		if ($tp[8]<$tp[9]){
			my $st=$tp[6]-$tp[8]+1-$p;
			my $ed=$tp[7]+$tp[12]-$tp[9]+$p;
			print OUT "$st..$ed\n";
		}else{
			my $ed=$tp[7]-$tp[9]+1+$p;
			my $st=$tp[6]-($tp[12]-$tp[8])-$p;
			print OUT "complement($st..$ed)\n";
		}		
    }
	print "\n";
}


sub BlastResultCombin{
	#RLOth_323	chr10_13110	100.00	2655	0	0	1	2655	554894	557548	0.0	5263	4437
	#add length to last line
	#if distance < total length, combin
	my $count=0;
	my @save;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		my $tpd;
		if($tp[8]>$tp[9]){
			$tpd=-1;
		}else{
			$tpd=1;
		}
		if($save[0] eq $tp[0]){
			#1-2 before corse 1 
			#2-1 before corse 2 
			#3-2 in 1 
			#4-1 in 2 
			#5-same loc 
			#6-no corse 1 before 2 
			#7-no corse 2 before 1
			my @c1=&Comparelocal($save[1],$save[6],$save[7],$tp[1],$tp[6],$tp[7]);
			my @c2=&Comparelocal($save[1],$save[8],$save[9],$tp[1],$tp[8],$tp[9]);			
			if($save[14]==1 && $tpd==1){
				if($c1[0]==$c2[0] && $c2[1]<$c1[1]*2){
					if($c1[0]==2 || $c1[0]==6){
						$save[7]=$tp[7];
						$save[9]=$tp[9];
						$save[3]=$save[9]-$save[8];
					}
					if($c1[0]==1 || $c1[0]==7){
						$save[6]=$tp[6];
						$save[8]=$tp[8];
					}
					if($c1[0]==4){
						@save=@tp;				
					}	
					$save[3]=$save[9]-$save[8];	
					$save[13]=$save[3]/$save[12];		
				}
			}elsif($save[14]==-1 && $tpd==-1){
				if($c2[1]<$c1[1]*2){
					if(($c1[0]==2 && $c2[0]==1) ||($c1[0]==6 && $c2[0]==7)){
						$save[7]=$tp[7];
						$save[9]=$tp[9];
					}
					if(($c1[0]==1 && $c2[0]==2) ||($c1[0]==7 && $c2[0]==6)){
						$save[6]=$tp[6];
						$save[8]=$tp[8];
					}
					if($c1[0]==4){
						@save=@tp;				
					}
					$save[3]=$save[8]-$save[9];	
					$save[13]=$save[3]/$save[12];						
				}			
			}else{
				print ".";
			}
			
		}else{
			if($count>0){
				#output
				print OUT "@save\t\n";
			}
			@save=@tp;
			$save[14]=$tpd;
			$save[13]=$save[3]/$save[12];

		}
		$count++;
		if($count%10000==0){
			print $count."\n";
		}
		
    }
    print OUT "@save\t\n";
	print "\n";
}

sub BlastDeleteRedundance{
	#if overlap, leave high e-value 
	my $count=0;
	my @save;
	while($lineIN=<IN>){
		if($count ==0){
			chomp($lineIN);
			@save= &Getheadinfo($lineIN);
			$lineIN=<IN>;
			$count++;
		}
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $o=&FindnosortOverlap($save[8],$save[9],$tp[8],$tp[9]);
		print $o;
		if($o==1){
			if($save[2]<$tp[2]){
				@save=@tp;
			}
		}else{
			print OUT "@save\t\n";
			@save=@tp;
		}
    }
    print OUT "@save\t\n";
    print "\n";
}

sub FindOverlap{ #s1 sorted
	#0-s1, 1-e1, 2-s2, 3-e2
	if($_[1]>=$_[2] && $_[1]<=$_[3]){
		return 1;
	}	
	elsif($_[0]>=$_[2] && $_[0]<=$_[3]){
		return 1;
	}
	elsif($_[2]>=$_[0] && $_[2]<=$_[1]){
		return 1;
	}
	elsif($_[3]>=$_[0] && $_[3]<=$_[1]){
		return 1;
	}	
	return 0;
}

sub OverlapLength{ #sorted  s2 > s1
	#0-s1, 1-e1, 2-s2, 3-e2
	my $l=-1;
	if($_[1]>=$_[2]){
		if($_[3] >= $_[1]){
			$l=$_[1]-$_[2] +1;
		}else{
			$l=$_[3]-$_[2] +1;
		}
	}
	return $l;
}

sub FindnosortOverlap{ #no sort
	#0-s1, 1-e1, 2-s2, 3-e2
	if($_[1]>=$_[2] && $_[1]<=$_[3]){
		return 1;
	}
	if($_[0]>=$_[2] && $_[0]<=$_[3]){
		return 1;
	}
	return 0;
}

sub ChangeRowbyCompair{
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		my $tp=@tp;
		if($tp[$_[1]] <= $tp[$_[2]]){
			for(my $i=0;$i<$tp;$i++){
				print OUT $tp[$i]."\t";
			}
			print OUT "+\n";
		}else{
			for(my $i=0;$i<$tp;$i++){
				if($i ==$_[1]){
					print OUT $tp[$_[2]]."\t";
				}elsif($i ==$_[2]){
					print OUT $tp[$_[1]]."\t";
				}else{
					print OUT $tp[$i]."\t";
				}
			}
			print OUT "-\n";		
		}
    }
}


sub SeperateFilebySortedRowType{  #sorted
	my $name='';
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if( $name eq ''){
			$name=$tp[$_[1]];	
			open(OUTN,'>',"$name");	
		}elsif ($tp[$_[1]] ne $name){		
			close(OUTN);
			print "$name\n";
			$name=$tp[$_[1]];
			open(OUTN,'>',"$name");	
		}
		print OUTN $lineIN."\n";
    }
    print "$name\n";
	close(OUTN);
}

sub SeperateFilebyRowType{  
	my %name;
	my %fname;
	my %oname;
	my $count=0;
	
	while($lineIN=<IN>){
		chomp($lineIN);
		#my @tp= &Getheadinfo($lineIN);
		my @tp = split(/\t/,$lineIN);
		if(!$name{$tp[$_[1]]}){
			$count++;
			$name{$tp[$_[1]]}=$tp[$_[1]];
			#$fname{$tp[$_[1]]}=$ARGV[0].'.'.$tp[$_[1]];
			$fname{$tp[$_[1]]}=$tp[$_[1]];
			$oname{$tp[$_[1]]}="OUT".$tp[$_[1]];
			open ($oname{$tp[$_[1]]},'>',"$fname{$tp[$_[1]]}");	
			print "$name{$tp[$_[1]]}\t$oname{$tp[$_[1]]}\t$fname{$tp[$_[1]]}\n";
			print OUT "$fname{$tp[$_[1]]}\n";
		}else{
			"aaa";		
		}
		my $a=$oname{$tp[$_[1]]};
		print $a "$lineIN\n";
		#print $a "$tp[1]\n";
  }
	my @onamek=keys(%oname);
  for(my $i=0;$i<$count;$i++){
		close($oname{$onamek[$i]});
	}
	print "$count filed was savered\n";
}




sub Seperate_cut_FastAbyLength{
	my $cut=500000;
	if($_[1]>0){
		$cut=$_[1];
	}
	my $seq="";
	my $name="";
	my $st=1;
	while($lineIN=<IN>){
		chomp($lineIN);

		if($lineIN=~/>/){
			#print "1";
			my $len=length($seq);

			if($len>0){
				for(my $i=0;$i<$len;$i+=$cut){
					my $o=substr($seq,$i,$cut);
					$o =~ tr/[a-z]/[A-Z]/;
					my $ed=$st+$cut-1;
					$ed=$len<$ed?$len:$ed;
					print OUT "$name:$st-$ed\n$o\n";
					$st=$st+$cut;
				}
			}
			$seq="";
			$st=1;
			$name=$lineIN;
				
		}else{
			$seq=$seq.$lineIN;
		}
  }
  my $len=length($seq);
	for(my $i=0;$i<$len;$i+=$cut){
		my $o=substr($seq,$i,$cut);
		$o =~ tr/[a-z]/[A-Z]/;
		my $ed=$st+$cut-1;
		$ed=$len<$ed?$len:$ed;
		print OUT "$name:$st..$ed\n$o\n";
		$st=$st+$cut;
	}
}


sub SeperateFastAbyLength{
	#group seperate by |
	my $n=500000000;
	if($_[1] > 0){
		$n=$_[1];
	}	

	my $count=0;
	my $len=0;
	print "each file contain large than $n bp\n";	
	while($lineIN=<IN>){		
		chomp($lineIN);
		if ($lineIN=~/>/){
			$tpn++;
			if($len >= $n || $count == 0){
				$count++;
				print "$count..";
				$len=0;
				if($count > 0){
					close(OUTN);
				}
				my $outn=$ARGV[0].".".$count;
				open(OUTN,'>',"$outn");
				print OUT "$outn\n";	
				print "$outn\n";
						
			}	
			print OUT "\t$lineIN\n";		
		}else{
			my $tplen=length($lineIN);
			$len+=$tplen;
		}
		print OUTN "$lineIN\n";
    }
	print "\ndone\n";
	close(OUTN);
}


sub SeperatFastAbyNumber{
	#group seperate by |
	my $n=1;
	if($_[1] > 0){
		$n=$_[1];
	}	

	my $count=0;
	my $tpn=0;
	print "each file contain $n seq\n";	
	while($lineIN=<IN>){		
		chomp($lineIN);
		if ($lineIN=~/>/){
			$tpn++;
			if($tpn >= $n || $count == 0){
				$count++;
				print "$count..";
				$tpn=0;
				if($count > 0){
					close(OUTN);
				}
				my $outn=$ARGV[0].".".$count;
				open(OUTN,'>',"$outn");	
						
			}			
		}
		print OUTN "$lineIN\n";
    }
	print "\ndone\n";
	close(OUTN);
}


sub SeperateGroupFastAtoFile{
	#group seperate by

	my $type="";
	my $seq='';
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN, "|", "=");
		if ($tp[0]=~/>/){
			if($type ne '' && $type ne $tp[0]){				
				my $name=substr($type,1,length($type)).".fas";
				if(-e $name){
					open(OUTN,'>>',"$name");	
				}else{
					print "$type\n";
					open(OUTN,'>',"$name");	
				}
				print OUTN $seq;
				close(OUTN);
				$seq="";			
			}	
			$type=$tp[0];		
		}
		$seq=$seq."$lineIN\n";
    }
	print "$type\n";
	my $name=substr($type,1,length($type)).".fas";
	if(-e $name){
		open(OUTN,'>>',"$name");	
	}else{
		open(OUTN,'>',"$name");	
	}
	print OUTN $seq;
	close(OUTN);
}





sub GetEnsemblGeneAnnCDS{
	my @gene;
	#GTF
	#1	ensembl	gene	29141326	29149072	.	-	.	gene_id "ENSSSCG00000004150"; gene_version "2"; gene_name "HEBP2"; gene_source "ensembl"; gene_biotype "protein_coding";
	#1	ensembl	transcript	29141326	29149072	.	-	.	gene_id "ENSSSCG00000004150"; gene_version "2"; transcript_id "ENSSSCT00000004587"; transcript_version "2"; gene_name "HEBP2"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "HEBP2-201"; transcript_source "ensembl"; transcript_biotype "protein_coding";
	#1	ensembl	exon	29148452	29149072	.	-	.	gene_id "ENSSSCG00000004150"; gene_version "2"; transcript_id "ENSSSCT00000004587"; transcript_version "2"; exon_number "1"; gene_name "HEBP2"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "HEBP2-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSSSCE00000034388"; exon_version "2";
	#1	ensembl	CDS	29148452	29148553	.	-	0	gene_id "ENSSSCG00000004150"; gene_version "2"; transcript_id "ENSSSCT00000004587"; transcript_version "2"; exon_number "1"; gene_name "HEBP2"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "HEBP2-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSSSCP00000004484"; protein_version "2";
	#0-chrom 2-gene 3-st 4-ed 6-drect 9-gene_id 13-gene_name/transcript_id 17-gene_biotype
	#0-chrom 2-CDS/transcript 3-st 4-ed 6-drect 9-gene_id 13-transcript_id
	#out
	#0-chrom 1-gene 2-st 3-ed 4-direct 5-name 6-type 7-transcript_no 8-transcript_id 9-st 10-ed 11-CDS_no 12-CDS_st1 13-CDS_ed1 14... 	
		
	my $gene="NA";
	my @transcript;   #0-transcript_id 1-cds_no 2-cds_st1 cds_ed2 ...
	my $tc=0;
	my $count=0;

	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, "\"",";");
		my $tp=@tp;
		if($tp[2]=~/gene/i){			
			$count++;
			if($count%100 == 0){
				print ".";
			}
			if($gene ne "NA"){
				$tc++;
				for(my $i=0;$i<$tc;$i++){
					print OUT "$gene$tc\t$transcript[$i][0]$transcript[$i][1]$transcript[$i][2]\n";
				}
			}
			#$gene="$tp[0]\t$tp[3]\t$tp[4]\t$tp[6]\t$tp[9]\t$tp[13]\t$tp[17]\t";
			$gene="$tp[0]\t$tp[3]\t$tp[4]\t$tp[6]\t$tp[9]\t";
			$tc=-1;
			@transcript=();
		}elsif($tp[2] eq "transcript"){
			$tc++;
			#add new
			$transcript[$tc][0]="$tp[13]\t$tp[3]\t$tp[4]\t";
			$transcript[$tc][1]=0;
			$transcript[$tc][2]="";
		}elsif($tp[2] eq "CDS"){
			$transcript[$tc][1]++;
			$transcript[$tc][2]=$transcript[$tc][2]."\t$tp[3]\t$tp[4]";
		}
		else{  			
			#
		}
	}
	
	$tc++;
	for(my $i=0;$i<$tc;$i++){
		print OUT "$gene$tc\t$transcript[$i][0]$transcript[$i][1]$transcript[$i][2]\n";
	}	
	print "\nTotal $count gene\n";

}



sub GetGffanno_by_index{
	#NC_000932.1	RefSeq	CDS	109405	110436	.	+	0	ID=cds-NP_051105.1;Parent=gene-ArthCp070;Dbxref=Genbank:NP_051105.1,GeneID:844706;Name=NP_051105.1;gbkey=CDS;gene=ycf1;locus_tag=ArthCp070;product=hypothetical protein;protein_id=NP_051105.1;transl_table=11
	#for example get ID=cds- gene= Parent=rna-
	
	#index	
	my @syno = @_;
	splice(@syno,0,1);
	my $count=@syno;
	print OUT join("\t",@syno)."\n";
	#print join ("\t",@syno);
	#hash used to save outputed info
	my %ids;

	while($lineIN=<IN>){
		chomp($lineIN);	
		#my @tp= &Getheadinfo($lineIN,";");
		#my @tp = split(/\t,;/,$lineIN);
		my @out;
		for(my $i=0; $i<$count;$i++){
			my $out="NA";
			if(index($lineIN,$syno[$i],)>0){
				my $tp=substr($lineIN,(index($lineIN,$syno[$i],)+length($syno[$i])),length($lineIN));
				my @info = split(/;/,$tp);				
				if($i==0){
					if($ids{$info[0]}==1){
						#print ".";
						last;
					}
				}
				$ids{$info[0]}=1;
				$out[$i]=$info[0]; 
			}else{
				last;
			}			 
		}
		if($out[0]){
			print OUT join("\t",@out)."\n";
		}
	}	
	print "\n";
}



sub GetGffanno_by_name{
	#NW_002477246.1	RefSeq	mRNA	56004	57442	.	+	.	ID=rna20;Parent=gene20;Dbxref=Genbank:XM_002381477.1,GeneID:7920930;Name=XM_002381477.1;Note=transcript AFLA_123240A;end_range=57442,.;gbkey=mRNA;partial=true;product=conserved hypothetical protein;start_range=.,56004;transcript_id=XM_002381477.1
	my $syno=@_;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN,";");
		#my @tp = split(/\t,;/,$lineIN);
		my $tp=@tp;
		print OUT "$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]";
		for(my $i=1; $i<$syno;$i++){
			my $add=0;
			for(my $j=8; $j<=$tp;$j++){
				if(index($tp[$j],$_[$i],0)==0){
					$add=$j;
					my $st=length($_[$i])+1;
					my $o=substr($tp[$j],$st,length($tp[$j]));
					print OUT "\t$o";
					last;
				}
			}
			if($add==0){
				print OUT "\tNA";
			}
		}
		print OUT "\n";
	}	
}




sub GetGff3Gene_mRNA_CDS_name{

#	NC_008394.4	RefSeq	gene	2167	9751	.	+	.	ID=gene0;Dbxref=GeneID:4326813;Name=Os01g0100100;gbkey=Gene;gene=Os01g0100100;gene_biotype=protein_coding
#	NC_008394.4	RefSeq	mRNA	2167	9751	.	+	.	ID=rna0;Parent=gene0;Dbxref=Genbank:NM_001048268.2,GeneID:4326813;Name=NM_001048268.2;Note=supported by CI317886 and CI567119;gbkey=mRNA;gene=Os01g0100100;transcript_id=NM_001048268.2
#	NC_008394.4	RefSeq	exon	2167	2268	.	+	.	ID=id3;Parent=rna0;Dbxref=Genbank:NM_001048268.2,GeneID:4326813;Note=supported by CI317886 and CI567119;gbkey=mRNA;gene=Os01g0100100;transcript_id=NM_001048268.2
#	NC_008394.4	RefSeq	CDS	2449	2616	.	+	0	ID=cds0;Parent=rna0;Dbxref=Genbank:NP_001041733.2,GeneID:4326813;Name=NP_001041733.2;Note=RabGAP/TBC domain containing protein~contains InterPro domain(s): IPR000195~supported by CI317886 and CI567119;gbkey=CDS;gene=Os01g0100100;product=hypothetical protein;protein_id=NP_001041733.2

	my $count=0;
	my $gene_name="";
	my $mRNA_name="";
	my $CDS_name="";
	my $add=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		my $name="";
		for(my $i=0;$i<$tp;$i++){
			if($tp[$i]=~/Name/){
				$name=$tp[$i+1];
				#print "$tp[$i] $name\n";
				last;
			}
		}		
		#print "$tp[2]";
		if($tp[2] eq "gene"){
			#nocds gene
			#if($info[$count] ne "" && $add==0){
			#	print OUT "$gene_name\t$mRNA_name\n";
			#}
			#add new
			$count++;
			$add=0;
			$gene_name=$name;
			$mRNA_name="";
			$CDS_name="";

		}elsif($tp[2] eq "mRNA"){
			$mRNA_name=$name;
			
		}elsif($tp[2] eq "CDS" && $add==0){
			$CDS_name=$name;
			print OUT "$gene_name\t$mRNA_name\t$CDS_name\n";
			$add=1;
		}
	}

	
}




sub GetGff3GeneInfo{
	my @gene;
	#0-scaf 1-start 2-end 3-ID 4-Direct 5-Exon_No 6-s1 7-e1 8-s2 9-e2 ...
	#gff3
	#scaffold619	.	gene	368802	373922	.	-	.	ID=LOC_Os01g37120.1.gene.63518;Name=LOC_Os01g37120.1.63518

	my $count=-1;
	my @info;
	my @exon;
	my $name;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		if($tp[2] eq "gene"){
			#output
			print "+";
			if($info[$count] ne ""){
				#add exon info
				my $exon=@exon;
				my $exonno=$exon/2;
				$info[$count]=$info[$count]."\t".$exonno;
				for(my $i=0; $i<$exon;$i++){
					$info[$count]=$info[$count]."\t".$exon[$i];
				}
			}
			#add new
			$count++;
			@exon=();
			$info[$count]="$tp[0]\t$tp[3]\t$tp[4]\t$tp[6]\t";
			for(my $i=0; $i<$tp;$i++){
				if($tp[$i] eq "ID"){
					$info[$count]=$info[$count].$tp[$i+1];
					$name=$tp[$i+1];
					last;
				}
			}
			
		}else{
			if($tp[2] eq "CDS" && index($lineIN,$name)>=0){	
				push(@exon,($tp[3],$tp[4]));
				print ".";
			}
		}
	}
	my $exon=@exon;
	my $exonno=$exon/2;
	$info[$count]=$info[$count]."\t".$exonno;
	for(my $i=0; $i<$exon;$i++){
		$info[$count]=$info[$count]."\t".$exon[$i];
	}
	$count++;
	print "\nTotal $count gene\n";
	
	@info=sort(@info);
	for(my $i=0; $i<$count;$i++){
		print OUT $info[$i]."\n";
	}
	
}



sub GetGff3Gene_exon_length_jv6{
	#gff3
	#scaffold619	.	gene	368802	373922	.	-	.	ID=LOC_Os01g37120.1.gene.63518;Name=LOC_Os01g37120.1.63518
	#if have alternative splicing, only use first record
	
	my $count=-1;
	my $name="NA";
	my $len;
	my $exon;
	my $canadd=2;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		if($tp[2] eq "gene"){
			#output
			print "+";
			if($name ne "NA"){
				#add exon info
				print OUT "$name\t$len\t$exon\n";
			}
			#add new
			$count++;
			$exon=0;
			$canadd=0;
			$len=$tp[4]-$tp[3]+1;
			for(my $i=0; $i<$tp;$i++){
				if($tp[$i] eq "ID"){
					$name=$tp[$i+1];
					last;
				}
			}			
		}else{
			if($tp[2] eq "CDS"){	
				print ".";
				if($canadd==0){
					$canadd=1;
				}
				if($canadd==1){
					$exon=$exon+$tp[4]-$tp[3]+1;
				}					
				
			}else{
				if($canadd==1){
					$canadd=2;
				}	
			}	
		}
	}
	print OUT "$name\t$len\t$exon\n";
	print "$count\n";
}

sub GetGff3out_to_gdb{
	#input #0-scaf 1-start 2-end 3-ID 4-Direct 5-5UTR_No 6-Exon_No 7-3UTR_No 8-s1 9-e1 10-s2 11-e2 ...
	#Chr1	11218	12435	LOC_Os01g01019	LOC_Os01g01019.1	+	1	2	1	11218	11797	11798	12060	12152	12317	12318	12435	
	#Chr1|LOC_Os01g01019|+|11218|12435|12060||4|11797, 12060, 12317, 12435, |11798, 12152, 12318, , |LOC_Os01g01019.1

	## output gdb
	#CREATE TABLE GeneTable (chrom,name,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,name2);
	#chrI|NM_058260|-|4123|10232|4220|10148|5|4123, 5194, 6036, 9726, 10094,|4358, 5296, 6327, 9846, 10232,|Y74C9A.3
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;		
		my $exonStarts="";
		my $exonEnds="";
		my $exonCount=$tp[6]+$tp[7]+$tp[8];
		#+
		my $cdsStart=$tp[9+($tp[6]*2)];		
		my $cdsEnd=$tp[8+(($tp[6]+$tp[7])*2)];
		#-
		if($tp[5] eq "-"){
			$cdsStart=$tp[9+($tp[8]*2)];	
			$cdsEnd=$tp[8+(($tp[8]+$tp[7])*2)];
		}
		
		for(my $i=9;$i<$tp;$i+=2){
			$exonStarts=$exonStarts.$tp[$i].", ";
			$exonEnds=$exonEnds.$tp[$i+1].", ";
		}
		
		print OUT "$tp[0]|$tp[3]|$tp[5]|$tp[1]|$tp[2]|$cdsStart|$cdsEnd|$exonCount|$exonStarts|$exonEnds|$tp[4]\n";
	}	
}


sub GetGff3CDSInfo_jv7{
	my @gene;
	#0-scaf 1-start 2-end 3-ID 4-Direct 5-5UTR_No 6-Exon_No 7-3UTR_No 8-s1 9-e1 10-s2 11-e2 ...
	#if +
	#if -
	#8-9-10-11- change direct
	#gff3
	#scaffold619	.	gene	368802	373922	.	-	.	ID=LOC_Os01g37120.1.gene.63518;Name=LOC_Os01g37120.1.63518
	#Chr1	MSU_osa1r6	gene	1903	9817	.	+	.	ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
	#if have alternative splicing, only use first record

	
	my $count=-1;
	my @gene;
	my @exon;
	my @UTR5;
	my @UTR3;
	my $name;
	my $father;
	my $d;
	#count
	my $fl;
	my $el;
	my $il;
	my $u5l;
	my $u3l;
	my $gl;
	my $canadd=2;

	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		if($tp[2] eq "gene"){
			$father=$tp[9];
		}elsif($tp[2] eq "mRNA"){
			#output
			print "+";
			if($count>=0){
				#add exon info
				my $exon=@exon;
				my $exonno=$exon/2;
				my $UTR5=@UTR5;
				#print "$UTR5";
				my $UTR5no=0;
				if($UTR5>0){
					$UTR5no=$UTR5/2;
				}
				my $UTR3=@UTR3;
				#print "$UTR3";
				my $UTR3no=0;
				if($UTR3>0){
					$UTR3no=$UTR3/2;
				}
				$info[$count]=$info[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
				if($d eq "+"){
					#length
					if($UTR5>0){
						$gl=$gl+($UTR5[0]-$gene[0])+($exon[0]-$UTR5[$UTR5-1]-1);
						$u5l=$UTR5[1]-$UTR5[0]+1;
						for(my $i=2; $i<$UTR5;$i+=2){
							$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
							$g+=$UTR5[$i]-$UTR5[$i-1]-1;}					
					}else{  
						$gl+=$exon[0]-$gene[0];
					}
					$el=$exon[1]-$exon[0];		
					for(my $i=2; $i<$exon;$i+=2){
						$el+=$exon[$i+1]-$exon[$i+1];
						$il+=$exon[$i]-$exon[$i-1]-1;
					}
					if($UTR3>0){
						$gl=$gl+($gene[1]-$UTR3[$UTR3-1])+($UTR3[0]-$exon[$exon-1]-1);
						$u3l=$UTR3[1]-$UTR3[0]+1;
						for(my $i=2; $i<$UTR3;$i+=2){
							$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
							$g+=$UTR3[$i]-$UTR3[$i-1]-1;}					
					}else{  
						$gl+=$gene[1]-$exon[$exon-1];
					}
										
					#info
					for(my $i=0; $i<$UTR5;$i++){
						$info[$count]=$info[$count]."\t".$UTR5[$i];
					}
					for(my $i=0; $i<$exon;$i++){
						$info[$count]=$info[$count]."\t".$exon[$i];
					}
					for(my $i=0; $i<$UTR3;$i++){
						$info[$count]=$info[$count]."\t".$UTR3[$i];
					}
				}else{
					@exon = sort{$a <=> $b} @exon;	
					@UTR3 = sort{$a <=> $b} @UTR3;
					@UTR5 = sort{$a <=> $b} @UTR5;
					
					#length
					if($UTR3>0){
						$gl=$gl+($UTR3[0]-$gene[0])+($exon[0]-$UTR3[$UTR3-1]-1);
						$u3l=$UTR3[1]-$UTR3[0]+1;
						for(my $i=2; $i<$UTR3;$i+=2){
							$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
							$g+=$UTR3[$i]-$UTR3[$i-1]-1;}					
					}else{  
						$gl+=$exon[0]-$gene[0];
					}
					$el=$exon[1]-$exon[0];		
					for(my $i=2; $i<$exon;$i+=2){
						$el+=$exon[$i+1]-$exon[$i]+1;
						$il+=$exon[$i]-$exon[$i-1]-1;
					}
					if($UTR5>0){
						$gl=$gl+($gene[1]-$UTR5[$UTR5-1])+($UTR5[0]-$exon[$exon-1]-1);
						$u5l=$UTR5[1]-$UTR5[0]+1;
						for(my $i=2; $i<$UTR5;$i+=2){
							$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
							$g+=$UTR5[$i]-$UTR5[$i-1]-1;}					
					}else{  
						$gl+=$gene[1]-$exon[$exon-1];
					}
										
					#info					
					for(my $i=0; $i<$UTR3;$i++){
						$info[$count]=$info[$count]."\t".$UTR3[$i];
					}	
					for(my $i=0; $i<$exon;$i++){
						$info[$count]=$info[$count]."\t".$exon[$i];
					}
					for(my $i=0; $i<$UTR5;$i++){
						$info[$count]=$info[$count]."\t".$UTR5[$i];
					}				
				}
				print OUT $info[$count]."\n";
				$fl=$gene[1]-$gene[0]+1;
			}			
			
			#add new
			$d=$tp[6];
			$count++;	 
			@exon=();
			@UTR5=();
			@UTR3=();
			$el=0;
			$il=0;
			$gl=0;
			$u5l=0;
			$u3l=0;
			$canadd=0;		
			$gene[0]=$tp[3];
			$gene[1]=$tp[4];
			for(my $i=0; $i<$tp;$i++){
				if($tp[$i] eq "ID"){
					$info[$count]="$tp[0]\t$tp[3]\t$tp[4]\t$father\t$tp[$i+1]\t$tp[6]";
					$name=$tp[$i+1];
					#print "$tp[$i+1] $name ";
					last;
				}
			}
			
		}
		else{  
			if($tp[2] eq "CDS" && index($lineIN,$name)>=0){	
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@exon,($tp[3],$tp[4]));
				}
			}elsif($tp[2] eq "five_prime_UTR" && index($lineIN,$name)>=0){
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@UTR5,($tp[3],$tp[4]));
					my $UTR5=@UTR5;
					print "$UTR5";
				}				
			}elsif($tp[2] eq "three_prime_UTR" && index($lineIN,$name)>=0){
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@UTR3,($tp[3],$tp[4]));
					my $UTR3=@UTR3;
					print "$UTR3";
				}				
			}
			else{
				if($canadd==1){
					$canadd=2;
				}	
			}
		}
	}
	print "\nTotal $count gene\n";
	
}





sub GetGff3GeneInfo_jv7{
	my @gene;
	#0-scaf 1-start 2-end 3-ID 4-Direct 5-Exon_No 6-s1 7-e1 8-s2 9-e2 ...
	#gff3
	#scaffold619	.	gene	368802	373922	.	-	.	ID=LOC_Os01g37120.1.gene.63518;Name=LOC_Os01g37120.1.63518
	#if have alternative splicing, only use first record
	
	my $count=-1;
	my @info;
	my @exon;
	my $name;
	my $canadd=2;
	
	my $outd=$ARGV[0].'.detail';
	open(OUTD,'>',"$outd");	
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		if($tp[2] eq "gene"){
			#output
			#print "+";
			if($info[$count] ne ""){
				#add exon info
				my $exon=@exon;
				my $exonno=$exon/2;
				$info[$count]=$info[$count]."\t".$exonno;
				for(my $i=0; $i<$exon;$i++){
					$info[$count]=$info[$count]."\t".$exon[$i];
				}
			}
			#add new
			$count++;
			@exon=();
			$canadd=0;
			$info[$count]="$tp[0]\t$tp[3]\t$tp[4]\t$tp[6]\t";
			for(my $i=0; $i<$tp;$i++){
				if($tp[$i] eq "Name"){
					#print "a";
					$info[$count]=$info[$count].$tp[$i+1];
					$name=$tp[$i+1];
					last;
				}
			}
			
		}else{
			if(($tp[2] eq "five_prime_UTR" || $tp[2] eq "three_prime_UTR" || $tp[2] eq "CDS") && index($lineIN,$name)>=0){	
				if($canadd==0){
					$canadd=1;
				}
				if($canadd==1){
					push(@exon,($tp[3],$tp[4]));
					if($tp[2] eq "five_prime_UTR"){
						print OUTD "$tp[0]\t$tp[3]\t$tp[4]\t$name\tfive_prime_UTR\n";
					}elsif($tp[2] eq "three_prime_UTR"){
						print OUTD "$tp[0]\t$tp[3]\t$tp[4]\t$name\tthree_prime_UTR\n";
					}elsif($tp[2] eq "CDS"){
						print OUTD "$tp[0]\t$tp[3]\t$tp[4]\t$name\tCDS\n";
					}
					
					
					#print ".";
				}
			}else{
				if($canadd==1){
					$canadd=2;
				}	
			}
		}
	}
	my $exon=@exon;
	my $exonno=$exon/2;
	$info[$count]=$info[$count]."\t".$exonno;
	for(my $i=0; $i<$exon;$i++){
		$info[$count]=$info[$count]."\t".$exon[$i];
	}
	$count++;
	print "\nTotal $count gene\n";
	
	@info=sort(@info);
	for(my $i=0; $i<$count;$i++){
		print OUT $info[$i]."\n";
	}
	close(OUTD);
	
}


sub GetGff3GeneInfo_AlternativeSplicing_jv6{
	my @gene;
	#0-scaf 1-start 2-end 3-ID 4-Direct 5-5UTR_No 6-Exon_No 7-3UTR_No 8-s1 9-e1 10-s2 11-e2 ...
	#if +
	#if -
	#8-9-10-11- change direct
	#gff3
	#scaffold619	.	gene	368802	373922	.	-	.	ID=LOC_Os01g37120.1.gene.63518;Name=LOC_Os01g37120.1.63518
	#Chr1	MSU_osa1r6	gene	1903	9817	.	+	.	ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
	#if have alternative splicing, only use first record

	my $count=-1;
	my @gene;
	my @exon;
	my @UTR5;
	my @UTR3;
	my $name;
	my $father;
	my $d;
	#count
	my $fl;
	my $el;
	my $il;
	my $u5l;
	my $u3l;
	my $gl;
	my $canadd=2;

	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		if($lineIN eq "###"){
			
			#output
			print "+";
			if($count>=0){
				#add exon info
				my $exon=@exon;
				my $exonno=$exon/2;
				my $UTR5=@UTR5;
				#print "$UTR5";
				my $UTR5no=0;
				if($UTR5>0){
					$UTR5no=$UTR5/2;
				}
				my $UTR3=@UTR3;
				#print "$UTR3";
				my $UTR3no=0;
				if($UTR3>0){
					$UTR3no=$UTR3/2;
				}
				$info[$count]=$info[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
				if($d eq "+"){
					#length
					if($UTR5>0){
						$gl=$gl+($UTR5[0]-$gene[0])+($exon[0]-$UTR5[$UTR5-1]-1);
						$u5l=$UTR5[1]-$UTR5[0]+1;
						for(my $i=2; $i<$UTR5;$i+=2){
							$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
							$g+=$UTR5[$i]-$UTR5[$i-1]-1;}					
					}else{  
						$gl+=$exon[0]-$gene[0];
					}
					$el=$exon[1]-$exon[0];		
					for(my $i=2; $i<$exon;$i+=2){
						$el+=$exon[$i+1]-$exon[$i+1];
						$il+=$exon[$i]-$exon[$i-1]-1;
					}
					if($UTR3>0){
						$gl=$gl+($gene[1]-$UTR3[$UTR3-1])+($UTR3[0]-$exon[$exon-1]-1);
						$u3l=$UTR3[1]-$UTR3[0]+1;
						for(my $i=2; $i<$UTR3;$i+=2){
							$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
							$g+=$UTR3[$i]-$UTR3[$i-1]-1;}					
					}else{  
						$gl+=$gene[1]-$exon[$exon-1];
					}
										
					#info
					for(my $i=0; $i<$UTR5;$i++){
						$info[$count]=$info[$count]."\t".$UTR5[$i];
					}
					for(my $i=0; $i<$exon;$i++){
						$info[$count]=$info[$count]."\t".$exon[$i];
					}
					for(my $i=0; $i<$UTR3;$i++){
						$info[$count]=$info[$count]."\t".$UTR3[$i];
					}
				}else{
					@exon = sort{$a <=> $b} @exon;	
					@UTR3 = sort{$a <=> $b} @UTR3;
					@UTR5 = sort{$a <=> $b} @UTR5;
					
					#length
					if($UTR3>0){
						$gl=$gl+($UTR3[0]-$gene[0])+($exon[0]-$UTR3[$UTR3-1]-1);
						$u3l=$UTR3[1]-$UTR3[0]+1;
						for(my $i=2; $i<$UTR3;$i+=2){
							$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
							$g+=$UTR3[$i]-$UTR3[$i-1]-1;}					
					}else{  
						$gl+=$exon[0]-$gene[0];
					}
					$el=$exon[1]-$exon[0];		
					for(my $i=2; $i<$exon;$i+=2){
						$el+=$exon[$i+1]-$exon[$i]+1;
						$il+=$exon[$i]-$exon[$i-1]-1;
					}
					if($UTR5>0){
						$gl=$gl+($gene[1]-$UTR5[$UTR5-1])+($UTR5[0]-$exon[$exon-1]-1);
						$u5l=$UTR5[1]-$UTR5[0]+1;
						for(my $i=2; $i<$UTR5;$i+=2){
							$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
							$g+=$UTR5[$i]-$UTR5[$i-1]-1;}					
					}else{  
						$gl+=$gene[1]-$exon[$exon-1];
					}
										
					#info					
					for(my $i=0; $i<$UTR3;$i++){
						$info[$count]=$info[$count]."\t".$UTR3[$i];
					}	
					for(my $i=0; $i<$exon;$i++){
						$info[$count]=$info[$count]."\t".$exon[$i];
					}
					for(my $i=0; $i<$UTR5;$i++){
						$info[$count]=$info[$count]."\t".$UTR5[$i];
					}				
				}
				print OUT $info[$count]."\n";
				$fl=$gene[1]-$gene[0]+1;
			}
		}	
		elsif($tp[2] eq "gene"){
			$father=$tp[9];
		}
		elsif($tp[2] eq "mRNA"){
			#add new
			$d=$tp[6];
			$count++;	 
			@exon=();
			@UTR5=();
			@UTR3=();
			$el=0;
			$il=0;
			$gl=0;
			$u5l=0;
			$u3l=0;
			$canadd=0;
			$gene[0]=$tp[3];
			$gene[1]=$tp[4];
			for(my $i=0; $i<$tp;$i++){
				if($tp[$i] eq "ID"){
					$info[$count]="$tp[0]\t$tp[3]\t$tp[4]\t$tp[$i+1]\_$father\t$tp[6]";
					$name=$tp[$i+1];
					#print "$tp[$i+1] $name ";
					last;
				}
			}
			
		}
		else{  
			if($tp[2] eq "CDS" && index($lineIN,$name)>=0){	
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@exon,($tp[3],$tp[4]));
				}
			}elsif($tp[2] eq "five_prime_UTR" && index($lineIN,$name)>=0){
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@UTR5,($tp[3],$tp[4]));
					my $UTR5=@UTR5;
					print "$UTR5";
				}				
			}elsif($tp[2] eq "three_prime_UTR" && index($lineIN,$name)>=0){
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@UTR3,($tp[3],$tp[4]));
					my $UTR3=@UTR3;
					print "$UTR3";
				}				
			}
			else{
				if($canadd==1){
					$canadd=2;
				}	
			}
		}
	}
	print "\nTotal $count gene\n";

}


sub GetGff3GeneInfo_detail_our{
	my @gene;
	#0-scaf 1-start 2-end 3-ID 4-Direct 5-5UTR_No 6-Exon_No 7-3UTR_No 8-s1 9-e1 10-s2 11-e2 ...
	#if +
	#if -
	#8-9-10-11- change direct
	#gff3
	#scaffold619	.	gene	368802	373922	.	-	.	ID=LOC_Os01g37120.1.gene.63518;Name=LOC_Os01g37120.1.63518
	#Chr1	MSU_osa1r6	gene	1903	9817	.	+	.	ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
	#scaffold623	Genewise	gene	32921	35248	.	-	.	ID=LOC_Omer.gene.21524;Name=LOC_Omer.21524
	#scaffold623	Genewise	mRNA	32921	35248	.	-	.	ID=LOC_Omer.21524;Parent=LOC_Omer.gene.21524
	#scaffold623	Genewise	exon	35095	35248	.	-	.	ID=LOC_Omer.21524.exon1;Parent=LOC_Omer.21524
	#scaffold623	Genewise	CDS	35095	35248	.	-	0	ID=cds.LOC_Omer.21524;Parent=LOC_Omer.21524

	#if have alternative splicing, only use first record
	
	my $count=-1;
	my @gene;
	my @gene1;  #after output then add gene, so use gene1 hold
	my @info;
	my @infoD;
	my @exon;
	my @UTR5;
	my @UTR3;
	my $name;
	my $d;
	my $mRNAlen;
	my $CDSlen;
	my $UTR5len;
	my $UTR3len;
	my $intronlen;
	
	#count
	my $otd=$ARGV[0].".detail";
	my $otd1=$ARGV[0].".0st";
	open(OUTD,'>',"$otd");
	open(OUTD1,'>',"$otd1");
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		#if($tp[2] eq "gene" || $tp[2] eq "mRNA"){
		if($tp[2] eq "gene" || $tp[2] eq "pseudogene"){
			my @tpn= &Getheadinfo($lineIN, ";","=");
			#if($tpn[2] eq "mRNA"){				
			for(my $i=0; $i<$tp;$i++){
				if($tp[$i] eq "ID"){
					$name=$tpn[$i+1];
					#print "$tp[$i+1] $name ";
					last;
				}
			}
			#add gene
			$gene1[0]=$tp[3];	
			$gene1[1]=$tp[4];
			$gene1[2]=$gene1[1]-$gene1[0]+1;
		}
		if($tp[2] eq "mRNA"){
			#output
			print "+";
			if($count>=0){
				$intronlen=$mRNAlen-$CDSlen-$UTR5len-$UTR3len;
				$info[$count]="$mRNAlen\t$CDSlen\t$UTR5len\t$UTR3len\t$intronlen\t".$info[$count];
				$infoD[$count]="$mRNAlen\t$CDSlen\t$UTR5len\t$UTR3len\t$intronlen\t".$infoD[$count];
				$infoD1[$count]="$mRNAlen\t$CDSlen\t$UTR5len\t$UTR3len\t$intronlen\t".$infoD1[$count];
				#add exon info
				my $exon=@exon;
				my $exonno=$exon/2;
				my $UTR5=@UTR5;
				#print "$UTR5";
				my $UTR5no=0;
				if($UTR5>0){
					$UTR5no=$UTR5/2;
				}
				my $UTR3=@UTR3;
				#print "$UTR3";
				my $UTR3no=0;
				if($UTR3>0){
					$UTR3no=$UTR3/2;
				}
				$info[$count]=$info[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
				$infoD[$count]=$infoD[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
				$infoD1[$count]=$infoD1[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
				#print OUTD "$info[$count]\t";		
				if($d eq "+"){
					#info
					for(my $i=0; $i<$UTR5;$i+=2){
						$info[$count]=$info[$count]."\t".$UTR5[$i]."\t".$UTR5[$i+1];
						my $l1=$UTR5[$i]-$gene[0];
						my $l2=$UTR5[$i+1]-$gene[0];
						my $l=$l2-$l1+1;
						$infoD[$count]=$infoD[$count]."\t$l1\t$l2\t[$l]";
						$infoD1[$count]=$infoD1[$count]."\t$l1\t$l2";	
					}
					for(my $i=0; $i<$exon;$i+=2){
						$info[$count]=$info[$count]."\t".$exon[$i]."\t".$exon[$i+1];
						my $l1=$exon[$i]-$gene[0];
						my $l2=$exon[$i+1]-$gene[0];
						my $l=$l2-$l1+1;
						$infoD[$count]=$infoD[$count]."\t$l1\t$l2\t[$l]";
						$infoD1[$count]=$infoD1[$count]."\t$l1\t$l2";	
					}
					for(my $i=0; $i<$UTR3;$i+=2){
						$info[$count]=$info[$count]."\t".$UTR3[$i]."\t".$UTR3[$i+1];
						my $l1=$UTR3[$i]-$gene[0];
						my $l2=$UTR3[$i+1]-$gene[0];
						my $l=$l2-$l1+1;
						$infoD[$count]=$infoD[$count]."\t$l1\t$l2\t[$l]";	
						$infoD1[$count]=$infoD1[$count]."\t$l1\t$l2";						
					}
				}else{  #"-"
					#infoD
					#info
					for(my $i=$UTR3-2; $i>=0;$i-=2){
						my $l1=$UTR3[$i+1]-$gene[0];
						my $l2=$UTR3[$i]-$gene[0];
						my $l=$UTR3[$i+1]-$UTR3[$i]+1;
						$infoD[$count]=$infoD[$count]."\t$l2\t$l1\t[$l]";	
						$infoD1[$count]=$infoD1[$count]."\t$l2\t$l1";						
					}	
					for(my $i=$exon-2; $i>=0;$i-=2){
						my $l1=$exon[$i+1]-$gene[0];
						my $l2=$exon[$i]-$gene[0];
						my $l=$exon[$i+1]-$exon[$i]+1;
						$infoD[$count]=$infoD[$count]."\t$l2\t$l1\t[$l]";	
						$infoD1[$count]=$infoD1[$count]."\t$l2\t$l1";		
					}
					for(my $i=$UTR5-2; $i>=0;$i-=2){
						my $l1=$UTR5[$i+1]-$gene[0];
						my $l2=$UTR5[$i]-$gene[0];
						my $l=$UTR5[$i+1]-$UTR5[$i]+1;
						$infoD[$count]=$infoD[$count]."\t$l2\t$l1\t[$l]";	
						$infoD1[$count]=$infoD1[$count]."\t$l2\t$l1";	
					}				
					
					#info	
					@exon = sort{$a <=> $b} @exon;	
					@UTR3 = sort{$a <=> $b} @UTR3;
					@UTR5 = sort{$a <=> $b} @UTR5;				
					for(my $i=0; $i<$UTR3;$i+=2){
						$info[$count]=$info[$count]."\t".$UTR3[$i]."\t".$UTR3[$i+1];
					}	
					for(my $i=0; $i<$exon;$i+=2){
						$info[$count]=$info[$count]."\t".$exon[$i]."\t".$exon[$i+1];
					}
					for(my $i=0; $i<$UTR5;$i+=2){
						$info[$count]=$info[$count]."\t".$UTR5[$i]."\t".$UTR5[$i+1];
					}				
				}
				print OUT $info[$count]."\n";
				print OUTD $infoD[$count]."\n";
				print OUTD1 $infoD1[$count]."\n";
			}
			#add new
			$d=$tp[6];	$count++;	@exon=();	@UTR5=();	@UTR3=();	@info=();	@infoD=();	@infoD1=();
			$CDSlen = 0;	$UTR5len = 0;	$UTR3len = 0;	$intronlen = 0;	
			@gene=@gene1;
			$mRNAlen = $tp[4]-$tp[3]+1;
			
			#next line
			#$lineIN=<IN>;
			my @tpn= &Getheadinfo($lineIN, ";","=");
			#if($tpn[2] eq "mRNA"){	
				for(my $i=0; $i<$tp;$i++){
					if($tp[$i] eq "ID"){
						$info[$count]="$name\t$gene[0]\t$gene[1]\t$tp[0]\t$tp[3]\t$tp[4]\t$tp[$i+1]\t$tp[6]";
						my $l=$tp[4]-$tp[3]+1;
						$infoD[$count]="$name\t$gene[0]\t$gene[1]\t$tp[0]\t$l\t$tp[$i+1]\t$tp[6]";
						my $mst=$tp[3]-$gene[0];
						my $med=$tp[4]-$gene[0];
						$infoD1[$count]="$name\t$gene[0]\t$gene[1]\t$gene[2]\t$tp[0]\t$l\t$mst\t$med\t$tp[$i+1]\t$tp[6]";
						#$name=$tpn[$i+1];
						#print "$tp[$i+1] $name ";
						last;
					}
				}
			#}
			
		}else{  
			if($tp[2] eq "CDS"){	
			#if($tp[2] eq "CDS" && index($lineIN,$name)>=0){	
				push(@exon,($tp[3],$tp[4]));
				$CDSlen+=($tp[4]-$tp[3]+1);
			}elsif($tp[2] eq "five_prime_utr" || $tp[2] eq "five_prime_UTR"){
				push(@UTR5,($tp[3],$tp[4]));
				$UTR5len+=($tp[4]-$tp[3]+1);			
			}elsif($tp[2] eq "three_prime_utr" || $tp[2] eq "three_prime_UTR"){
				push(@UTR3,($tp[3],$tp[4]));
				$UTR3len+=($tp[4]-$tp[3]+1);			
			}
		}
	}
	
	#last out
	#output
	$intronlen=$mRNAlen-$CDSlen-$UTR5len-$UTR3len;
	$info[$count]="$mRNAlen\t$CDSlen\t$UTR5len\t$UTR3len\t$intronlen\t".$info[$count];
	$infoD[$count]="$mRNAlen\t$CDSlen\t$UTR5len\t$UTR3len\t$intronlen\t".$infoD[$count];
	$infoD1[$count]="$mRNAlen\t$CDSlen\t$UTR5len\t$UTR3len\t$intronlen\t".$infoD1[$count];
	print "+";

	my $exon=@exon;
	my $exonno=$exon/2;
	my $UTR5=@UTR5;
	#print "$UTR5";
	my $UTR5no=0;
	if($UTR5>0){
		$UTR5no=$UTR5/2;
	}
	my $UTR3=@UTR3;
	#print "$UTR3";
	my $UTR3no=0;
	if($UTR3>0){
		$UTR3no=$UTR3/2;
	}
	$info[$count]=$info[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
	$infoD[$count]=$infoD[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
	$infoD1[$count]=$infoD1[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
	#print OUTD "$info[$count]\t";		
	if($d eq "+"){
		#info
		for(my $i=0; $i<$UTR5;$i+=2){
			$info[$count]=$info[$count]."\t".$UTR5[$i]."\t".$UTR5[$i+1];
			my $l1=$UTR5[$i]-$gene[0];
			my $l2=$UTR5[$i+1]-$gene[0];
			my $l=$l2-$l1+1;
			$infoD[$count]=$infoD[$count]."\t$l1\t$l2\t[$l]";
			$infoD1[$count]=$infoD1[$count]."\t$l1\t$l2";	
		}
		for(my $i=0; $i<$exon;$i+=2){
			$info[$count]=$info[$count]."\t".$exon[$i]."\t".$exon[$i+1];
			my $l1=$exon[$i]-$gene[0];
			my $l2=$exon[$i+1]-$gene[0];
			my $l=$l2-$l1+1;
			$infoD[$count]=$infoD[$count]."\t$l1\t$l2\t[$l]";
			$infoD1[$count]=$infoD1[$count]."\t$l1\t$l2";	
		}
		for(my $i=0; $i<$UTR3;$i+=2){
			$info[$count]=$info[$count]."\t".$UTR3[$i]."\t".$UTR3[$i+1];
			my $l1=$UTR3[$i]-$gene[0];
			my $l2=$UTR3[$i+1]-$gene[0];
			my $l=$l2-$l1+1;
			$infoD[$count]=$infoD[$count]."\t$l1\t$l2\t[$l]";	
			$infoD1[$count]=$infoD1[$count]."\t$l1\t$l2";						
		}
	}else{  #"-"
		#infoD
		#info
		for(my $i=$UTR3-2; $i>=0;$i-=2){
			my $l1=$UTR3[$i+1]-$gene[0];
			my $l2=$UTR3[$i]-$gene[0];
			my $l=$UTR3[$i+1]-$UTR3[$i]+1;
			$infoD[$count]=$infoD[$count]."\t$l2\t$l1\t[$l]";	
			$infoD1[$count]=$infoD1[$count]."\t$l2\t$l1";						
		}	
		for(my $i=$exon-2; $i>=0;$i-=2){
			my $l1=$exon[$i+1]-$gene[0];
			my $l2=$exon[$i]-$gene[0];
			my $l=$exon[$i+1]-$exon[$i]+1;
			$infoD[$count]=$infoD[$count]."\t$l2\t$l1\t[$l]";	
			$infoD1[$count]=$infoD1[$count]."\t$l2\t$l1";		
		}
		for(my $i=$UTR5-2; $i>=0;$i-=2){
			my $l1=$UTR5[$i+1]-$gene[0];
			my $l2=$UTR5[$i]-$gene[0];
			my $l=$UTR5[$i+1]-$UTR5[$i]+1;
			$infoD[$count]=$infoD[$count]."\t$l2\t$l1\t[$l]";	
			$infoD1[$count]=$infoD1[$count]."\t$l2\t$l1";	
		}				
		
		#info	
		@exon = sort{$a <=> $b} @exon;	
		@UTR3 = sort{$a <=> $b} @UTR3;
		@UTR5 = sort{$a <=> $b} @UTR5;				
		for(my $i=0; $i<$UTR3;$i+=2){
			$info[$count]=$info[$count]."\t".$UTR3[$i]."\t".$UTR3[$i+1];
		}	
		for(my $i=0; $i<$exon;$i+=2){
			$info[$count]=$info[$count]."\t".$exon[$i]."\t".$exon[$i+1];
		}
		for(my $i=0; $i<$UTR5;$i+=2){
			$info[$count]=$info[$count]."\t".$UTR5[$i]."\t".$UTR5[$i+1];
		}				
	}
	print OUT $info[$count]."\n";
	print OUTD $infoD[$count]."\n";
	print OUTD1 $infoD1[$count]."\n";
	print "\nTotal $count gene\n";

	close(OUTD);
	close(OUTD1);
}


sub GetGff3GeneInfo_jv6{
	my @gene;
	#0-scaf 1-start 2-end 3-ID 4-Direct 5-5UTR_No 6-Exon_No 7-3UTR_No 8-s1 9-e1 10-s2 11-e2 ...
	#if +
	#if -
	#8-9-10-11- change direct
	#gff3
	#scaffold619	.	gene	368802	373922	.	-	.	ID=LOC_Os01g37120.1.gene.63518;Name=LOC_Os01g37120.1.63518
	#Chr1	MSU_osa1r6	gene	1903	9817	.	+	.	ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
	#if have alternative splicing, only use first record
	
	my $count=-1;
	my @gene;
	my @info;
	my @exon;
	my @UTR5;
	my @UTR3;
	my $name;
	my $d;
	#count
	my $fl;
	my $el;
	my $il;
	my $u5l;
	my $u3l;
	my $gl;
	my $canadd=2;
	my $otd=$ARGV[0].".detail";
	open(OUTD,'>',"$otd");
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN, ";","=");
		my $tp=@tp;
		if($tp[2] eq "gene"){
			#output
			print "+";
			if($count>=0){
				#add exon info
				my $exon=@exon;
				my $exonno=$exon/2;
				my $UTR5=@UTR5;
				#print "$UTR5";
				my $UTR5no=0;
				if($UTR5>0){
					$UTR5no=$UTR5/2;
				}
				my $UTR3=@UTR3;
				#print "$UTR3";
				my $UTR3no=0;
				if($UTR3>0){
					$UTR3no=$UTR3/2;
				}
				$info[$count]=$info[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;
				print OUTD "$info[$count]\t";		
				if($d eq "+"){
					#length
					if($UTR5>0){
						$u5il=$u5il+($UTR5[0]-$gene[0])+($exon[0]-$UTR5[$UTR5-1]-1);
						$u5l=$UTR5[1]-$UTR5[0]+1;
						for(my $i=2; $i<$UTR5;$i+=2){
							$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
							$u5il+=$UTR5[$i]-$UTR5[$i-1]-1;}					
					}else{  
						$u5il+=$exon[0]-$gene[0];
					}
					$el=$exon[1]-$exon[0]+1;	
					for(my $i=2; $i<$exon;$i+=2){
						$el+=$exon[$i+1]-$exon[$i]+1;
						$il+=$exon[$i]-$exon[$i-1]-1;
						
					}
					if($UTR3>0){
						$u3il=$u3il+($gene[1]-$UTR3[$UTR3-1])+($UTR3[0]-$exon[$exon-1]-1);
						$u3l=$UTR3[1]-$UTR3[0]+1;
						for(my $i=2; $i<$UTR3;$i+=2){
							$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
							$u3il+=$UTR3[$i]-$UTR3[$i-1]-1;}					
					}else{  
						$u3il+=$gene[1]-$exon[$exon-1];
					}
										
					#info
					for(my $i=0; $i<$UTR5;$i++){
						$info[$count]=$info[$count]."\t".$UTR5[$i];
					}
					for(my $i=0; $i<$exon;$i++){
						$info[$count]=$info[$count]."\t".$exon[$i];
					}
					for(my $i=0; $i<$UTR3;$i++){
						$info[$count]=$info[$count]."\t".$UTR3[$i];
					}
				}else{
					@exon = sort{$a <=> $b} @exon;	
					@UTR3 = sort{$a <=> $b} @UTR3;
					@UTR5 = sort{$a <=> $b} @UTR5;
					
					#length
					if($UTR3>0){
						$u3il=$u3il+($UTR3[0]-$gene[0])+($exon[0]-$UTR3[$UTR3-1]-1);
						$u3l=$UTR3[1]-$UTR3[0]+1;
						for(my $i=2; $i<$UTR3;$i+=2){
							$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
							$u3il+=$UTR3[$i]-$UTR3[$i-1]-1;}					
					}else{  
						$u3il+=$exon[0]-$gene[0];
					}
					$el=$exon[1]-$exon[0]+1;		
					for(my $i=2; $i<$exon;$i+=2){
						$el+=$exon[$i+1]-$exon[$i]+1;
						$il+=$exon[$i]-$exon[$i-1]-1;
					}
					if($UTR5>0){
						$u5il=$u5il+($gene[1]-$UTR5[$UTR5-1])+($UTR5[0]-$exon[$exon-1]-1);
						$u5l=$UTR5[1]-$UTR5[0]+1;
						for(my $i=2; $i<$UTR5;$i+=2){
							$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
							$u5il+=$UTR5[$i]-$UTR5[$i-1]-1;}					
					}else{  
						$u5il+=$gene[1]-$exon[$exon-1];
					}
										
					#info					
					for(my $i=0; $i<$UTR3;$i++){
						$info[$count]=$info[$count]."\t".$UTR3[$i];
					}	
					for(my $i=0; $i<$exon;$i++){
						$info[$count]=$info[$count]."\t".$exon[$i];
					}
					for(my $i=0; $i<$UTR5;$i++){
						$info[$count]=$info[$count]."\t".$UTR5[$i];
					}				
				}
				print OUT $info[$count]."\n";
				$fl=$gene[1]-$gene[0]+1;
				print OUTD "$fl\t$el\t$il\t$u5l\t$u3l\t$u5il\t$u3il\n";
			}
			#add new
			$d=$tp[6];
			$count++;	 
			@exon=();
			@UTR5=();
			@UTR3=();
			$el=0;
			$il=0;
			$u5l=0;
			$u3l=0;
			$u5il=0;
			$u3il=0;
			$canadd=0;
			$gene[0]=$tp[3];
			$gene[1]=$tp[4];
			#next line
			$lineIN=<IN>;
			my @tpn= &Getheadinfo($lineIN, ";","=");
			if($tpn[2] eq "mRNA"){	
				for(my $i=0; $i<$tp;$i++){
					if($tp[$i] eq "ID"){
						$info[$count]="$tp[0]\t$tp[3]\t$tp[4]\t$tp[$i+1]\t$tp[6]";
						$name=$tpn[$i+1];
						#print "$tp[$i+1] $name ";
						last;
					}
				}
			}
			
		}else{  
			if($tp[2] eq "CDS" && index($lineIN,$name)>=0){	
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@exon,($tp[3],$tp[4]));
				}
			}elsif($tp[2] eq "five_prime_UTR" && index($lineIN,$name)>=0){
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@UTR5,($tp[3],$tp[4]));
					my $UTR5=@UTR5;
					print "$UTR5";
				}				
			}elsif($tp[2] eq "three_prime_UTR" && index($lineIN,$name)>=0){
				if($canadd<2){
					$canadd=1;
				}
				if($canadd==1){
					push(@UTR3,($tp[3],$tp[4]));
					my $UTR3=@UTR3;
					print "$UTR3";
				}				
			}
			else{
				if($canadd==1){
					$canadd=2;
				}	
			}
		}
	}
	my $exon=@exon;
	my $exonno=$exon/2;
	my $UTR5=@UTR5;
	print "$UTR5";
	my $UTR5no=0;
	if($UTR5>0){
		$UTR5no=$UTR5/2;
	}
	my $UTR3=@UTR3;
	print "$UTR3";
	my $UTR3no=0;
	if($UTR3>0){
		$UTR3no=$UTR3/2;
	}
	$info[$count]=$info[$count]."\t".$UTR5no."\t".$exonno."\t".$UTR3no;	
	print OUTD "$info[$count]\t";	
	if($d eq "+"){
		#length
		if($UTR5>0){
			$u5il=$u5il+($UTR5[0]-$gene[0])+($exon[0]-$UTR5[$UTR5-1]-1);
			$u5l=$UTR5[1]-$UTR5[0]+1;
			for(my $i=2; $i<$UTR5;$i+=2){
				$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
				$u5il+=$UTR5[$i]-$UTR5[$i-1]-1;}					
		}else{  
			$u5il+=$exon[0]-$gene[0];
		}
		$el=$exon[1]-$exon[0]+1;		
		for(my $i=2; $i<$exon;$i+=2){
			$el+=$exon[$i+1]-$exon[$i]+1;
			$il+=$exon[$i]-$exon[$i-1]-1;
		}
		if($UTR3>0){
			$u3il=$u3il+($gene[1]-$UTR3[$UTR3-1])+($UTR3[0]-$exon[$exon-1]-1);
			$u3l=$UTR3[1]-$UTR3[0]+1;
			for(my $i=2; $i<$UTR3;$i+=2){
				$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
				$u3il+=$UTR3[$i]-$UTR3[$i-1]-1;}					
		}else{  
			$u3il+=$gene[1]-$exon[$exon-1];
		}
							
		#info
		for(my $i=0; $i<$UTR5;$i++){
			$info[$count]=$info[$count]."\t".$UTR5[$i];
		}
		for(my $i=0; $i<$exon;$i++){
			$info[$count]=$info[$count]."\t".$exon[$i];
		}
		for(my $i=0; $i<$UTR3;$i++){
			$info[$count]=$info[$count]."\t".$UTR3[$i];
		}
	}else{
		@exon = sort{$a <=> $b} @exon;	
		@UTR3 = sort{$a <=> $b} @UTR3;
		@UTR5 = sort{$a <=> $b} @UTR5;
		
		#length
		if($UTR3>0){
			$u3il=$u3il+($UTR3[0]-$gene[0])+($exon[0]-$UTR3[$UTR3-1]-1);
			$u3l=$UTR3[1]-$UTR3[0]+1;
			for(my $i=2; $i<$UTR3;$i+=2){
				$u3l+=$UTR3[$i+1]-$UTR3[$i]+1;
				$u3il+=$UTR3[$i]-$UTR3[$i-1]-1;}					
		}else{  
			$u3il+=$exon[0]-$gene[0];
		}
		$el=$exon[1]-$exon[0]+1;		
		for(my $i=2; $i<$exon;$i+=2){
			$el+=$exon[$i+1]-$exon[$i]+1;
			$il+=$exon[$i]-$exon[$i-1]-1;
		}
		if($UTR5>0){
			$u5il=$u5il+($gene[1]-$UTR5[$UTR5-1])+($UTR5[0]-$exon[$exon-1]-1);
			$u5l=$UTR5[1]-$UTR5[0]+1;
			for(my $i=2; $i<$UTR5;$i+=2){
				$u5l+=$UTR5[$i+1]-$UTR5[$i]+1;
				$u5il+=$UTR5[$i]-$UTR5[$i-1]-1;}					
		}else{  
			$u5il+=$gene[1]-$exon[$exon-1];
		}
							
		#info					
		for(my $i=0; $i<$UTR3;$i++){
			$info[$count]=$info[$count]."\t".$UTR3[$i];
		}	
		for(my $i=0; $i<$exon;$i++){
			$info[$count]=$info[$count]."\t".$exon[$i];
		}
		for(my $i=0; $i<$UTR5;$i++){
			$info[$count]=$info[$count]."\t".$UTR5[$i];
		}				
	}
	print OUT $info[$count]."\n";
	$fl=$gene[1]-$gene[0]+1;
	print OUTD "$fl\t$el\t$il\t$u5l\t$u3l\t$u5il\t$u3il\n";
	$count++;
	print "\nTotal $count gene\n";

	close(OUTD);
}


sub RegulateBESctg{
	#################################################################
	#
	#  Data analysys have not done
	#
	#################################################################
	
	#ctg99	Contig93469	1	339	Contig163595	1	432
	#put all ctg info into hash,name as key,if has same name within same region, combin them
	#then sort by loc in value,combin same name record, output
	my $tpctg="start";
	my %ctghash;
	my %sameloc;
	my %directmiss;
	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);	
		if($tpctg ne $tp[0]){
			if($tpctg ne "start"){
				my %lochash;
				my $locinfo;
				my $name;
				#combin loc
				while (($name, $locinfo) = each(%ctghash)){
					my @tploc=&Getheadinfo($locinfo);
					if ($lochash{$tploc[1]}){
						$lochash{$tploc[1]}=$lochash{$tploc[1]}."|".$name.",".$tploc[0];
						$sameloc{$tploc[1]}=$lochash{$tploc[1]};	
					}else{
						$lochash{$tploc[1]}=$name.",".$tploc[0];
					}
					delete($ctghash{$name});
				}
				#output
				my @loc = sort{$a <=> $b} keys(%lochash);
				my $loc=@loc;
				for(my $i=0; $i<$loc;$i++){
					$d=($loc[$i+1]-$loc[$i])*100;
					if($i != ($loc-1)){
						print OUT $lochash{$loc[$i]}.",".$d."\t";	
					}else{
						print OUT $lochash{$loc[$i]}.",-1\n";
					}
					delete($lochash{$loc[$i]});
				}
			}
			$tpctg = $tp[0]; 
		}
		#add record, same name 1,within same region-combin 2,contradiction-output
		#first group
		if($ctghash{$tp[1]}){
			my @tploc=&Getheadinfo($ctghash{$tp[1]});
			if($tploc[0]==$tp[2]){
				#get mid loc
				$tploc[1]=($tploc[1]+$tp[3])/2;
				$ctghash{$tp[1]}=$tploc[0]." ".$tploc[1];
			}else{
				#output to contradiction list
				$directmiss{$tp[1]}=2;
				$tploc[1]=($tploc[1]+$tp[3])/2;
				$ctghash{$tp[1]}="0 ".$tploc[1];
			}
		}else{
			$ctghash{$tp[1]}=$tp[2]." ".$tp[3];	
		}
		#second group
		if($ctghash{$tp[4]}){
			my @tploc=&Getheadinfo($ctghash{$tp[4]});
			if($tploc[0]==$tp[5]){
				#get mid loc
				$tploc[1]=($tploc[1]+$tp[6])/2;
				$ctghash{$tp[4]}=$tploc[0]." ".$tploc[1];
			}else{
				#output to contradiction list
				$directmiss{$tp[4]}=2;
				$tploc[1]=($tploc[1]+$tp[6])/2;
				$ctghash{$tp[4]}="0 ".$tploc[1];
			}
		}else{
			$ctghash{$tp[4]}=$tp[5]." ".$tp[6];	
		}		
		
	}
	#output last record
	my %lochash;
	my $locinfo;
	my $name;
	#combin loc
	while (($name, $locinfo) = each(%ctghash)){
		@tploc=&Getheadinfo($locinfo);
		if ($lochash{$tploc[1]}){
			$lochash{$tploc[1]}=$lochash{$tploc[1]}."|".$name.",".$tploc[0];	
		}else{
			$lochash{$tploc[1]}=$name.",".$tploc[0];
		}
	}
	#output
	my @loc = sort keys(%lochash);
	my $loc=@loc;
	for(my $i=0; $i<$loc;$i++){
		$d=($loc[$i+1]-$loc[$i])*100;
		if($i != ($loc-1)){
			print OUT $lochash{$loc[$i]}.",".$d."\t";	
		}else{
			print OUT $lochash{$loc[$i]}.",-1\n";
		}
		delete($lochash{$loc[$i]});
	}
	#output same loc info
	print OUT "name in Same Loc-----------------------------\n";
	@loc = sort keys(%sameloc);
	$loc=@loc;
	for(my $i=0; $i<$loc;$i++){
		print OUT $loc[$i]."\t".$sameloc{$loc[$i]}."\n";	
		delete($sameloc{$loc[$i]});
	}	
	#output direction contradiction info
	print OUT "Direct error info-----------------------------\n";
	@loc = sort keys(%directmiss);
	$loc=@loc;
	for(my $i=0; $i<$loc;$i++){
		print OUT $loc[$i]."\t".$directmiss{$loc[$i]}."\n";	
		delete($directmiss{$loc[$i]});
	}
}

sub REannotatedCount{
	#remove head
	#input REannotated file
	#1  377   24.1  6.8  6.8  chr01_13101     70889    71079 (43197800) C RL004_043-int      LTR/Copia          (19462)    932     742      t10111
	#1	0		1	2	 3		4				5		6		7		8		9				10					11		12		13		14
	
	#output 
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent
	#lines which begin with ## -- drop
	#different type save to different file
	my @rm; 
	my $count=0; 
	my @type;
	my $typecount=0;
	my @lasttp;
	my @tp;
	
	#read in REannotated file
	while($lineIN=<IN>){
		chomp($lineIN);
     	@tp= &Getheadinfo($lineIN,"/");	
		if(index($lineIN,'/')<0){
			$tp[11]="z-other";
		}
		
		#test overlap
		#0-type 1-family 2-name 3-scafold 4-direct 5-start 6-end 7-score
		if($count>0){
			if($tp[5]<$lasttp[6] && $tp[4] eq $lasttp[3] ){
				
				if($tp[6]<$lasttp[6]){
					print "+";
					next;								
				}else{
					print ".";
					if($tp[11] eq "z-other" || $tp[10] eq "Simple_repeat" || $tp[10] eq "Satellite" || $tp[10] eq "Unknown"){
						$tp[5]=	$lasttp[6]+1;
					}elsif($lasttp[2] eq "z-other" || $tp[10] eq "Simple_repeat" || $tp[10] eq "Satellite" || $tp[10] eq "Unknown"){
						$lasttp[6] = $tp[5] -1;
					}else{
						if($tp[0]>$lasttp[7]){
							$lasttp[6] = $tp[5] -1;	
							
						}else{
							$tp[5]=	$lasttp[6]+1;
						}
					}
				}
			}
			$rm[$count-1]="$lasttp[0]\t$lasttp[1]\t$lasttp[2]\t$lasttp[3]\t$lasttp[4]\t$lasttp[5]\t$lasttp[6]\t$lasttp[7]\n";
		}
		@lasttp=($tp[10],$tp[11],$tp[9],$tp[4],$tp[8],$tp[5],$tp[6],$tp[0]);
		
		$count++;
		my $canadd=1;
		for(my $i=0; $i<$typecount;$i++){
			if($type[$i] eq $tpf[0]){
				$canadd=0;
				last;
			}
		}
		if($canadd == 1){
			$type[$typecount]=$tpf[0];
			$typecount++;
		}

	}
	$rm[$count-1]="$lasttp[0]\t$lasttp[1]\t$lasttp[2]\t$lasttp[3]\t$lasttp[4]\t$lasttp[5]\t$lasttp[6]\t$lasttp[7]\n";
	print "\n";
	
	#classify result
	#0-type 1-family 2-name 3-scafold 4--direct 5-start 6-end 7-score
	@type = sort(@type);
	@rm=sort(@rm);
	#percent count
	my @p00; #total
	my $tptype="99999_start_symble";
	my $tpfamily;
	my @t00; #total
	$p00[0]=0;
	$p00[1]=0;
	$t00[0]=0;
	$t00[1]=0;	
	my @stype;
		
	for(my $i=0; $i<$count;$i++){
		my @tp= &Getheadinfo($rm[$i]);
		if($tptype ne $tp[0]){
			#new type
			if($tptype ne "99999_start_symble"){
				#output family count
				print OUT "\t$tpfamily\t$t00[0]\t$t00[1]\n";
				#output type count
				print OUT "\tCount\t$p00[0]\t$p00[1]\n";
				#clean old type
				close(OUTD);
				$p00[0]=0;
				$p00[1]=0;
				$t00[0]=0;
				$t00[1]=0;
			}
			
			#open new type
			$tptype=$tp[0];
			$tpfamily=$tp[1];		
			print OUT "$tptype\n";
			$ot=$ARGV[0].".".$tptype.".detail";
			open(OUTD,'>',"$ot");
			print OUTD "0-type 1-family 2-name 3-scafold 4-direct 5-start 6-end 7-score\n";		
		}else{			
			if($tpfamily ne $tp[1]){
				#output old family
				print OUT "\t$tpfamily\t$t00[0]\t$t00[1]\n";
				#new family
				$tpfamily=$tp[1];
				$t00[0]=0;
				$t00[1]=0;		
			}
		}
		#add info
		print OUTD $rm[$i]."\n";
		$stype[$i]="$tp[3]\t$tp[5]\t$tp[6]\t$tp[0]#$tp[1]/$tp[2]\t$tp[4]";
		my $len=$tp[6]-$tp[5]+1;
		$p00[0]++;
		$p00[1]=$p00[1]+$len;	
		$t00[0]++;
		$t00[1]=$t00[1]+$len;		
	}
	print OUT "\t$tpfamily\t$t00[0]\t$t00[1]\n";
	#output type count
	print OUT "\tCount\t$p00[0]\t$p00[1]\n";
	#clean old type
	
	
	my $otall=$ARGV[0].".all.detail";
	open(OUTA,'>',"$otall");
	@stype = sort(@stype);
	#0-scaf 1-start 2-end 3-type#family/name 4-direct
	for(my $i=0; $i<$count;$i++){
		print OUTA "$stype[$i]\n";	
	}
	close(OUTA);

}
sub RepeatMaskerCat_2_out{
	#cat
	#270 19.45 2.70 2.70 chromosome01 1242 1315 (45063454) rnd-6_family-1695#RC/Helitron 1049 1122 (902) m_b1s001i2
	#225 15.59 12.50 9.41 chromosome01 1247 1331 (45063438) C DTM_#DNAt_MULE_japo_Os2088:Chr2:17937451:17941272 (3455) 367 280 m_b1s001i3
	#out
	#   493    5.1  2.0  9.3  Sc0000000        6     113 (7596657) + DTC_                             DNA/En-Spm_rufi_CL19C_fasta_Contig6                     2035   2134  (2329)      1 *
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~m/(?:^)[0-9].*/){
			my @tp= &Getheadinfo($lineIN);
			my @name=();
			if($tp[8] eq "C"){
				@name = split(/\#/,$tp[9]);	
				print OUT "$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]\tC\t$name[0]\t$name[1]\t$tp[10]\t$tp[11]\t$tp[12]\t$tp[13]\n";
			}else{
				@name = split(/\#/,$tp[8]);	
				print OUT "$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]\t$tp[4]\t$tp[5]\t$tp[6]\t$tp[7]\t+\t$name[0]\t$name[1]\t$tp[9]\t$tp[10]\t$tp[11]\t$tp[12]\n";
			}	
			
		}
	}
}


sub RepeatMaskerCatDetailCount{
	#only test length for family, do not count actually total length
	#remove head
	#825 19.71 13.94 0.00 scaff_2259 277519 277726 (5553466) rnd-5_family-7520#RC/Helitron 1 237 (80) 5
	#1062 12.57 0.00 1.51 scaff_2259 277738 277939 (5553253) C rnd-2_family-228#DNA (5) 201 3 5
	
	@rm; #two-dimensional array
	$count=0; 
	@type;
	$typecount=0;
	#read in cat file
	while($lineIN=<IN>){
		chomp($lineIN);
		if(substr($lineIN, 0, 1) ne "#"){
     		my @tp= &Getheadinfo($lineIN);
			my $tp=@tp;
			my $len=$tp[6]-$tp[5];      
			my @name;
			my $Tlength;
			my $Tfrom;
			my $Tto;
			if ($tp[8] ne "C"){
				@name=&Getheadinfo($tp[8],"#","/");
				my $remain=substr(substr($tp[11],1,(length($tp[11])-1)),0,(length($tp[11])-2));
				$Tlength=$tp[10]+$remain;
				$Tfrom=$tp[9];
				$Tto=$tp[10];
			}else{
				#complement
				@name=&Getheadinfo($tp[9],"#","/");
				my $remain=substr(substr($tp[10],1,(length($tp[10])-1)),0,(length($tp[10])-2));
				$Tlength=$remain+$tp[11];
				$Tfrom=$tp[12];
				$Tto=$tp[11];				
			}	
			my $name=@name;
			if($name<3){
				for(my $j=$name;$j<3;$j++){
					$name[$j]="z_other";
				}
			}
			if($_[1] ne ""){
				my @tptp=&Getheadinfo($name[0],$_[1]);
				$name[0]=$tptp[0];	
			}
			
			$rm[$count]="$name[1]\t$name[2]\t$name[0]\t$tp[4]\t$tp[5]\t$tp[6]\t$len\t$Tfrom\t$Tto\t$Tlength\t$percent\n";
			$count++;
			my $canadd=1;
			for(my $i=0; $i<$typecount;$i++){
				if($type[$i] eq $name[1]){
					$canadd=0;
					last;
				}
			}
			if($canadd == 1){
				$type[$typecount]=$name[1];
				$typecount++;
			}
		}else{
			last;
		}
	}
	#classify result
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent
	@type = sort(@type);
	@rm=sort(@rm);

	my $tptype="99999_start_symble";
	my $tpfamily;
	my $tpname;
	my @p00; #type
	my @t00; #family
	my @n00; #name
		
	for(my $i=0; $i<$count;$i++){
		my @tp= &Getheadinfo($rm[$i]);
		if($tptype ne $tp[0]){
			#new type
			if($tptype ne "99999_start_symble"){
				#output family count
				print OUT "\t\t$tpname\t$n00[0]\t$n00[1]\n";
				#output family count
				print OUT "\t$tpfamily\t$t00[0]\t$t00[1]\n";
				#output type count
				print OUT "\tCount\t$p00[0]\t$p00[1]\n";
				#clean old type
				$p00[0]=0;
				$p00[1]=0;
				$t00[0]=0;
				$t00[1]=0;
				$n00[0]=0;
				$n00[1]=0;
			}
			
			#open new type
			$tptype=$tp[0];
			$tpfamily=$tp[1];
			$tpname=$tp[2];		
			print OUT "$tptype\n";
					
		}else{			
			if($tpfamily ne $tp[1]){
				#output family count
				print OUT "\t\t$tpname\t$n00[0]\t$n00[1]\n";
				#output old family
				print OUT "\t$tpfamily\t$t00[0]\t$t00[1]\n";
				#new family
				$tpfamily=$tp[1];
				$tpname=$tp[2];	
				$t00[0]=0;
				$t00[1]=0;
				$n00[0]=0;
				$n00[1]=0;	
			}else{
				if($tpname ne $tp[2]){
					#output family count
					print OUT "\t\t$tpname\t$n00[0]\t$n00[1]\n";	
					$tpname=$tp[2];	
					$n00[0]=0;
					$n00[1]=0;	
				}
			}
		}
		#add info
		$p00[0]++;
		$p00[1]=$p00[1]+$tp[6];	
		$t00[0]++;
		$t00[1]=$t00[1]+$tp[6];		
		$n00[0]++;
		$n00[1]=$n00[1]+$tp[6];	
			
	}
	print OUT "\t\t$tpname\t$n00[0]\t$n00[1]\n";
	print OUT "\t$tpfamily\t$t00[0]\t$t00[1]\n";
	#output type count
	print OUT "\tCount\t$p00[0]\t$p00[1]\n";

}


sub RepeatMaskerCatOUTtypetotalCombin{
	# sort by chrom and loc, add type
	#same family in region , combin
	#LTR	z_other	RM009_049	chr01	817122	818195	1074	2075	3172	3172	33.858764186633	-
	#LTR	z_other	RM009_049	chr01	818243	820273	2031	1	2074	3172	64.0290037831021	-
	#0-class 1-family 2-type 3-chrom 4-st 5-ed 6-len 7-Tst 8-Ted 9-Tlen 10-percent 11-direct 
	#out put same list
	my $range=0.5;   
	if($_[1]>0){
		$range=$_[1];
	}
	my $count=0;
	my @info;
	
	while($lineIN=<IN>){
		chomp($lineIN);
     		my @tp= &Getheadinfo($lineIN);
     		#my @name=$tp[2]."#".$tp[0]."/".$tp[1];
     		#$tp[0]=$name[0];
     		my $l=@info;
     		my $add=1;
     		#another chromsome	
		if($tp[0] eq $info[0] && $tp[1] eq $info[1] && $tp[3] eq $info[3]){
			my $d=$range * $tp[9];
			if($_[1]>0){
				$d=$_[1];
			}
			if(($tp[4]-$d) < $info[5]){
				#combin
				print ".";
				$add=0;
				$tp[2] = $tp[2]."#".$info[2];
				$tp[5]=$tp[5]>=$info[5]?$tp[5]:$info[5];
				$tp[9]=$tp[9]>=$info[9]?$tp[9]:$info[9];
				$tp[6]=$tp[5]-$tp[4]+1;
				$tp[10]=$tp[6]/$info[9];
			}
	
		}
		if($add==1 && $count>0){
			print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t$info[9]\t$info[10]\t$info[11]\n";			
		}
		@info=@tp;
		$count++;			
     	}
     	print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t$info[9]\t$info[10]\t$info[11]\n";
     	print "\n";
}


sub RepeatMaskerOutReadsType{
	#for each reads, find the best match
	
	#Tls_697 use full name, other use family name
	# 1163    4.0  0.0  0.7  HWI-D00105:203:H0R0FADXX:1:1101:10000:19529      1   151    (0) + TL001_3364        Unspecified         10000  10149  (1053)      1  
	#type
	#DNA/ LINE/ SINE/
	
	
	my %typenumber; 
	my @lastline;
	
	$lineIN=<IN>;
	$lineIN=<IN>;
	$lineIN=<IN>;
			
	while($lineIN=<IN>){
		chomp($lineIN);
     	my @tp= &Getheadinfo($lineIN);

		###### only for tea#################
		if($tp[9]=~/TL/){
			$tp[9] = substr($tp[9],0,index($tp[9],"_",));			
		}
		####################################
		elsif($tp[9]=~/Tls/){
			
		}else{
			#otehr use $tp[10]
			$tp[9] = $tp[10];
		}
		
		if($lastline[4] ne $tp[4]){
			if(length($lastline[4]) > 0){
				#add
				if($typenumber{$lastline[9]}>0){
					$typenumber{$lastline[9]}=$typenumber{$lastline[9]}+1;
				}else{
					$typenumber{$lastline[9]}=1;
				}
			}
			@lastline=@tp;
		}else{
			#leave lines with biggest SW score 
			if($lastline[0] < $tp[0]){
				@lastline=@tp;
			}
		}			
	}
	#add
	if($typenumber{$lastline[9]}>0){
		$typenumber{$lastline[9]}++;
	}else{
		$typenumber{$lastline[9]}=1;
	}	
		
	#print out
	my @lineskey=keys(%typenumber);
	my $count=@lineskey;
	for(my $i=0;$i<$count;$i++){
		print OUT "$lineskey[$i]\t$typenumber{$lineskey[$i]}\n";
	}
    print "\n";
}

sub RepeatMaskerOutPercent{

	#4-chrom 5-Qst 6-Qed 7-Qleft 8-direct 9-name 10-type 11-13 Tst Ted Tleft
	#   25    5.1  0.0  0.0  C6150428             1      39      (61) + AT_rich    Low_complexity      1     39  (261)      1  
	#  293   26.1  0.0  0.0  C6150458             3      94       (6) C J7L_1046_3 Unknown         (810)   2218   2127      2  

	@db; #0-sortkey 1-name 2-chrom 3-st 4-ed 5-Qlen  6-chrom_len 7-Tst 8-Ted  9-Tleft 10-Tlen 11-direct 12-percent
	@data; #2-d
	print "read in data...\n";
	my $count=0;
	my $countd=0;
	my $max=10000000000;
	if($_[1] ne "n"){
		$lineIN=<IN>;
		$lineIN=<IN>;
		$lineIN=<IN>;
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);
		$lineIN=~s/[\(\)]//g;
     	my @tp= &Getheadinfo($lineIN);
     	@db=();
     	if($tp[10] ne "Low_complexity" && $tp[10] ne "Simple_repeat"){
     		my $tpadd=$max+$tp[5]; 
     		$db[0]=$tp[4]."_$tpadd";#sort key
     		$db[1]=$tp[9];			#name
     		$db[2]=$tp[4];			#chrom
     		$db[3]=$tp[5];			#Qst
     		$db[4]=$tp[6];			#Qed
     		$db[5]=$tp[6]-$tp[5]+1;	#Qlen
     		$db[6]=$tp[6]+$tp[7];	#chrom_len
        	  	
     		if($tp[8] eq '+'){
     			$db[7]=$tp[11];				#Tst
     			$db[8]=$tp[12];				#Ted
     			$db[9]=$tp[13];				#Tleft
     			$db[10]=$tp[12]+$tp[13];	#Tlen
     			$db[11]="+";
     		}else{
     			$db[7]=$tp[13];				#Tst
     			$db[8]=$tp[12];				#Ted
     			$db[9]=$tp[11];				#Tleft
     			$db[10]=$tp[12]+$tp[11];	#Tlen    
     			$db[11]="-";		
     		}
     		$db[12]=$db[5]/$db[10];
     		
     		#cut $chrom_len << Tlen 80%
     		#cut $percent <<80%
     		if($db[6]/$db[10]>=0.8 && $db[12]>=0.8){
     			$data[$count]="$db[0]\t$db[1]\t$db[2]\t$db[3]\t$db[4]\t$db[5]\t$db[6]\t$db[7]\t$db[8]\t$db[9]\t$db[10]\t$db[11]\t$db[12]";
     			$count++;
				#print "+";	
			}else{
				#print "-";
				$countd++;
			}	
		}
		
    }
   	@data=sort(@data);
   	@db=();
    $data=@data;
	my $om=$ARGV[0].".mid";
	open(OUTM,'>',"$om"); 
	for(my $j=0;$j<$data;$j++){
		print OUTM "$data[$j]\n";
	}
	print "\n";
    close(OUTM);
    print "\nget $count\tremove $countd\n";
    
   	#remove overlap

   	##loop compair leave big one (if overrlap < 10% all leave)
    #loop clean
	my %o;  #count overlap
	my $ot=$ARGV[0].".detail";
	open(OUTD,'>',"$ot");  
    my $d=9999;
    my $loop=0;
    my $dc=0;
    #0-sortkey 1-name 2-chrom 3-st 4-ed 5-Qlen  6-chrom_len 7-Tst 8-Ted  9-Tleft 10-Tlen 11-direct 12-percent
    while($d>0){
    	$d=0;
    	$dc=0;
    	@db=@data;
    	@data=();
    	$db=$tp;
    	$loop++;
    	my $db=@db;    	
    	print "\nLoop $loop: $db\t";
    	print OUTD "##Loop $loop: $tp-----------------------\n";
   		for(my $i=1;$i<$db;$i++){
   			#print "$data[$i][0]\n";
   			@db1=&Getheadinfo($db[$i-1]);
   			@db2=&Getheadinfo($db[$i]);
    		if($db2[2] eq $db1[2]){	#same chrom
    			my $ov=$db1[4]-$db2[3];
    			if($ov>20){   #overlap #choice big one 
    				$d++; 
    				#print ".";	
    				if($db1[5]>$db2[5]){
    					print OUTD "$ov\t$db2[0],$db2[1],$db2[2],$db2[3],$db2[4],$db2[5],$db2[6],$db2[7],$db2[8],$db2[9],$db2[10],$db2[11],$db2[12] <- $db1[0],$db1[1],$db1[2],$db1[3],$db1[4],$db1[5],$db1[6],$db1[7],$db1[8],$db1[9],$db1[10],$db1[11],$db1[12]\n";				
						$db[$i]=$db[$i-1] 
												
    				}else{
    					print OUTD "$ov\t$db2[0],$db2[1],$db2[2],$db2[3],$db2[4],$db2[5],$db2[6],$db2[7],$db2[8],$db2[9],$db2[10],$db2[11],$db2[12] -> $db1[0],$db1[1],$db1[2],$db1[3],$db1[4],$db1[5],$db1[6],$db1[7],$db1[8],$db1[9],$db1[10],$db1[11],$db1[12]\n";
    				}

    			}else{
    				if($ov>0){
    					#print "-$ov "; #seperate
    				}
    				#print ".";
    				#$data[$dc]=[$db1];
    				$data[$dc]="$db1[0]\t$db1[1]\t$db1[2]\t$db1[3]\t$db1[4]\t$db1[5]\t$db1[6]\t$db1[7]\t$db1[8]\t$db1[9]\t$db1[10]\t$db1[11]\t$db1[12]";
    				#print "$data[$dc][0]\n";
    				$dc++;
    			}
    		}			
		
		}
	}
	print "\nremain $dc\n";
	my $ol=$ARGV[0].".loc";
	open(OUTL,'>',"$ol"); 	
	#0-sortkey 1-name 2-chrom 3-st 4-ed 5-Qlen  6-chrom_len 7-Tst 8-Ted  9-Tleft 10-Tlen 11-direct 12-percent
	for(my $j=0;$j<$dc;$j++){
		@tp=&Getheadinfo($data[$j]);
		$tp=@tp;
		my $newname=$tp[1]."=$j";
		my $st=$tp[3]-$tp[7];
		if($st<0){
			$st=0;
		}
		my $ed=$tp[4]+$tp[9];
		print OUTL ">$newname|$tp[10]|$tp[2]:$st..$ed\n";
		print OUT "$newname";
		for(my $i=1;$i<$tp;$i++){
			print OUT "\t$tp[$i]";
		}
		print OUT "\n";
	}
	print "\n";
	close(OUTD);
	close(OUTL);
}



sub RepeatMaskerOUTTypeLoopCombin{

	#use first value to cut
	#out
	#   897   10.3  7.3  3.2  scaffold1           1     151 (9561602) + rnd-3_family-127   Unknown               392    548   (382)      1 *
	#   596   25.4  1.1  0.6  scaffold1       26486   26663 (9535090) C rnd-5_family-1105  Unknown            (1967)    940     762     33 *
	#0-value 1-2-3- 4-chrom 5-st 6-ed 7- 8-c/+ 9-name 10-type
	#cat
	#677 18.65 1.32 5.50 scaffold39 9935 10161 (12274) C Single160 (1647) 1033 816 5
	#399 17.53 3.06 1.00 scaffold39 10078 10175 (12260) Single160 984 1083 (1597) 5	
	#0-value 1-2-3- 4-chrom 5-st 6-ed 7- 8-c/+name 9-name
	#DNA/RC/LINE/SINE score add 0 behand

	# @db; save data
	print "read in data...\n";
	my $count=0;
	my @row;
	my $max=100000000000000;
	
		#$lineIN=<IN>;
		#$lineIN=<IN>;
		#$lineIN=<IN>;

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
}

sub RepeatMaskerCatOUTTypeLoopCombin{

	#use first value to cut
	#677 18.65 1.32 5.50 scaffold39 9935 10161 (12274) C Single160 (1647) 1033 816 5
	#399 17.53 3.06 1.00 scaffold39 10078 10175 (12260) Single160 984 1083 (1597) 5	
	#0-value 1-2-3- 4-chrom 5-st 6-ed 7- 8-c/+name 9-name
	#DNA/RC/LINE/SINE score add 0 behand

	# @db; save data
	print "read in data...\n";
	my $count=0;
	my @row;
	my $max=100000000000000;
	
	if($_[1] eq "out"){
		$lineIN=<IN>;
		$lineIN=<IN>;
		$lineIN=<IN>;
	}else{
		#cat file
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);
     	my @tp= &Getheadinfo($lineIN);
     	if($tp[8] eq 'C' || $tp[8] eq '+'){
     		$tp[8]=$tp[9];
     	}
     	my @tpsb=&Getheadinfo($tp[8],"_");
     	my $tpsb=@tpsb;
     	my $name=$tpsb[0];
     	if($tpsb==1){
     		#$type="other";
     	}
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
		

		$row[$count]="$tp[4]-$tpadd\t$tp[4]\t$tp[5]\t$tp[6]\t$len\t$tp[10]\t$name\t$tp[0]\n";
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
    			$tpdb[$tpc]="$db[$i-1][1]-$tpadd\t$db[$i-1][1]\t$db[$i-1][2]\t$db[$i-1][3]\t$db[$i-1][4]\t$db[$i-1][5]\t$db[$i-1][6]\t$db[$i-1][7]\n"; 	    			
    			#print "length $l last line $db[$i-1][1]-$tpadd\n";
    			$tpc++;  			
    		}
    		if($add ==2){
				my $tpadd=$max+$db[$i-1][2]; 	
    			$tpdb[$tpc]="$db[$i-1][1]-$tpadd\t$db[$i-1][1]\t$db[$i-1][2]\t$db[$i-1][3]\t$db[$i-1][4]\t$db[$i-1][5]\t$db[$i-1][6]\t$db[$i-1][7]\n"; 	    			
    			$tpc++;
 				$tpadd=$max+$db[$i][2]; 	
    			$tpdb[$tpc]="$db[$i][1]-$tpadd\t$db[$i][1]\t$db[$i][2]\t$db[$i][3]\t$db[$i][4]\t$db[$i][5]\t$db[$i][6]\t$db[$i][7]\n"; 	    			
    			$tpc++;
    			$i++;
    		}
    	} 
		my $tpadd=$max+$db[$db-1][2]; 	
    	$tpdb[$tpc]="$db[$db-1][1]-$tpadd\t$db[$db-1][1]\t$db[$db-1][2]\t$db[$db-1][3]\t$db[$db-1][4]\t$db[$db-1][5]\t$db[$db-1][6]\t$db[$db-1][7]\n"; 	    			
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
}

sub RepeatMaskerCatOUTTypeSimpleCombin{
	#same family simple combin
	#677 18.65 1.32 5.50 scaffold39 9935 10161 (12274) C Single160 (1647) 1033 816 5
	#399 17.53 3.06 1.00 scaffold39 10078 10175 (12260) Single160 984 1083 (1597) 5
	#
	#LTR	z_other	RM009_049	chr01	817122	818195	1074	2075	3172	3172	33.858764186633	-
	#0-1-2-3- 4-chrom 5-st 6-ed 7- 8-c/name 9-name
	#mid 0-name 1-len 2-3- 4-chrom 5-st 6-ed 7-c/name 8-name
	#out put same list
	#search before 10 

	if($_[1] eq "out"){
		$lineIN=<IN>;
		$lineIN=<IN>;
		$lineIN=<IN>;
	}else{
		#cat file
	}
	
	#sample,all cut first one
	my $count=0;
	my $countD=0;
	my $countA=0;
	my @info;
	my %o;  #count overlap
	my $c=0;
	my $otp=$ARGV[0].".tp";
	my @cb; #0-chrom 1-st 2-ed
	my $cbc=-1;
	open(OUTT,'>',"$otp");
	while($lineIN=<IN>){
		chomp($lineIN);
     	my @tp= &Getheadinfo($lineIN);
     	if($tp[8] eq 'C' || $tp[8] eq '+'){
     		$tp[8]=$tp[9];
     	}
     	my $add=1;
     	my @tpsb=&Getheadinfo($tp[8],"_");
     	$tp[0]=$tpsb[0];
     	$tp[1]=$tp[6]-$tp[5]+1;
     	$countA+=$tp[1];
     	#another chromsome	
     	my $cbca=1;
		if($tp[4] eq $info[4]){
			#print "a";
			if($tp[5]<$info[6]){
				$cbca=0;
				$cb[$cbc][2]=$tp[6]>$cb[$cbc][2]?$tp[6]:$cb[$cbc][2];		
				if($tp[0] eq $info[0]){
					#combin
					#print ".";
					$C++;		
					$add=0;
					$tp[5]=$info[5];
					$tp[6]=$tp[6]>$info[6]?$tp[6]:$info[6];
					$tp[1]=$tp[6]-$tp[5]+1;
					print OUTT ".\t$lineIN\n";
				}else{
					#cut info
					if($tp[6]<$info[6]){#inclued	
						#print "i";	
						print OUTT "i\t$lineIN\n";
						$C++;	
						#discard tp
						$add=0;	
						my $len=$tp[6]-$tp[5]+1;
						if($o{$tp[$sb]}>0){
							$o{$tp[$sb]}+=$len;
						}	else{
							$o{$tp[$sb]}=$len;
						}
						@tp=@info;						
						
					}else{	
						if($tp[5]==$info[5]){
							print OUTT "i\t$lineIN\n";
							#print "i";
							$C++;
							#discard info
							$add=0;	
							#@tp=@info;
							my $len=$tp[6]-$tp[5]+1;
							if($o{$tp[$sb]}>0){
								$o{$tp[$sb]}+=$len;
							}	else{
								$o{$tp[$sb]}=$len;
							}
						}else{
							#print "d";
							print OUTT "d\t$lineIN\n";
							my $d=$info[6]-$tp[5]-1;		
							$info[6]=$tp[5]-1;	
							$info[1]=$info[6]-$info[5]+1;									
							if($o{$info[$sb]}>0){
								$o{$info[$sb]}+=$d;
							}	else{
								$o{$info[$sb]}=$d;
							}
						}	
					}					
				}			
			}
		}
		if($add==1 && $info[1]>0){
			print OUT "$info[0]\t$info[1]\t$info[4]\t$info[5]\t$info[6]\t$info[8]\n";	
			$count+=$info[1];
		}
		if($cbca==1){
			$cbc++;
			$cb[$cbc][0]=$tp[4];
			$cb[$cbc][1]=$tp[5];
			$cb[$cbc][2]=$tp[6];
		}
		@info=@tp;				
    }
    print OUT "$info[0]\t$info[1]\t$info[4]\t$info[5]\t$info[6]\t$info[8]\n";	
    $count+=$info[1];
	my $ot=$ARGV[0].".delete";
	open(OUTD,'>',"$ot");
    my @t=keys(%o);
    my $t=@t;
    for(my $i=0;$i<$t;$i++){
    	print OUTD "$t[$i]\t$o{$t[$i]}\n";
    	$countD+=$o{$t[$i]};
    	
	}
	print OUTD "---------------\n";
	my $cbt=0;
	for(my $i=0;$i<=$cbc;$i++){
		my $l=$cb[$i][2]-$cb[$i][1]+1;
		$cbt+=$l;
		print OUTD "$cb[$i][0]\t$cb[$i][1]\t$cb[$i][2]\n";
	}
	print OUTD "---------------\nTotal length\t$count\nDeleat length\t$countD\nTotal Length Count\t$cbt\n";
	print "\n---------------\nTotal length\t$count\nDeleat length\t$countD\n$countA\nTotal Length Count\t$cbt\n";
    close(OUTD);
    close(OUTT);
}

sub RepeatMaskerCatOUTfamilyCombin{
	# sort by chrom and loc, add type
	#same family in region , combin
	#LTR	z_other	RM009_049	chr01	817122	818195	1074	2075	3172	3172	33.858764186633	-
	#LTR	z_other	RM009_049	chr01	818243	820273	2031	1	2074	3172	64.0290037831021	-
	#old 0-1- 2-family 3-type 4-chrom 5-st 6-ed 7-len 8-Tst 9-Ted 10-Tlen 11-percent 12-direct 
	#0-class 1-family 2-type 3-chrom 4-st 5-ed 6-len 7-Tst 8-Ted 9-Tlen 10-percent 11-direct 
	#out put same list
	my $range=0.4;   
	if($_[1]>0){
		$range=$_[1];
	}
	my $overlap=0.2;	
	if($_[2]>0){
		$overlap=$_[2];
	}
	my $count=0;
	my @info;
	
	while($lineIN=<IN>){
		chomp($lineIN);
     		my @tp= &Getheadinfo($lineIN);
     		#my @name=$tp[2]."#".$tp[0]."/".$tp[1];
     		#$tp[0]=$name[0];
     		my $l=@info;
     		my $add=1;
     		#another chromsome	
		if($tp[0] eq $info[0] && $tp[1] eq $info[1] && $tp[2] eq $info[2] && $tp[3] eq $info[3] && $tp[11] eq $info[11]){
			print "a";
			if( $tp[11] eq "+"){
				my $d=$tp[7]-$info[8];
				my $o=-1*$overlap*$tp[9];
				my $r=$range*$tp[9];
				#my $r=$range;
				if($d>$o && $d<$r){
					if($r>($tp[4]-$info[5])){
						#combin
						print ".";
						$add=0;
						$tp[4]=$info[4];
						$tp[7]=$info[7];
						$tp[6]=$tp[5]-$info[4]+1;
						$tp[10]=$tp[6]/$info[9];
					}
				}
			}else{
				my $d=$info[7]-$tp[8];
				my $o=-1*$overlap*$tp[9];
				my $r=$range*$tp[9];
				#my $r=$range;
				if($d>$o && $d<$r){
					if($r>($tp[4]-$info[5])){
						#combin
						print ".";
						$add=0;
						$tp[4]=$info[4];
						$tp[8]=$info[8];
						$tp[6]=$tp[5]-$info[4]+1;
						$tp[10]=$tp[6]/$info[9];
					}
				}				
			}		
		}
		if($add==1 && $count>0){
			print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t$info[9]\t$info[10]\t$info[11]\n";			
		}
		@info=@tp;
		$count++;			
     	}
     	print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t$info[9]\t$info[10]\t$info[11]\n";
     	print "\n";
}

sub RepeatMaskerCatLTRCombin{
	# sort by chrom and loc
	#same family in region , combin
	#RM004	LTR	chr01	42171	42588	417	1	418	418	99.7607655502392	-
	#RM004	int	chr01	42589	53014	10425	1	10450	10450	99.7607655502392	+
	#Singl	int	chr01	42589	42588	208	9226	9500	18279	1.13791782920291	-
	#RM004	LTR	chr01	53015	53453	438	1	418	418	104.784688995215	-
	#0-family 1-type 2-chrom 3-st 4-ed 5-len 6-Tst 7-Ted 8-Tlen 9-percent 10-direct
	my $LTRrange=20000;   #two LTR max distance
	if($_[1]>0){
		$LTRrange=$_[1];
	}
	print "LTR combin range $LTRrange\n";
    my $range=200;
	if($_[2]>0){
		$range=$_[2];
	}
	print "family combin range $range\n";	
	my $st=0;
	if($_[3] eq "S" || $_[3] eq "s"){
		#strict  combin only L+I+L
		$st=1;
		print "Use strect\n";
	}	     	
			
	
	my @info;  #two-dimensional array
	my $pcount=0;
	my $l=0;
	#combin full LTR 
	#output
	#0-family 1-typr 2-chrom 3-st1 4-ed1 5-st2 6-ed2 7-len 7-direct
	while($lineIN=<IN>){
		chomp($lineIN);
     		my @tp= &Getheadinfo($lineIN);
     		$l=@info;
     		#print "$l ";
     		my $add=1;
			#RM001	LTR	Chr01	30315942	30315956	15	1	16	1195	1.25523012552301	+
     		#0-family 1-type 2-chrom 3-st 4-ed 5-len 6-Tst 7-Ted 8-Tlen 9-percent 10-direct
     		for(my $i=$l-1; $i>0;$i--){
     			#print $info[$i][0]." ";
     			if($tp[2] eq $info[$i][2] && ($info[$i][4]+$LTRrange) > $tp[3]){
     				if($tp[1] eq $info[$i][1] && $tp[0] eq $info[$i][0] && $tp[10] eq $info[$i][10] && $tp[1] eq "LTR"){
     					#combin
     					my $sb="L+L";
     					for(my $j=$i;$j<$l;$j++){
     						if($tp[0] eq $info[$j][0] && $info[$j][1] eq "int"){
     							if($st ==0){
     								if($tp[10] eq $info[$j][10]){
     									$sb="L+I+L";
     								}else{
     									$sb="L-I-L";
     								}
     								last;
     							}else{
      								if($tp[10] eq $info[$j][10]){
     									$sb="L+I+L";
     									last;
     								}    								
     							}
     						}
     					}
     					if($sb=="L+I+L"){
     						my $len=$tp[4]-$info[$i][3]+1;
						print OUT "$tp[0]\t$sb\t$tp[2]\t$info[$i][3]\t$info[$i][4]\t$tp[3]\t$tp[4]\t$len\t$tp[8]\t$tp[10]\n";
						splice(@info,$i,$l-$i);	
						print ".";
						$pcount++;
						$add=0;
						last;
					}
     				}
     			}else{
     				last;				
     			}
     		}
     		if($add==1){
     			$l=@info;
     			$info[$l]=[@tp];
     			#print $info[$l][0]." ";
     		}
     	}
     	$l=@info;
     	my $Acount=0;
     	my $Ocount=0;
     	#combin fragment same family
     	if($info[0][1] eq "LTR"){
     		$info[0][1]="L";	
     	}else{
     		$info[0][1]="I";
     	}	
     	for(my $i=1;$i<$l;$i++){
     		my $add=1;
     		if($info[$i][1] eq "LTR"){
     			$info[$i][1]="L";	
     		}else{
     			$info[$i][1]="I";
     		}
     		if($info[$i][0] eq $info[$i-1][0] && $info[$i][2] eq $info[$i-1][2]){
     			if(($info[$i-1][4]+$range)<$info[$i][3]){
     				$add=0;
     				$Acount++;
     				$info[$i][1]=$info[$i][1].$info[$i-1][1];
     				$info[$i][3]=$info[$i-1][3];
				$info[$i][5]=$info[$i][4]-$info[$i][3];
				print "c";
     			}
     		}
     		if($add==1){
     			#0-family 1-type 3-st 4-ed 5-len 6-direct
     			$Ocount++;
     			print OUT "$info[$i-1][0]\t$info[$i-1][1]\t$info[$i-1][2]\t$info[$i-1][3]\t$info[$i-1][4]\t$info[$i-1][5]\t$info[$i-1][10]\n";	
     		}	
     	} 
     	my $f= $Ocount-$Acount;
     	print "\nL-I-L	$pcount\n";
     	print "Combin	$Acount\n";
     	print "Partai	$Ocount\n";
}



sub RepeatMaskerCatCombinCount{
	#input cat file
	#must sorted by 3-scafold 4-from
	#LTR remove family member name
	#e.g. RM001-001 to RM001 
	#825 19.71 13.94 0.00 scaff_2259 277519 277726 (5553466) rnd-5_family-7520#RC/Helitron 1 237 (80) 5
	#1062 12.57 0.00 1.51 scaff_2259 277738 277939 (5553253) C rnd-2_family-228#DNA (5) 201 3 5
	#output 
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct
	#lines which begin with ## -- drop
	#different type save to different file
	#overlap cut to high score
	@rm; #two-dimensional array
	$count=0; 
	@type;
	$typecount=0;
	@lastline;
	$lastline[0]= "NA";
	my $n=0;
	#read in cat file
	while($lineIN=<IN>){
		$n++;
		chomp($lineIN);
		if(substr($lineIN, 0, 1) ne "#"){
     			my @tp= &Getheadinfo($lineIN);
			my $tp=@tp;
			my @name;
			my $direct;
			my $drop=0;
			#add info    
			my $Tlength;
			my $Tfrom;
			my $Tto;
			if ($tp[8] ne "C"){
				$direct="+";
				@name=&Getheadinfo($tp[8],"#","/");
				my $name=@name;
				if($name<3){
					for(my $j=$name;$j<3;$j++){
						$name[$j]="z_other";
					}
				}
				if($name[1] eq "LTR"){
					my @tpname=&Getheadinfo($name[0],"_", "-");
					my $tpname=@tpname;
					if($tpname==3){
						$name[0]=$tpname[0]."-".$tpname[2];
					}
				}
				my $remain=substr(substr($tp[11],1,(length($tp[11])-1)),0,(length($tp[11])-2));
				$Tlength=$tp[10]+$remain;
				$Tfrom=$tp[9];
				$Tto=$tp[10];				
			}else{
				#complement
				$direct="-";
				@name=&Getheadinfo($tp[9],"#","/");
				my $name=@name;
				if($name<3){
					for(my $j=$name;$j<3;$j++){
						$name[$j]="z_other";
					}
				}				
				#remove LTR family member name
				if($name[1] eq "LTR"){
					my @tpname=&Getheadinfo($name[0],"_", "-");
					my $tpname=@tpname;
					if($tpname==3){
						$name[0]=$tpname[0]."-".$tpname[2];
					}
				}
				#print "@name\n";				
				my $remain=substr(substr($tp[10],1,(length($tp[10])-1)),0,(length($tp[10])-2));
				$Tlength=$remain+$tp[11];
				$Tfrom=$tp[12];
				$Tto=$tp[11];									
			}			
			#overlap cut
			#print "$name[0] $lastline[2]\n";
			#@lastline: 0-1-2-name 3-chrom 4-start 5-end 6-len 7-Tfrom 8-Tto 9-Tlen 10-percent 11-direct 12-score
			#@tp:       0-score 4-chrom 5-st 6-end
			#if new start before old end, cut new
			if($lastline[5]>=$tp[5] && $lastline[0] ne "NA" && $tp[4] eq $lastline[3]){
				if($lastline[5]>=$tp[6]){
					#drop new
					$drop=1;
					print "d";
				}elsif($lastline[4]>=$tp[5]){
					#drap last line
					print "D";
				}else{
					#same name, same direct combin
					if($name[0] eq $lastline[2] && $direct eq $lastline[11]){						
						$tp[5]=$lastline[4];
						if ($direct eq "+"){
							$Tfrom=$lastline[7];	
						}else{
							#complement
							$remain=$lastline[9];
							$Tto=$lastline[8];									
						}
						print "c";						
					}else{
						#LTR  combin to nest, leave nest info in name
						if($lastline[0] eq $name[1] && $lastline[0] eq "LTR"){
							print "n";
							$drop=1;
							#print "\ntest $lastline[0]\n";
							$lastline[2]="@".$lastline[2]."#$name[0]=$name[1]=$name[2]";
							$lastline[5]=$tp[6];
							$lastline[6]=$tp[6]-$lastline[4]+1;
						}else{  #else no combin
							#do nothing now
							$rm[$count]="$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";
							print OUT "$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";;
							$count++;
							print "+";							
						}
						#old
						#cut by length, no score
						#my $add=$lastline[5]-$tp[5];
						#print " $add ";
						##cut @tp
						#if($lastline[12] > $tp[0]){								
						#	$tp[5]=$tp[5]+$add+1;
						#	if ($direct eq "+"){
						#		$Tfrom=$Tfrom+$add+1;
						#	}else{
						#		$Tto=$Tto-$add-1;
						#	}
						#}#cut @lastline
						#else{
						#	$lastline[5]=$lastline[5]-$add-1;
						#	if ($direct eq "+"){
						#		$lastline[8]=$lastline[8]-$add-1;
						#	}else{
						#		$lastline[7]=$lastline[7]+$add+1;	
						#	}
						#	$lastline[6]=$lastline[6]-$add-1;
						#	$lastline[10]=$lastline[6]*100/$lastline[9];														
						#}
						#$rm[$count]="$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";
						#print OUT "$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";;
						#$count++;
						#print ".";							
					}					
				}
						
			}else{
					if($lastline[0] ne "NA"){
					$rm[$count]="$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";
					print OUT "$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";;
					$count++;
					print "+";
				}
			}
						
			my $len=$tp[6]-$tp[5]+1;
			my $percent=$len*100/$Tlength;
			
			#$rm[$count]="$name[1]\t$name[2]\t$name[0]\t$tp[4]\t$tp[5]\t$tp[6]\t$len\t$Tfrom\t$Tto\t$Tlength\t$percent\t$direct\n";
			
			if($drop==0){
				@lastline=($name[1],$name[2],$name[0],$tp[4],$tp[5],$tp[6],$len,$Tfrom,$Tto,$Tlength,$percent,$direct,$tp[0]);
				my $canadd=1;
				#print "\n@lastline\n";
			}
			for(my $i=0; $i<$typecount;$i++){
				if($type[$i] eq $name[1]){
					$canadd=0;
					last;
				}
			}
			if($canadd == 1){
				$type[$typecount]=$name[1];
				$typecount++;
			}
		}else{
			last;
		}
	}
	$rm[$count]="$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";
	print OUT "$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";;
	print OUT "-------------------------------------------\n";
	print "+";
	$count++;
	print "\n$count\n";
	#classify result
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct
	@type = sort(@type);
	@rm=sort(@rm);
	#percent count
	my @p75; #0-number 1-length
	my @p50;
	my @p25;
	my @p00; #total
	my $tptype="99999_start_symble";
	my $tpfamily;
	my @t75; #0-number 1-length
	my @t50;
	my @t25;
	my @t00; #total
	$p75[0]=0;
	$p50[0]=0;
	$p25[0]=0;
	$p00[0]=0;
	$p75[1]=0;
	$p50[1]=0;
	$p25[1]=0;
	$p00[1]=0;
	$t75[0]=0;
	$t50[0]=0;
	$t25[0]=0;
	$t00[0]=0;
	$t75[1]=0;
	$t50[1]=0;
	$t25[1]=0;
	$t00[1]=0;	
	
	my @stype;
	my $ot;
	for(my $i=0; $i<$count;$i++){
		my @tp= &Getheadinfo($rm[$i]);
		if($tptype ne $tp[0]){
			#new type
			if($tptype ne "99999_start_symble"){
				#output family count
				print OUT "\t$tpfamily\t$t00[0]($t00[1])\t$t25[0]($t25[1])\t$t50[0]($t50[1])\t$t75[0]($t75[1])\n";
				#output type count
				print OUT "\tCount\t$p00[0]($p00[1])\t$p25[0]($p25[1])\t$p50[0]($p50[1])\t$p75[0]($p75[1])\n";
				#clean old type
				close(OUTD);
				$p75[0]=0;
				$p50[0]=0;
				$p25[0]=0;
				$p00[0]=0;
				$p75[1]=0;
				$p50[1]=0;
				$p25[1]=0;
				$p00[1]=0;
				$t75[0]=0;
				$t50[0]=0;
				$t25[0]=0;
				$t00[0]=0;
				$t75[1]=0;
				$t50[1]=0;
				$t25[1]=0;
				$t00[1]=0;
			}
			
			#open new type
			$tptype=$tp[0];
			$tpfamily=$tp[1];		
			$ot=$ARGV[0].".".$tptype.".detail";
			open(OUTD,'>',"$ot");
			print OUT "$tptype\n";
					
		}else{			
			if($tpfamily ne $tp[1]){
				#output old family
				print OUT "\t$tpfamily\t$t00[0]($t00[1])\t$t25[0]($t25[1])\t$t50[0]($t50[1])\t$t75[0]($t75[1])\n";
				#new family
				$tpfamily=$tp[1];
				$t75[0]=0;
				$t50[0]=0;
				$t25[0]=0;
				$t00[0]=0;
				$t75[1]=0;
				$t50[1]=0;
				$t25[1]=0;
				$t00[1]=0;		
			}
		}
		#add info
		#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct
		$stype[$i]="$tp[3]\t$tp[4]\t$tp[5]\t$tp[0]#$tp[1]/$tp[2]";
		print OUTD $rm[$i]."\n";	
		if($tp[10] > 75){
			$p75[0]++;
			$p75[1]=$p75[1]+$tp[6];
			$t75[0]++;
			$t75[1]=$t75[1]+$tp[6];
		}
		if($tp[10] > 50){
			$p50[0]++;
			$p50[1]=$p50[1]+$tp[6];
			$t50[0]++;
			$t50[1]=$t50[1]+$tp[6];
		}
		if($tp[10] > 25){
			$p25[0]++;
			$p25[1]=$p25[1]+$tp[6];
			$t25[0]++;
			$t25[1]=$t25[1]+$tp[6];
		}			
		$p00[0]++;
		$p00[1]=$p00[1]+$tp[6];	
		$t00[0]++;
		$t00[1]=$t00[1]+$tp[6];		
	}
	print OUT "\t$tpfamily\t$t00[0]($t00[1])\t$t25[0]($t25[1])\t$t50[0]($t50[1])\t$t75[0]($t75[1])\n";
	#output type count
	print OUT "\tCount\t$p00[0]($p00[1])\t$p25[0]($p25[1])\t$p50[0]($p50[1])\t$p75[0]($p75[1])\n";
	#clean old type
	close(OUTD);
	
	my $otall=$ARGV[0].".all.detail";
	open(OUTA,'>',"$otall");
	@stype = sort(@stype);
	#0-scaf 1-start 2-end 3-type#family/name
	for(my $i=0; $i<$count;$i++){
		print OUTA "$stype[$i]\n";	
	}
	close(OUTA);
	print "\n";
}


sub RepeatMaskerCatnoCombinCount{
	#input cat file
	#must sorted by 3-scafold 4-from
	#LTR remove family member name
	#e.g. RM001-001 to RM001 
	#825 19.71 13.94 0.00 scaff_2259 277519 277726 (5553466) rnd-5_family-7520#RC/Helitron 1 237 (80) 5
	#1062 12.57 0.00 1.51 scaff_2259 277738 277939 (5553253) C rnd-2_family-228#DNA (5) 201 3 5
	#output 
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct 12-div
	#lines which begin with ## -- drop
	#different type save to different file
	#overlap cut to high score
	@rm; #two-dimensional array
	$count=0; 
	@type;
	$typecount=0;
	my $n=0;
	my %family;
	my @name;
	my $cut="_";
	if ($_[1]){
		$cut=$_[1];
	}
	
	my $direct; 
	my $Tlength;
	my $Tfrom;
	my $Tto;
	my @tpf;
	my @tp;
	my $tp;	
	
	#read in cat file
	while($lineIN=<IN>){
		$n++;
		chomp($lineIN);
		@name=();
		if(substr($lineIN, 0, 1) ne "#"){
     		@tp= &Getheadinfo($lineIN);
			$tp=@tp;

			if ($tp[8] ne "C"){
				$direct="+";
				@name=&Getheadinfo($tp[8],"#","/");
				my $name=@name;
				if($name<3){
					for(my $j=$name;$j<3;$j++){
						$name[$j]="z_other";
					}
				}
				
				my $remain=substr(substr($tp[11],1,(length($tp[11])-1)),0,(length($tp[11])-2));
				$Tlength=$tp[10]+$remain;
				$Tfrom=$tp[9];
				$Tto=$tp[10];				
			}else{
				#complement
				$direct="-";
				@name=&Getheadinfo($tp[9],"#","/");
				my $name=@name;
				if($name<3){
					for(my $j=$name;$j<3;$j++){
						$name[$j]="z_other";
					}
				}				

				my $remain=substr(substr($tp[10],1,(length($tp[10])-1)),0,(length($tp[10])-2));
				$Tlength=$remain+$tp[11];
				$Tfrom=$tp[12];
				$Tto=$tp[11];									
			}	
			
									
			my $len=$tp[6]-$tp[5]+1;
			my $percent=$len*100/$Tlength;
			
			@tpf=&Getheadinfo($name[0],$cut);		
			if($family{$tpf[0]}>0){
			$family{$tpf[0]}+=$percent;
			}else{
				$family{$tpf[0]}=$percent;
				print "+";
			}
			
			
			#print "$tpf[0]\t$family[$tpf[0]]\n";
			$rm[$count]="$name[1]\t$name[2]\t$name[0]\t$tp[4]\t$tp[5]\t$tp[6]\t$len\t$Tfrom\t$Tto\t$Tlength\t$percent\t$direct\t$tp[0]\t$tp[1]\n";
			print OUT "$name[1]\t$name[2]\t$name[0]\t$tp[4]\t$tp[5]\t$tp[6]\t$len\t$Tfrom\t$Tto\t$Tlength\t$percent\t$direct\t$tp[0]\t$tp[1]\n";
			$count++;
			print "+";
		}else{
			last;
		}
	}

	print "\n$count\n";
	
	print OUT "\n\nTE type Count------------------\n";
	
	#classify result
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct
	@type = sort(@type);
	@rm=sort(@rm);
	#percent count
	my @p00; #total
	my $tptype="99999_start_symble";
	my $tpfamily;
	my @t00; #0-number 1-length
	$p00[0]=0;
	$p00[1]=0;
	$t00[0]=0;
	$t00[1]=0;	
	
	my @stype;
	my $ot;
	for(my $i=0; $i<$count;$i++){
		my @tp= &Getheadinfo($rm[$i]);
		if($tptype ne $tp[0]){
			#new type
			if($tptype ne "99999_start_symble"){
				#output family count
				print OUT "\t$tpfamily\t$t00[0]($t00[1])\n";
				#output type count
				print OUT "\tCount\t$p00[0]($p00[1])\n";
				#clean old type
				close(OUTD);
				$p00[0]=0;
				$p00[1]=0;
				$t00[0]=0;
				$t00[1]=0;
			}
			
			#open new type
			$tptype=$tp[0];
			$tpfamily=$tp[1];		
			$ot=$ARGV[0].".".$tptype.".detail";
			open(OUTD,'>',"$ot");
			print OUT "$tptype\n";
					
		}else{			
			if($tpfamily ne $tp[1]){
				#output old family
				print OUT "\t$tpfamily\t$t00[0]($t00[1])\n";
				#new family
				$tpfamily=$tp[1];
				$t00[0]=0;
				$t00[1]=0;		
			}
		}
		#add info
		#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct
		$stype[$i]="$tp[3]\t$tp[4]\t$tp[5]\t$tp[0]#$tp[1]/$tp[2]";
		print OUTD $rm[$i]."\n";	
		$p00[0]++;
		$p00[1]=$p00[1]+$tp[6];	
		$t00[0]++;
		$t00[1]=$t00[1]+$tp[6];		
	}
	print OUT "\t$tpfamily\t$t00[0]($t00[1])\n";
	#output type count
	print OUT "\tCount\t$p00[0]($p00[1])\n";
	#clean old type
	close(OUTD);
	
	my $otall=$ARGV[0].".all.detail";
	open(OUTA,'>',"$otall");
	@stype = sort(@stype);
	#0-scaf 1-start 2-end 3-type#family/name
	for(my $i=0; $i<$count;$i++){
		print OUTA "$stype[$i]\n";	
	}
	close(OUTA);
	
	#print family
	$of=$ARGV[0].".family";
	open(OUTF,'>',"$of");
	my @o=keys(%family);
	@o=sort(@o);
	my $o=@o;
	print "$o families counted\n";
	for(my $i=0;$i<$o;$i++){
		print OUTF "$o[$i]\t$family{$o[$i]}\n";
	}	
	close(OUTF);
	print "\n";
}

	
sub RepeatMaskerCatCount{
	#input cat file
	#must sorted by 3-scafold 4-from
	#LTR remove family member name
	#e.g. RM001-001 to RM001 
	#825 19.71 13.94 0.00 scaff_2259 277519 277726 (5553466) rnd-5_family-7520#RC/Helitron 1 237 (80) 5
	#1062 12.57 0.00 1.51 scaff_2259 277738 277939 (5553253) C rnd-2_family-228#DNA (5) 201 3 5
	#output 
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct 12-div
	#lines which begin with ## -- drop
	#different type save to different file
	#overlap cut to high score
	@rm; #two-dimensional array
	$count=0; 
	@type;
	$typecount=0;
	@lastline;
	$lastline[0]= "NA";
	my $n=0;
	#read in cat file
	while($lineIN=<IN>){
		$n++;
		chomp($lineIN);
		if(substr($lineIN, 0, 1) ne "#"){
     			my @tp= &Getheadinfo($lineIN);
			my $tp=@tp;
			my @name;
			my $direct;
			my $drop=0;
			#add info    
			my $Tlength;
			my $Tfrom;
			my $Tto;
			if ($tp[8] ne "C"){
				$direct="+";
				@name=&Getheadinfo($tp[8],"#","/");
				my $name=@name;
				if($name<3){
					for(my $j=$name;$j<3;$j++){
						$name[$j]="z_other";
					}
				}

				my $remain=substr(substr($tp[11],1,(length($tp[11])-1)),0,(length($tp[11])-2));
				$Tlength=$tp[10]+$remain;
				$Tfrom=$tp[9];
				$Tto=$tp[10];				
			}else{
				#complement
				$direct="-";
				@name=&Getheadinfo($tp[9],"#","/");
				my $name=@name;
				if($name<3){
					for(my $j=$name;$j<3;$j++){
						$name[$j]="z_other";
					}
				}				
			
				my $remain=substr(substr($tp[10],1,(length($tp[10])-1)),0,(length($tp[10])-2));
				$Tlength=$remain+$tp[11];
				$Tfrom=$tp[12];
				$Tto=$tp[11];									
			}			
			#overlap cut
			#print "$name[0] $lastline[2]\n";
			#@lastline: 0-1-2-name 3-chrom 4-start 5-end 6-len 7-Tfrom 8-Tto 9-Tlen 10-percent 11-direct 12-score
			#@tp:       0-score 4-chrom 5-st 6-end
			#if new start before old end, cut new
			if($lastline[5]>=$tp[5] && $lastline[0] ne "NA" && $tp[4] eq $lastline[3]){
				if($lastline[5]>=$tp[6]){
					#drop new
					$drop=1;
					print "d";
				}elsif($lastline[4]>=$tp[5]){
					#drap last line
					print "D";
				}else{
					#same name, same direct combin
					if($name[0] eq $lastline[2] && $direct eq $lastline[11]){						
						$tp[5]=$lastline[4];
						if ($direct eq "+"){
							$Tfrom=$lastline[7];	
						}else{
							#complement
							$remain=$lastline[9];
							$Tto=$lastline[8];									
						}
						print "c";						
					}else{

							#do nothing now
							$rm[$count]="$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";
							print OUT "$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";;
							$count++;
							print "+";							
						
					}					
				}
						
			}else{
					if($lastline[0] ne "NA"){
					$rm[$count]="$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";
					print OUT "$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";;
					$count++;
					print "+";
				}
			}
						
			my $len=$tp[6]-$tp[5]+1;
			my $percent=$len*100/$Tlength;
			
			#$rm[$count]="$name[1]\t$name[2]\t$name[0]\t$tp[4]\t$tp[5]\t$tp[6]\t$len\t$Tfrom\t$Tto\t$Tlength\t$percent\t$direct\n";
			
			if($drop==0){
				@lastline=($name[1],$name[2],$name[0],$tp[4],$tp[5],$tp[6],$len,$Tfrom,$Tto,$Tlength,$percent,$direct,$tp[0]);
				my $canadd=1;
				#print "\n@lastline\n";
			}
			for(my $i=0; $i<$typecount;$i++){
				if($type[$i] eq $name[1]){
					$canadd=0;
					last;
				}
			}
			if($canadd == 1){
				$type[$typecount]=$name[1];
				$typecount++;
			}
		}else{
			last;
		}
	}
	$rm[$count]="$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";
	print OUT "$lastline[0]\t$lastline[1]\t$lastline[2]\t$lastline[3]\t$lastline[4]\t$lastline[5]\t$lastline[6]\t$lastline[7]\t$lastline[8]\t$lastline[9]\t$lastline[10]\t$lastline[11]\n";;
	print OUT "-------------------------------------------\n";
	print "+";
	$count++;
	print "\n$count\n";
	#classify result
	#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct
	@type = sort(@type);
	@rm=sort(@rm);
	#percent count
	my @p75; #0-number 1-length
	my @p50;
	my @p25;
	my @p00; #total
	my $tptype="99999_start_symble";
	my $tpfamily;
	my @t75; #0-number 1-length
	my @t50;
	my @t25;
	my @t00; #total
	$p75[0]=0;
	$p50[0]=0;
	$p25[0]=0;
	$p00[0]=0;
	$p75[1]=0;
	$p50[1]=0;
	$p25[1]=0;
	$p00[1]=0;
	$t75[0]=0;
	$t50[0]=0;
	$t25[0]=0;
	$t00[0]=0;
	$t75[1]=0;
	$t50[1]=0;
	$t25[1]=0;
	$t00[1]=0;	
	
	my @stype;
	my $ot;
	for(my $i=0; $i<$count;$i++){
		my @tp= &Getheadinfo($rm[$i]);
		if($tptype ne $tp[0]){
			#new type
			if($tptype ne "99999_start_symble"){
				#output family count
				print OUT "\t$tpfamily\t$t00[0]($t00[1])\t$t25[0]($t25[1])\t$t50[0]($t50[1])\t$t75[0]($t75[1])\n";
				#output type count
				print OUT "\tCount\t$p00[0]($p00[1])\t$p25[0]($p25[1])\t$p50[0]($p50[1])\t$p75[0]($p75[1])\n";
				#clean old type
				close(OUTD);
				$p75[0]=0;
				$p50[0]=0;
				$p25[0]=0;
				$p00[0]=0;
				$p75[1]=0;
				$p50[1]=0;
				$p25[1]=0;
				$p00[1]=0;
				$t75[0]=0;
				$t50[0]=0;
				$t25[0]=0;
				$t00[0]=0;
				$t75[1]=0;
				$t50[1]=0;
				$t25[1]=0;
				$t00[1]=0;
			}
			
			#open new type
			$tptype=$tp[0];
			$tpfamily=$tp[1];		
			$ot=$ARGV[0].".".$tptype.".detail";
			open(OUTD,'>',"$ot");
			print OUT "$tptype\n";
					
		}else{			
			if($tpfamily ne $tp[1]){
				#output old family
				print OUT "\t$tpfamily\t$t00[0]($t00[1])\t$t25[0]($t25[1])\t$t50[0]($t50[1])\t$t75[0]($t75[1])\n";
				#new family
				$tpfamily=$tp[1];
				$t75[0]=0;
				$t50[0]=0;
				$t25[0]=0;
				$t00[0]=0;
				$t75[1]=0;
				$t50[1]=0;
				$t25[1]=0;
				$t00[1]=0;		
			}
		}
		#add info
		#0-type 1-family 2-name 3-scafold 4-from 5-to 6-length 7-Tfrom 8-Tto 9-Tlength 10-percent 11-direct
		$stype[$i]="$tp[3]\t$tp[4]\t$tp[5]\t$tp[0]#$tp[1]/$tp[2]";
		print OUTD $rm[$i]."\n";	
		if($tp[10] > 75){
			$p75[0]++;
			$p75[1]=$p75[1]+$tp[6];
			$t75[0]++;
			$t75[1]=$t75[1]+$tp[6];
		}
		if($tp[10] > 50){
			$p50[0]++;
			$p50[1]=$p50[1]+$tp[6];
			$t50[0]++;
			$t50[1]=$t50[1]+$tp[6];
		}
		if($tp[10] > 25){
			$p25[0]++;
			$p25[1]=$p25[1]+$tp[6];
			$t25[0]++;
			$t25[1]=$t25[1]+$tp[6];
		}			
		$p00[0]++;
		$p00[1]=$p00[1]+$tp[6];	
		$t00[0]++;
		$t00[1]=$t00[1]+$tp[6];		
	}
	print OUT "\t$tpfamily\t$t00[0]($t00[1])\t$t25[0]($t25[1])\t$t50[0]($t50[1])\t$t75[0]($t75[1])\n";
	#output type count
	print OUT "\tCount\t$p00[0]($p00[1])\t$p25[0]($p25[1])\t$p50[0]($p50[1])\t$p75[0]($p75[1])\n";
	#clean old type
	close(OUTD);
	
	my $otall=$ARGV[0].".all.detail";
	open(OUTA,'>',"$otall");
	@stype = sort(@stype);
	#0-scaf 1-start 2-end 3-type#family/name
	for(my $i=0; $i<$count;$i++){
		print OUTA "$stype[$i]\n";	
	}
	close(OUTA);
	print "\n";
}


sub Selectfastqbylength{
	# lt 80
	my $cutlen=100;
	my $operate="gt";
	if(length($_[1])>0){
		$operate=$_[1];
		$cutlen=$_[2];
	}
	print "select fastq $operate $cutlen bp\n";
	my $name;
	while($lineIN=<IN>){
		chomp($lineIN);	
		$name=$lineIN;
		$lineIN=<IN>;
		chomp($lineIN);
		my $len=length($lineIN);
		if($operate eq "gt"){
			if($len>=$cutlen){
				print OUT $name."\n".$lineIN."\n";
				$lineIN=<IN>;
				print OUT $lineIN;
				$lineIN=<IN>;
				print OUT $lineIN;
			}else{
				$lineIN=<IN>;
				$lineIN=<IN>;
			}
		}else{
			if($len<=$cutlen){
				print OUT $name."\n".$lineIN."\n";
				$lineIN=<IN>;
				print OUT $lineIN;
				$lineIN=<IN>;
				print OUT $lineIN;
			}else{
				$lineIN=<IN>;
				$lineIN=<IN>;
			}		
		}		
    }
}


sub selectsequencebylength{
	my $name;
	my $seq;
	my $outd=$ARGV[0].".delete";
	open(OUTD,'>',"$outd");
	while($lineIN=<IN>){
		chomp($lineIN);	
     	if ($lineIN=~/>/){
     		my $len=length($seq);	
     		if ($len >= $_[1]){
     			print OUT $name."\n";
     			print OUT $seq."\n";
     		}else{
     			print OUTD $name."\n";
     			print OUTD $seq."\n";
     		}
     		$name=$lineIN;
     		$seq="";
     	}else{
     		$seq=$seq.$lineIN;
     	}
    }
    close(OUTD);
}


sub GetLineinEverageRegion{
	#1-class 2-value 3-region
	#0 2 3000
	my @r;
	my $name;
	my $min;
	my $avg=0;
	my $max;
	my $count=0;
	my $len=0;
	
	my $outd=$ARGV[0].".delete";
	open(OUTD,'>',"$outd");
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		if($len == 0){
			$len=@tp;
		}
		if($name ne $tp[$_[1]]){
				if($count>0){
				#compair region
				$avg= $avg/$count;
				$min=$avg-$_[3]/2;
				$max=$avg+$_[3]/2;
				print "$name\tAvg=$avg Min=$min Max=$max\n";
				#output
				for(my$i=0;$i<$count;$i++){
					if($min<$r[$i][$_[2]] && $max>$r[$i][$_[2]]){
						print OUT "$r[$i][0]";
						for(my $j=1;$j<$len;$j++){
							print OUT "\t$r[$i][$j]";
						}
						print OUT "\n";
					}else{
						print OUTD "$r[$i][0]";
						for(my $j=1;$j<$len;$j++){
							print OUTD "\t$r[$i][$j]";
						}
						print OUTD "\n";		
					}
				}
				#clean
				$count=0;
				@r=();
				$avg=0;
			}
			$name = $tp[$_[1]];
		}
		$avg+=$tp[$_[2]];
		$r[$count]=[@tp];
		$count++;
    }
	$avg= $avg/$count;
	$min=$avg-$_[3]/2;
	$max=$avg+$_[3]/2;
	print "$name\tAvg=$avg Min=$min Max=$max\n";
	#output
	for(my$i=0;$i<$count;$i++){
		if($min<$r[$i][$_[2]] && $max>$r[$i][$_[2]]){
			print OUT "$r[$i][0]";
			for(my $j=1;$j<$len;$j++){
				print OUT "\t$r[$i][$j]";
			}
			print OUT "\n";
		}else{
			print OUTD "$r[$i][0]";
			for(my $j=1;$j<$len;$j++){
				print OUTD "\t$r[$i][$j]";
			}
			print OUTD "\n";			
		}
	} 
	close(OUTD); 
}


sub Countline_value_No{

	while($lineIN=<IN>){
		chomp($lineIN);	
		#my @tp= &Getheadinfo($lineIN);
		my @tp=split(/\t/,$lineIN);
		my $tp=@tp;
		my $n=0;
		for(my $i=0;$i<$tp;$i++){
			if($_[1] eq "gt"){
				if($tp[$i] >= $_[2]){
					$n++;
				}
			}elsif($_[1] eq "lt"){
				if($tp[$i] <= $_[2]){
					$n++;
				}				
			}elsif($_[1] eq "eq"){
				if($tp[$i] == $_[2]){
					$n++;
				}				
			}
		}			
		print OUT "$tp[0]\t$n\n";

	}	
}


sub count_seq_Symble_No{
	my @seq;
	my @name;	
	my $count=-1;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){	
			$count++;
			$name[$count]=$lineIN;	
			$seq[$count]="";	
		}else{
			$seq[$count]=$seq[$count].$lineIN;
		}	
	}	
	
#	#check
#	$count=0;
#	foreach(@seq)
#	{
#		my %word=();
#	  while(m/([a-zA-Z])/g)
#	  {
#	    if(!exists $word{$1})
#	     {
#	        $word{$1}=1;
#	     }
#	     else
#	     {
#	         $word{$1}+=1;
#	     }
#	  }
#		#print
#		print OUT "$name[$count]";
#    foreach( sort keys %word)  #sort hash by keys letters ASC 
#    {           
#    	print OUT "\t$_=$word{$_}" ;  
#    }
#		print OUT "\n";
#		$count++;
#	}	


	#check
	$count=0;
	foreach(@seq)
	{
		print ".";
		my $az=0;
	  while(m/([a-z])/g)
	  {
			$az++;
	  }
	  my $len=length($_);
	  my $p=$az/$len;
		#print
		print OUT "$name[$count]\t$p\t$len\n";
		$count++;
	}	
	print "done\n";
}


sub CountlineSymbleNo{
	my @sc;	
	$sc[0]=@_;
	while($lineIN=<IN>){
		chomp($lineIN);
		for(my $j=1;$j<$sc[0];$j++){
			$sc[$j]=0;
		}
		#my @tp= &Getheadinfo($lineIN);
		my @tp=split(/\t/,$lineIN);
		my $tp=@tp;
		for(my $i=0;$i<$tp;$i++){
			for(my $j=1;$j<$sc[0];$j++){
				if(index($tp[$i],$_[$j],0)>=0){
					$sc[$j]++;
				}
			}
		}
		for(my $j=1;$j<$sc[0];$j++){
			print OUT "$sc[$j]\t"; 
		}
		print OUT "\n";
	}	
}



sub FastqToFasta{
	my @name;
	my @seq;
	my $len;
	my $count=-1;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my $name=substr($lineIN,1,length($lineIN));
		print OUT ">$name\n";
		$lineIN=<IN>;
		chomp($lineIN);	
		print OUT "$lineIN\n";
		$lineIN=<IN>;
		$lineIN=<IN>;
    }
}


sub FasToPhy{
	my @name;
	my @seq;
	my $len;
	my $maxnamelen=0;
	my $count=-1;
	while($lineIN=<IN>){
		chomp($lineIN);	
     	if ($lineIN=~/>/){
     		$count++;
     		my @tp= &Getheadinfo($lineIN,"|");	
 			$name[$count]=substr($tp[0],1,length($tp[0]));
			#$name[$count]=$tp[1];
 			my $tpnamelen=length($name[$count]);
 			if($tpnamelen>$maxnamelen){
 				$maxnamelen=$tpnamelen;
 			}
     	}else{
     		$seq[$count]=$seq[$count].$lineIN;
     	}
    }
    print "max name len = $maxnamelen\n";
    $count++;
    $len=length($seq[0]);
    print OUT " $count $len\n";
    for (my $i=0;$i<$count;$i++){
    	my $add=$maxnamelen-length($name[$i])+10;
    	#print "$add\n";
    	print OUT "$name[$i]";
    	for(my $j=0;$j<$add;$j++){
    		print OUT " ";
    	}
    	
    	print OUT "$seq[$i]\n";
    }
}

sub DndToTrees{
	# : start delete
	# , ) delete end
	while($lineIN=<IN>){
		chomp($lineIN);	
		@line=split('',$lineIN);
		my $char;
		my $canadd=1;
		foreach $char (@line){
			if($char eq ":"){
				$canadd=0;
			}
			if($char eq ")"){
				$canadd=1;	
				
			}
			if($char eq ","){
				$canadd=1;	
			}
			if($canadd==1){
				print OUT $char;
			}			
		}
		#print OUT "\n";
	}
}

sub StopCodenCount{
	my $name;
	my $seq;
	my $aa;
	my $symble="*";
	while($lineIN=<IN>){
		chomp($lineIN);	
     	if ($lineIN=~/>/){
     		my $len=length($seq);	
     		if ($len > 0){
				#translate
				$aa=&TranslateSeq($seq);
				#find stop coden
				my @r=&CountSymbleLocNumber($symble,$aa);
				#print "@r\n";
				print OUT "$name\t@r\n";
     		}
     		$name=$lineIN;
     		$seq="";
     		$aa="";
     	}else{
     		$seq=$seq.$lineIN;
     	}
    }
	#translate
	$aa=&TranslateSeq($seq);
	#find stop coden    
	my @r=&CountSymbleLocNumber($symble,$aa);
	print OUT "$name\t@r\n";
}

sub SoapDepthFastCount{
	#SRR077425.3903628       1       b       76      -       Chr1    1181    2       C->7G-27        C->11A-27
	#@dba save loc depth
	my $wlen=100;
	if($_[1]>0){
		$wlen=$_[1];
	}
	my @add;
	my $add=0;
	my $chrom='';
	my $loc=1;
	my $tploc=0;
	my $max=0;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($tp[5] ne $chrom){
			#output
			print "$max ";
			for(;$loc<=$max;$loc++){
				for(my $i=0;$i<$add;$i++){
					if($add[$i]<=$loc){
						splice(@add,$i,1);
					}
				}
				$add=@add;
				$count+=$add;
				if($loc%$wlen==0){
					my $o=$count/$wlen;
					my $oo=sprintf "%4.0f", $o;
					print OUT "$chrom\t$loc\t$oo\n";
					$count=0;
				}				
			}
			#output by window length
			if ($chrom ne ""){
				my $o=$count/$wlen;
				my $oo=sprintf "%4.0f", $o;
				print OUT "$chrom\t$loc\t$oo\n";
				$count=0;
			}
			@add=();
			$loc=1;
			#new chrom
			$chrom=$tp[5];
			$max=0;

		}
		#add
		$add=@add;
		for(my $i=$loc;$i<=$tp[6];$i++){
			#print "$add ";
			for(my $i=0;$i<$add;$i++){
				if($add[$i]<=$loc){
					splice(@add,$i,1);
				}
			}
			$add=@add;
			$count+=$add;
			if($loc%$wlen==0){
				my $o=$count/$wlen;
				my $oo=sprintf "%4.0f", $o;
				print OUT "$chrom\t$loc\t$oo\n";
				$count=0;
			}
			$loc++;
		}
		$add[$add]=$tp[6]+$tp[3];	
		if($max<$add[$add]){
			$max=$add[$add];
		}

		
	}
	print "$max\n";
	for(;$loc<=$max;$loc++){
		for(my $i=0;$i<$add;$i++){
			if($add[$i]<=$loc){
				splice(@add,$i,1);
			}
		}
		$add=@add;
		$count+=$add;
		if($loc%$wlen==0){
			my $o=$count/$wlen;
			my $oo=sprintf "%4.0f", $o;
			print OUT "$chrom\t$loc\t$oo\n";
		}
	}
	my $o=$count/$wlen;
	my $oo=sprintf "%4.0f", $o;
	print OUT "$chrom\t$loc\t$oo\n";

}

sub SoapDepthCount{
	#SRR077425.3903628       1       b       76      -       Chr1    1181    2       C->7G-27        C->11A-27
	#@dba save loc depth
	my $wlen=100;
	if($_[1]>0){
		$wlen=$_[1];
	}
	my @add;
	my $add=0;
	my $chrom='';
	my $loc=0;
	my $tploc=0;
	my $max=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= &Getheadinfo($lineIN);
		if($tp[5] ne $chrom){
			#output
			print "$max ";
			for(;$loc<$max;$loc++){
				$dba[$loc]=$add;
				for(my $i=0;$i<$add;$i++){
					if($add[$i]==$loc){
						splice(@add,$i,1);
					}
				}
				$add=@add;
			}
			#output by window length
			if ($chrom ne ""){
				print "\n$chrom\t$loc\n";
				&OutPutSoapDepth($wlen,$chrom,@dba);
			}
			@dba=(); 
			@add=();
			$loc=0;
			#new chrom
			$chrom=$tp[5];
			$max=0;

		}
		#add
		$add=@add;
		for(my $i=$loc;$i<$tp[6];$i++){
			$dba[$loc]=$add;
			#print "$add ";
			for(my $i=0;$i<$add;$i++){
				if($add[$i]==$loc){
					splice(@add,$i,1);
				}
			}
			$add=@add;
			$loc++;
		}
		$add[$add]=$tp[6]+$tp[3];	
		if($max<$add[$add]){
			$max=$add[$add];
		}
		
	}
	print "$max ";
	for(;$loc<$max;$loc++){
		$dba[$loc]=$add;
		for(my $i=0;$i<$add;$i++){
			if($add[$i]==$loc){
				splice(@add,$i,1);
			}
		}
		$add=@add;
	}
	#output by window length
	print "\n$chrom\t$loc\n";
	&OutPutSoapDepth($wlen,$chrom,@dba);	

}

sub OutPutSoapDepth{
	my $wlen=$_[0];
	my $count=1;
	my $len=@_;
	my $total=0;
	for($i=2;$i<$len;$i++){
		$total=$total+$_[$i];
		my $oloc=$wlen*$count;		
		if($i==$oloc+1){
			#output
			#print "$oloc\n";
			my $avg=$total/$wlen;
			
			print OUT "$_[1]\t$oloc\t$avg\n";
			$count++;
			$total=0;
		}	
	}
	my $avg=$total/$wlen;
	my $oloc=$wlen*$count;
	print OUT "$_[1]\t$oloc\t$avg\n";
}


sub SoapDirectCount{
	my $data=0;
	my $mate=0;
	my $same=0;
	my $nosame=0;
	while($lineIN=<IN>){
		my @tp1= &Getheadinfo($lineIN);
		$lineIN=<IN>;
		my @tp2= &Getheadinfo($lineIN);

		if($tp1[6] eq $tp2[6]){
			$same++;
		}else{
			$nosame++;
			if($tp1[6] eq "-"){
				if($tp1[8]<$tp2[8]){
					$mate++;	
				}else{
					$data++;
				}
			}else{
				if($tp1[8]>$tp2[8]){
					$mate++;	
				}else{
					$data++;
				}			
			}			
		}
		
	}
    print "-> <- data $data\n";
    print "<- -> data $mate\n";
    print "Same direct $same\n";
    print "Opposite direct $nosame\n";    
}

sub BlastHitLengthCount{
	#output querry and target length
	my $qname="";
	my $qlen=0;
	my @qloc;
	my @tloc;
	my $tlen=0;
	my $qcount=0;
	my $tcount=0;
	my $tname="start_9999999";
	my $info;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		#Jap_0257-LTR	Jap_1535	93.61	438	27	1	1	437	438	1	0.0	 638
		if ($tname ne $tp[1] || $qname ne $tp[0]){
			if($tname ne "start_9999999"){
				for(my $i=0;$i<$qcount;$i++){
					my $add=$qloc[$i][1]-$qloc[$i][0]+1;
					$qlen=$qlen+$add;
					#print "$qloc[$i][1] $qloc[$i][0] $add\n";
				}
				for(my $i=0;$i<$tcount;$i++){
					my $add=$tloc[$i][1]-$tloc[$i][0]+1;
					$tlen=$tlen+$add;
				}
				print OUT "$qname\t$qcount\t$qlen\t$tname\t$tcount\t$tlen\t$info\n";
				 
				
			}
			$qname=$tp[0];
			$tname=$tp[1];
			$count=0;
			$qlen=0;
			@qloc=();
			$tlen=0;
			@tloc=();
			$qcount=0;
			$tcount=0;
			$info="+"
		}
		#add loc		
		my $canadd=0;
		my @d;
		my $dc=0;
		#querry
		my $min=$tp[6];
		my $max=$tp[7];
		for(my $i=0;$i<$qcount;$i++){
			my @c=&Comparelocal($tname,$min,$max,$tname,$qloc[$i][0],$qloc[$i][1]);
			
			#print "$c[0] $min $max $qloc[$i][0] $qloc[$i][1]\n";
			if($c[0] != 6 && $c[0] != 7){
				#print ".";
				if($min>$qloc[$i][0]){
					$min=$qloc[$i][0];	
				}
				if($max<$qloc[$i][1]){
					$max=$qloc[$i][1];	
				}
				$qloc[$i][0]=0;
				$qloc[$i][1]=0;
			}
		}
		$qloc[$qcount][0]=$min;
		$qloc[$qcount][1]=$max;
		$qcount++;
		#print "$qloc[$count][0] $qloc[$count][1]\n";
		
		#target
		$min=$tp[8];
		$max=$tp[9];
		if($tp[8]>$tp[9]){
			$info="-";
			$min=$tp[9];
			$max=$tp[8];
		}
		for(my $i=0;$i<$tcount;$i++){
			my @c=&Comparelocal($tname,$min,$max,$tname,$tloc[$i][0],$tloc[$i][1]);			
			#print "$c[0] $min $max $tloc[$i][0] $tloc[$i][1]\n";
			if($c[0] != 6 && $c[0] != 7){
				#print ".";
				if($min>$tloc[$i][0]){
					$min=$tloc[$i][0];	
				}
				if($max<$tloc[$i][1]){
					$max=$tloc[$i][1];	
				}
				$tloc[$i][0]=0;
				$tloc[$i][1]=0;
			}
		}
		$tloc[$tcount][0]=$min;
		$tloc[$tcount][1]=$max;
		$tcount++;		
		 
    }
   	 #print "$tcount $qcount\n";
	for(my $i=0;$i<$qcount;$i++){
		my $add=$qloc[$i][1]-$qloc[$i][0]+1;
		$qlen+=$add;
	}
	for(my $i=0;$i<$tcount;$i++){
		my $add=$tloc[$i][1]-$tloc[$i][0]+1;
		$tlen=$tlen+$add;
	}
	print OUT "$qname\t$qcount\t$qlen\t$tname\t$tcount\t$tlen\t$info\n";
	print "\n";
}




sub BlastGetFirstHit{
	#@db  0-name1 1-name2
	
	my $count=0;
	my $lastline="NA";
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		#Jap_0257-LTR	Jap_1535	93.61	438	27	1	1	437	438	1	0.0	 638
		#remove self hit
		if($tp[0] ne $tp[1]){
			#get first hit
			my $tpline="$tp[0]\t$tp[1]";
			if($tpline ne $lastline){
				$lastline=$tpline;
				#check querry/target again
				my $canadd=1;
				for(my $i=0;$i<$count;$i++){
					if($db[$i][0] eq $tp[1] && $db[$i][1] eq $tp[0]){
						$canadd=0;
						last;
					}
				}				
				if($canadd==1){
					print OUT $lineIN."\n";
					$db[$count][0]=$tp[0];
					$db[$count][1]=$tp[1];
					$count++;
				}	
			}
  	}
  }
  print "$count first hit\n";
}

sub BlastGetDirect{
	#use first to change
	$name="";
	if($_[1]<0){
		$_[1]=1;
	}
	
	print OUT "Wrong direct----------------------\n";
	#@db; save ok direct
	my @mid; #save temp dir 0-q 1-t 2-direct
	my $add=1; 
	my $cm=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);
		#Jap_0257-LTR	Jap_1535	93.61	438	27	1	1	437	438	1	0.0	 638
		#8 > 9 verse
		my $tpd=$tp[8]<$tp[9]?1:-1;	
		my $tno=0;
		if($add==1){
			$add=2;
			$db{$tp[0]}='1';
			$db{$tp[1]}=$tpd;
			print "+";
		}else{
			if($db{$tp[0]}==1 ||$db{$tp[0]}==-1){
				#print "1";			
				if($db{$tp[1]}==1 ||$db{$tp[1]}==-1){					
					if($db{$tp[1]} ne ($db{$tp[0]}*$tpd)){
						print OUT "$lineIN\n";
						print "$lineIN\n";
					}
				}else{
					$db{$tp[1]}=$db{$tp[0]}*$tpd;
					print "+";
				}				
			}elsif($db{$tp[1]}==1 ||$db{$tp[1]}==-1){
				$db{$tp[0]}=$db{$tp[1]}*$tpd;
				print "+";
			}else{
				$mid[$cm][0]=$tp[0];
				$mid[$cm][1]=$tp[1];
				$mid[$cm][2]=$tp[8]<$tp[9]?1:-1;
				$cm++;
				#print "n";
				$tno=1;
			}
		}
		if($tno==0){
			for(my $i=0;$i<$cm;$i++){
				my $ok=0;
				if($db{$mid[$i][0]}==1 ||$db{$mid[$i][0]}==-1){
					$ok=1;
					if($db{$mid[$i][1]}==1 ||$db{$mid[$i][1]}==-1){
						if($db{$mid[$i][1]} ne ($db{$mid[$i][0]}*$mid[$i][2])){
							print OUT "$mid[$i][0]\t$mid[$i][1]\t$mid[$i][2]\n";
						}
					}else{
						$db{$mid[$i][1]}==$db{$mid[$i][0]}*$mid[$i][2];
					}
				}elsif($db{$mid[$i][1]}==1 ||$db{$mid[$i][1]}==-1){
					$ok=1;
					$db{$mid[$i][2]}==$db{$mid[$i][1]}*$mid[$i][2];
				}
				if($ok==1){
					$cm--;
					$mid[$i][0]=$mid[$cm][0];
					$mid[$i][1]=$mid[$cm][1];
					$mid[$i][2]=$mid[$cm][2];
					print ".";
				}
			}
		}

    }
    
    my @o=keys(%db);
    my $o=@o;
    print "$db directed\n$cm no linked\n";
    print OUT "$db directed----------------------\n";
    for(my $i=0;$i<$o;$i++){
    	print OUT "$o[$i]\t$db{$o[$i]}\n";
    }
    print OUT "$cm no linked----------------------\n";
    for(my $i=0;$i<$cm;$i++){
    	print OUT "$cm[$i][0]\t$cm[$i][1]\t$cm[$i][2]\n";
    }    
}



sub ReplaceFastAOtherLetter{
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($lineIN=~/>/){
			print OUT $lineIN."\n";
		}else{
			for (my $i = 0; $i < length($lineIN); $i++)
			{
				my $tp = substr($lineIN,$i,1);	
				if($tp eq 'A' || $tp eq 'a' || $tp eq 'T' || $tp eq 't' || $tp eq 'G' || $tp eq 'g' || $tp eq 'C' || $tp eq 'c'){
					print OUT $tp;
				}else{
					print OUT "-";
				}
    		}
    		print OUT "\n";
    	}
	}
}

sub GetGapInfo{
	my $char="-";
	if ($_[1]){
		$char=$_[1];
	}
	my $st;
	my $add=0;
	my $chrom;
	my $count=0;
	my $countlen=0;
	my $loc=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($lineIN=~/>/){
			my @tp= &Getheadinfo($lineIN);
			$chrom = substr($tp[0],1,length($tp[0]));
			$loc=0;
		}else{
			for (my $i = 0; $i < length($lineIN); $i++)
			{
				my $tp = substr($lineIN,$i,1);	
				if($tp eq $char){
					if($add==0){
						$st=$loc;
						$add=1;
					}
				}else{
    				if ($add == 1){
    					$count++;
    					my $len = $loc - $st +1;
    					$countlen += $len;
    					print OUT "$count\t$chrom\t$st\t$loc\t$len\n";
    					$add=0;
    				}					
				}
    			$loc++;
    		}
    	}
	}
	print OUT "Gap Number = $count\tTotal Gap Length = $countlen\n";
}


sub AddLineSameSymble{
	#(f for head) (b for behand)
	my $loc="f";
	if($_[2] eq "b"){
		$loc="b";
	}
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($loc eq "f"){
			print OUT "$_[1]\t$lineIN\n";
		}else{
			print OUT "$lineIN\t$_[1]\n";
		}
    }
}

sub AddLineSN{

	my $count=0;
	my $sb="";
	if($_[2]){
		$sb=$_[2];
	}
	while($lineIN=<IN>){
		chomp($lineIN);	
		$count++;
		print OUT "$sb$count\t$lineIN\n";
    }
	print "Total $count SN added\n";
}


sub AddLineASymble{
	my $count=0;
	my $n=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
     		my @tp= &Getheadinfo($lineIN);
     		my $tp=@tp;
     		if($n==0){
     			$n=length($tp);
     		}
     		my $l=length($count);
     		my $add=$count;
     		for(my $i=$l; $i<$_[2];$i++){
     			$add="0".$add;		
     		} 
     		for(my $i=0; $i<$tp;$i++){
     			$tp[$i]=$_[1].$add."|".$tp[$i];
     		}
     		$count++;
     		print OUT "@tp\t\n";
    	}
}

#each line add same symble for every memble
sub AddLineandMemberASymble{
	my $count=0;
	my $n=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
     		my @tp= &Getheadinfo($lineIN);
     		my $tp=@tp;
     		if($n==0){
     			$n=length($tp);
     		}
     		my $l=length($count);
     		my $add=$count;
     		for(my $i=$l; $i<$_[2];$i++){
     			$add="0".$add;		
     		} 
     		for(my $i=0; $i<$tp;$i++){
     			$tp[$i]=$_[1].$add."|".$tp[$i];
     		}
     		$count++;
     		print OUT "@tp\t\n";
    	}
}


sub Renameclass{
	#Renameclass TL 1 3
	#start from 1
	#length 3
	my $st=0;
	if($_[2]>0){
		$st=$_[2];
	}
	my $name="New";
	if(length($_[1])>0){
		$name=$_[1];
	}
	my $len=3;
	if($_[3]>0){
		$len=$_[3];
	}	
	my $count=1;

	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
		for(my $i=$st; $i<$tp; $i++){
			my $add=$count;
			my $tplen=length($count);
			for(my $j=$tplen;$j<$len;$j++){
				$add="0".$add;
			}
			print OUT "$name$add\_$i\t$tp[$i]\t$tp[0]\n";
		}		
		$count++;
	}
}


sub RenamebyTypelist{
	#start from 1
	my $count;
	my $type="NA";
	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
		if($tp[$_[1]] ne $type){
			if($type ne "NA"){
				print "$type have $count member\n";
			}
			$count=1;
			$type=$tp[$_[1]];
		}	
		my $add=$count;
		my $len=length($count);
		if ($len < $_[2]){
			for(my $j=$len; $j<$_[2]; $j++){
				$add="0".$add;
			}
		}
		print OUT "$lineIN\t$type\_$add\n";
		$count++;
	}
	print "$type have $count member\n";
}



sub RenameFamilysequence{
	my $cut="_";
	my $f4=$ARGV[0]."changedetail";
	open(OUTD,'>',"$f4");
	my $count = $_[1];
	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
		if ($lineIN=~/>/){
			my @tpf=&Getheadinfo($tp[0],$cut);	
			my $len=length($count);
			my $add=$count;
			if ($len < $_[3]){
				for(my $j=$len; $j<$_[3]; $j++){
					$add="0".$add;
				}
			}
			if($_[4]){
				$add=$add.$_[4];	
			}
			print OUT "$tpf[0]_$_[2]$add\n";
			print OUTD "$lineIN\t$tpf[0]_$_[2]$add\n";
			$count++;
		}else{
			print OUT "$lineIN\n";	
		}
	}
	print "$count sequence changed name\n";
	close(OUTD);	
}
	
sub renamesequence{
	my $f4=$ARGV[0]."changedetail";
	open(OUTD,'>',"$f4");
	my $count = $_[1];
	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
		if ($lineIN=~/>/){
			my $len=length($count);
			my $add=$count;
			if ($len < $_[3]){
				for(my $j=$len; $j<$_[3]; $j++){
					$add="0".$add;
				}
			}
			if($_[4]){
				$add=$add.$_[4];	
			}
			print OUT ">$_[2]_$add\n";
			print OUTD "$lineIN\t>$_[2]_$add\n";
			$count++;
		}else{
			print OUT "$lineIN\n";	
		}
	}
	print "$count sequence changed name\n";
	close(OUTD);	
}
		
sub GetBesInfo{
	#BAC : "a0001A12"
	#map "ctg612" Ends Left 58.000
	#map "ctg612" Ends Right 150.000 Oldctg 236
	my $o;
	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
     	if ($lineIN=~/BAC :/){
     		$o=substr($tp[2],1,(length($tp[2])-2)); 
     	}
     	elsif ($lineIN=~/Ends Left/){
     		$o=$o."\t".substr($tp[1],1,(length($tp[1])-2))."\t".$tp[4]; 	
     	}
     	elsif ($lineIN=~/Ends Right/){
     		$o=$o."\t".$tp[4]; 	
     		print OUT $o."\n";
     	}	
	}
}


sub CombinLinestoOne{
	my $count=1;
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($count<$_[1]){
			print OUT $lineIN."\t";
			$count++;
		}else{
			$count=1;
			print OUT $lineIN."\n";
		}
	}
}

sub Replace_Low_High_to_cutoff
{
	my $cut=$_[2];
	#first row is ID
	$lineIN=<IN>;
	chomp($lineIN);
	print OUT "$lineIN\n";
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp= split(/\t/,$lineIN);
		my $tp =@tp;
		print OUT "$tp[0]";
		for(my $i=1;$i<$tp;$i++){
			if($_[1] eq "gt"){
				if($tp[$i] > $cut){
					$tp[$i]=$cut;
				}
			}elsif($_[1] eq "lt"){
				if($tp[$i] < $cut){
					$tp[$i]=$cut;
				}	
			}
			print OUT "\t$tp[$i]";
		}
		print OUT "\n";
	}
}

sub ReplaceSymble1by2{
	while($lineIN=<IN>){
		chomp($lineIN);
		#change here
		$lineIN =~ s/#0\/2/.r/;
		print OUT $lineIN."\n";
	}
}


sub countFileLength{
	my $count=0;
	while($lineIN=<IN>){
		$count++;
	}
	print OUT $count."\n";
	print $count."\n";
}


sub countGCcontant{
	my $count=0;
	my $countGC=0;
	my $countN=0;
	my @c;
	my @t;
	my $tc=0;
	my $start=0;
	my $GC=0;
	my $AT=0;
	my $N=0;
	my $gap=0;
	my $other=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($lineIN =~/>/){	
			if($start==0){
				$start=1;
			}else{
				my $Totallen=$GC+$AT;
				print OUT "$GC\t$AT\t$N\t$Totallen\t$gap\t$other\n";
				$GC=0;
				$AT=0;
				$N=0;
				$gap=0;
				$other=0;
			}
			print OUT $lineIN."\t";
		}else{
			for(my $i=0;$i<length($lineIN);$i++){
				my $aa=substr($lineIN,$i,1);
				my $canadd=0;
				for(my $j=0;$j<$tc;$j++){
					if($t[$j] eq $aa){
						$c[$j]++;
						$canadd++;
						last;
					}
				}
				if($canadd == 0){
					$c[$tc]=1;
					$t[$tc]=$aa;
					$tc++;
					print "$aa";
				}
				
				if($aa eq 'G' || $aa eq 'g' || $aa eq 'C' || $aa eq 'c'){
					$countGC++;
					$count++;
					$GC++;
				}elsif($aa eq 'A' || $aa eq 'a' || $aa eq 'T' || $aa eq 't'){
					$count++;	
					$AT++;
				}elsif($aa eq 'N' || $aa eq 'n'){
					$countN++;
					$N++;	
				}elsif($aa eq '-'){
					$gap++;	
				}else{
					$other++;
				}
				
			}
		}
	}
	print OUT "$GC\t$AT\t$N\t$Totallen\t$gap\t$other\n";
	print OUT "-------------------------------------\n";	
	for (my $i=0;$i<$tc;$i++){
		print OUT "$t[$i]\t$c[$i]\n";
	}
	print OUT "-------------------------------------\n";
	my $p = 100*$countGC/$count;
	print OUT "$countGC\t$count\t$p\n";

}

sub float2int{
	#only get head before .
	my $r=$_[0];
	if(index($_[$i],".")>0){
		$r=substr($r,0,index($r,"."));
	}	
	return $r;
}	


sub CountSourtedRegionAvgValue{
	#input is sorted
	#chrom site value
	#para
	my $r=$_[1];
	my $n=0;
	my $totalvalue=0;
	my $cut=$r;
	my $count=0;
	my $avg=0;

	#print "region wideth $_[1]\n";
	#read
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($lineIN=~/#/){	
		}else{
			my @tp = split(/\t/,$lineIN);
			#my @tp = &Getheadinfo($lineIN);	

			#if($tp[0]>$cut){
			if($tp[1]>$cut){
				$count++;
				if($n>0){
					$avg=$totalvalue/$n;
				}
				print OUT "$cut\t$avg\n";					
				$cut+=$r;
				$n=1;
				#$totalvalue=$tp[1];
				$totalvalue=$tp[2];
			}else{
				$n++;
				#$totalvalue+=$tp[1];
				$totalvalue+=$tp[2];
			}
		}
	}
	$count++;
	if($n>0){
		$avg=$totalvalue/$n;
	}	
	print OUT "$cut\t$avg\n";		
	print "$count\ndone\n";

}



sub CountfixedRegionAvgValue{

	#inpute: 0-region 1-value
	#para: 1-region 10000, 2-max 100000
	
	#creat region
	#print "region wideth $_[1] \n";
	my $max=$_[2]/$_[1];
	my %value;
	my %number;
	for(my $i=0;$i<$max;$i++){
		$value{$i}=0;
		$number{$i}=0;
	}
	
	my $count=0;
	
	#read
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp = split(/\t/,$lineIN);
		#my @tp = &Getheadinfo($lineIN);	
		#print "$tp[1]\n";
		#rage
		if($tp[1] ne "NaN" && $tp[1] ne "R^2"){
			my $r = $tp[0]/$_[1];
			my $cutloc=index($r,".",);
			if($cutloc>0){
				$r=substr($r,0,$cutloc);
			}
			#print "$r @tp\n";
			#add
			$value{$r}=$value{$r}+$tp[1];
			$number{$r}=$number{$r}+1;
		}
		$count++;
	}

	
	#output
	my @hkeys=keys(%value);
	my $hkeys=@hkeys;	
	if($hkeys ne $hkeys){
		print "$hkeys ne $max\n";		
	}
	print "region number is $max\n$count useable data\n";	
	for(my $i=0;$i<$max;$i++){
		my $avg=0;
		if($number{$i}>0){
			$avg=$value{$i}/$number{$i};
		}
		print OUT "$i\t$number{$i}\t$value{$i}\t$avg\n";
	}
	
	print "done\n";

}


sub CountchromRegionTotalValue{

	#inpute: 0-chr 1-loc 2-value
	#para: 1-step, 
	#input region no redent
	
	my $step = $_[1];
	
	#creat storage
	my %data;			#2d array
	my %chrom_max;	#
	my $count=0;
	
	#read
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp = split(/\t/,$lineIN);
		#my @tp = &Getheadinfo($lineIN);	
		#print "$tp[0]\n";
		
		my $a1=int($tp[1]/$step);
		if($a1>$chrom_max{$tp[0]}){
			$chrom_max{$tp[0]}=$a1;
		}		

		my $add=$tp[2];
		if(!$data{$tp[0]}->[$a1]){
			$data{$tp[0]}->[$a1]=$add;
		}else{
			$data{$tp[0]}->[$a1]+=$add;
		}
		
		$count++;
	}
	
	print "$count lines data readed\n";
	
	#output
	my @ckeys=keys(%chrom_max);
	my $ckeys=@ckeys;	

	for(my $i=0;$i<$ckeys;$i++){
		for(my $j=0;$j<$chrom_max{$ckeys[$i]};$j++){
			if(!$data{$ckeys[$i]}->[$j]){
				$data{$ckeys[$i]}->[$j]=0;
			}
			my $st=$j*$step+1;
			my $ed=$j*$step+$step;

			print OUT "$ckeys[$i]\t$st\t$ed\t$data{$ckeys[$i]}->[$j]\n";
		}
	}
	
	print "done\n";

}


sub CountchromRegionCoverage{

	#inpute: 0-chr 1-st 2-ed
	#para: 1-step, 
	#input region no redent
	
	my $step = $_[1];
	
	#creat storage
	my %data;			#2d array
	my %chrom_max;	#
	my $count=0;
	
	#read
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp = split(/\t/,$lineIN);
		#my @tp = &Getheadinfo($lineIN);	
		#print "$tp[0]\n";
		
		my $a1=int($tp[1]/$step);
		my $a2=int($tp[2]/$step);
		if($a2>$chrom_max{$tp[0]}){
			$chrom_max{$tp[0]}=$a2;
		}
		
		if($a1==$a2){
			my $add=$tp[2]-$tp[1]+1;
			if(!$data{$tp[0]}->[$a1]){
				$data{$tp[0]}->[$a1]=$add;
			}else{
				$data{$tp[0]}->[$a1]+=$add;
			}
		}else{
			my $add=$a2*$step-$tp[1]+1;
			if(!$data{$tp[0]}->[$a1]){
				$data{$tp[0]}->[$a1]=$add;
			}else{
				$data{$tp[0]}->[$a1]+=$add;
			}	
			#	
			$add=$tp[2]-$a2*$step;
			if(!$data{$tp[0]}->[$a2]){
				$data{$tp[0]}->[$a2]=$add;
			}else{
				$data{$tp[0]}->[$a2]+=$add;
			}					
		}
		$count++;
	}
	
	print "$count lines data readed\n";
	
	#output
	my @ckeys=keys(%chrom_max);
	my $ckeys=@ckeys;	

	for(my $i=0;$i<$ckeys;$i++){
		for(my $j=0;$j<$chrom_max{$ckeys[$i]};$j++){
			if(!$data{$ckeys[$i]}->[$j]){
				$data{$ckeys[$i]}->[$j]=0;
			}
			my $st=$j*$step+1;
			my $ed=$j*$step+$step;

			print OUT "$ckeys[$i]\t$st\t$ed\t$data{$ckeys[$i]}->[$j]\n";
		}
	}
	
	print "done\n";

}


sub CountTypefixedRegionNumber{
	#input is sorted
	#para
	my $r=$_[1];
	my $type="Type";
	my $n=0;
	my $cut=$r;
	my $count=0;

	#print "region wideth $_[1]\n";
	#read
	while($lineIN=<IN>){
		chomp($lineIN);	
		if($lineIN=~/#/){	
		}else{
			my @tp = split(/\t/,$lineIN);
			#my @tp = &Getheadinfo($lineIN);	
			if($tp[0] eq $type){
				if($tp[1]>$cut){
					$count++;
					my $st=$cut-$r+1;
					print OUT "$type\t$st\t$cut\t$n\n";				
					$cut+=$r;
					$n=1;
				}else{
					$n++;
				}
			}else{
				if($type!="Type"){
					$count++;
					my $st=$cut-$r+1;
					print OUT "$type\t$st\t$cut\t$n\n";		
				}			
				$type=$tp[0];
				print "$type\n";
				$cut=$r;
				$n=1;			
			}	
		}
	}
	$count++;
	my $st=$cut-$r+1;
	print OUT "$type\t$st\t$cut\t$n\n";		
	print "\ndone\n";

}


sub CountfixedRegionNumber{
	#para
	my $n=@_;
	$n=$n-2;
	my @row;
	for(my $i=0;$i<$n;$i++){
		$row[$i]=$_[$i+2];
	}
	
	my @array; #2d
	my $ac=0;
	my $max=0;
	print "region wideth $_[1]\n";
	#read
	while($lineIN=<IN>){
		chomp($lineIN);	
		#my @tp = split(/\t/,$lineIN);
		my @tp = &Getheadinfo($lineIN);	
		for(my $i=0;$i<$n;$i++){
			$array[$ac][$i]=$tp[$row[$i]];
			#print "$tp[$row[$i]];";
			if($array[$ac][$i]>$max){
				$max=$array[$ac][$i];
			}
		}
		$ac++;		
	}
	#@array=sort {$a <=> $b} @array;
	print "$ac line and $n row readed\n";
	#init
	my @Rcount; #2d	
	my $rc=&float2int($max/$_[1]);

	for(my $i=0;$i<$rc;$i++){
		for(my $j=0;$j<$n;$j++){
			$Rcount[$i][$j]=0;
		}
	}
	#add
	for(my $i=0;$i<$ac;$i++){
		for(my $j=0;$j<$n;$j++){
			my $tp=&float2int($array[$i][$j]/$_[1]);
			$Rcount[$tp][$j]++;
			
		}
	}
	#output
	for(my $i=0;$i<$rc;$i++){
		my $st=$i*$_[1];
		my $ed=($i+1)*$_[1];
		print OUT "$st_$ed";
		for(my $j=0;$j<$n;$j++){			
			print OUT "\t$Rcount[$i][$j]";
		}
		print OUT "\n";
	}
	print "$rc region done\n";
}


sub CountRegionDetail{
	my @array;
	my @r;		#region
	my @rl;		#region length
	my @rn;		#region number
	my $tn=0;   #total number
	my $tl=0;	#total length
	my $region = @_;
	$_=@_;

	for (my $i=1; $i<$_; $i++){
		$r[$i-1]=$_[$i];
		$rn[$i-1]=0;
		$rl[$i-1]=0;
	}
	my $r=@r;
	while($lineIN=<IN>){
		chomp($lineIN);	
		$array[$tn]=$lineIN;
		$tn++;
		$tl+=$lineIN;
	}
	
	print "total $tn number seperate to $r region";
	
	for(my $i=0;$i<$tn;$i++)
	{	
		for(my $j=0;$j<$r;$j++){				
			if($array[$i] >= $r[$j] && $array[$i] < $r[$j+1]){
				$rn[$j]++;
				$rl[$j]+=$array[$i];
				last;
			}
		}	
	}	
	
	#output
	print OUT "Region\tNumber\tLength\tAvgLength\tPercent\n";
	print "Region\tNumber\tLength\tAvgLength\tPercent\n";
	for(my $i=0;$i<$r;$i++){
		my $avg=0;
		if($rn[$i]>0){$avg=$rl[$i]/$rn[$i];}
		my $p=$rl[$i]*100/$tl;
		print OUT "$r[$i]-$r[$i+1]\t$rn[$i]\t$rl[$i]\t$avg\t$p\n";
		print "$r[$i]-$r[$i+1]\t$rn[$i]\t$rl[$i]\t$avg\t$p\n";
	}
}

	
sub CounttypeRegionNumberAndTotalvalue{
	my %arraytype;  #key-type-0 value-1 length-2
	my @arrayno; #2d
	my @arrayvalue; #2d
	my $max=0;
	my $typecount=0;
	
	my $region = 10;
	if($_[1]>0){
		$region=$_[1];
	}

	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		$max=$max>$tp[1]?$max:$tp[1];
		
		if($arraytype{$tp[0]}>0){
		}else{
			$typecount++;
			$arraytype{$tp[0]}=$typecount;
			print "$tp[0] $typecount\n";			
		}
		my $tptype=$arraytype{$tp[0]};
		
		#add to array
		if($arrayno[$tptype][$tp[1]]>0){
			$arrayno[$tptype][$tp[1]]++;
			$arrayvalue[$tptype][$tp[1]]+=$tp[2];
		}else{
			$arrayno[$tptype][$tp[1]]=1;
			$arrayvalue[$tptype][$tp[1]]=$tp[2];
		}		
	}
	print "\n$typecount types and $max detail\n\nnnumber count ...\n\n";
	
	#output	
	my @hashkeys=keys(%arraytype);
	@hashkeys=sort(@hashkeys);

	print OUT "NumberCount";
	for(my $j=1;$j<=$max;$j++){
		print OUT "\t$j";
	}	
	for(my $i=0;$i<$typecount;$i++){
		print OUT "\n$hashkeys[$i]";
		for(my $j=1;$j<=$max;$j++){
			if($arrayno[$arraytype{$hashkeys[$i]}][$j]>0){
				print OUT "\t$arrayno[$arraytype{$hashkeys[$i]}][$j]";
			}else{
				print OUT "\t0";
			}
		}
	}

	print "value ... \n";
	print OUT "\n\nValueCount";
	for(my $j=1;$j<=$max;$j++){
		print OUT "\t$j";
	}	
	for(my $i=0;$i<$typecount;$i++){
		print OUT "\n$hashkeys[$i]";
		for(my $j=1;$j<=$max;$j++){
			if($arrayno[$arraytype{$hashkeys[$i]}][$j]>0){
				print OUT "\t$arrayvalue[$arraytype{$hashkeys[$i]}][$j]";
			}else{
				print OUT "\t0";
			}
		}
	}	
	
	print OUT "\n\nCombinNumber";
	for(my $j=1;$j<=$max;$j+=$region){
		my $pp=$j+$region-1;	
		print OUT "\t$j..$pp";
	}
	for(my $i=0;$i<$typecount;$i++){
		print OUT "\n$hashkeys[$i]";
		for(my $j=1;$j<=$max;$j+=$region){
			my $tv=0;
			for(my $k=0;$k<$region;$k++){
				my $v=0;
				if($arrayno[$arraytype{$hashkeys[$i]}][$j+$k]>0){
					$v = $arrayno[$arraytype{$hashkeys[$i]}][$j+$k];
				}
				$tv+=$v;
			}
			print OUT "\t$tv";
		}
	}
	print OUT "\n\nCombinValue";
	for(my $j=1;$j<=$max;$j+=$region){
		my $pp=$j+$region-1;	
		print OUT "\t$j..$pp";
	}
	for(my $i=0;$i<$typecount;$i++){
		print OUT "\n$hashkeys[$i]";
		for(my $j=1;$j<=$max;$j+=$region){
			my $tv=0;
			for(my $k=0;$k<$region;$k++){
				my $v=0;
				if($arrayno[$arraytype{$hashkeys[$i]}][$j+$k]>0){
					$v = $arrayvalue[$arraytype{$hashkeys[$i]}][$j+$k];
				}
				$tv+=$v;
			}
			print OUT "\t$tv";
		}
	}		
		
}		
	
sub CountRegionNumber{
	my @array;
	my $ac=0;
	my @r;
	my @rl; 
	my $region = @_[1];
	my $first=1;

	for (my $i=0; $i<=$_; $i++){
		$r[$i]=0;
		$rl[$i]=0;
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		$array[$ac]=$lineIN;
		$ac++;
	}
	@array=sort {$a <=> $b} @array;
	my $rc=0; #type
	my $nc=0; #number
	for(my $i=0;$i<$ac;$i++)
	{	
		if($array[$i]<$rc*$region){
			$nc++;
			print 
		}else{
			if($first==1){
				$first=2;
			}else{
				print OUT "$rc\t$nc\n";
				print "$rc\t$nc\n";
			}
			$rc++;
			$nc=1;
		}		
	}	
}	

	
sub CountRegionLength{
	#musted sorted
	#input
	my @array;
	my $count=0;
	my $cut=0;
	my $len=0;
	my $nu=0;  
	my $avg=0;
	my $region = @_[1];	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		if($tp[0]>$cut){
			while($tp[0]>$cut){
				$count++;  
				if($nu==0){
					$avg=0;
				}else{
					$avg=$len/$nu;
				}
				print OUT "$cut\t$nu\t$len\t$avg\n";
				print "$cut\t$nu\t$len\t$avg\n";
				$nu=0;
				$len=0;
				$cut=$count*$region;
			}
		}
		$len+=$tp[0];
		$nu++;		
	} 
				if($nu==0){
		$avg=0;
	}else{
		$avg=$len/$nu;
	}
	print OUT "$cut\t$nu\t$len\t$avg\n";
	print "$cut\t$nu\t$len\t$avg\n";
}

sub CombinbyLine{
	my $count = 0;
	my $combin = 0;
	while($lineIN=<IN>){
		chomp($lineIN);	
     	if ($count < $_[1]){
     		$combin+=$lineIN;
     		$count++;
     	}else{
			print OUT $combin."\n";
			$count=1;
			$combin=$lineIN;     	
     	}
	}
	print OUT $combin."\n";	
}	

	
sub AddLosstoLine{
	my $count = 0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp= &Getheadinfo($lineIN);	
		if ($count == 0){
			$count = $tp[0];
			print OUT $count."\t".$tp[1]."\n";
		}else
		{
			$count+=$_[1];
			if ($count < $tp[0]){
				while($count < $tp[0]){
					print OUT $count."\t".$_[2]."\n";
					$count+=$_[1];
				}	
			}
			print OUT $count."\t".$tp[1]."\n";
		}
	}
	print OUT $combin."\n";	
}
		
sub CombinRepeatmaskerResult{
	#combin record form littletool
	#come same family record in Certain Scope\
	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
	}
}
		
		
sub RreturnLocTorows{
	while($lineIN=<IN>){
		chomp($lineIN);	
     	my @tp= &Getheadinfo($lineIN);		
		my $tp=@tp;
		for (my $i=0;$i<$tp;$i++){
			if ($tp[$i] =~/:/){
				@w = split(/:/,$tp[$i]);
				print OUT "$w[0]\t";
				my $p= index($w[1],'..',);
				if ($w[1] =~/complement/){					
					my $r=substr($w[1],11,$p-11)."\t".substr($w[1],$p+2,length($w[1])-$p-3);
					print OUT $r;
				}else{
					my $r=substr($w[1],0,$p)."\t".substr($w[1],$p+2,length($w[1])-$p-2);
					print OUT $r;					
				}
			}else{
				print OUT "$tp[$i]\t";
			}			
		}
		print OUT "\n";
	}
}
		
		
sub BlatSearchFliter{
	#1,hit length >original lentth X 80% £¨0+1 > 10*0.8£©
	#2,for each hit,calculate total length by start and end distance,15-11 to 16+(10-12)
	#3,hit total length <(9)original lentth X 2
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $scop=$tp[0]+$tp[1];
		if ($scop>($tp[10]*$_[1])){
			my $l=($tp[16]+$tp[10]-$tp[12])-($tp[15]-$tp[11]);
			my $a =$tp[10]*$_[2];
			if ($l<($tp[10]*$_[2])){
				print OUT $lineIN."\n";
			}		
		}		
	}			
}



sub FindAllORF{
	#read in record,translate from 1,2,3,-1,-2,-3
	#count stop codon number,find best translation
	my $name='start';
	my $seq;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN =~/>/){
			if($name ne 'start'){
				my @o=&GetallORF($seq);
				for(my $i=0;$i<6;$i++){
					print OUT "$name=$i\n$o[$i]\n";
				}
			}
			$name = $lineIN;
			$seq = '';
		}else{
			$seq=$seq.$lineIN;
		}
	
	}
	my @o=&GetallORF($seq);	
	for(my $i=0;$i<6;$i++){
		print OUT "$name=$i\n$o[$i]\n";
	}		
}


sub GetallORF{
	my @stopcount;
	my @transleted;
	my $min=9999;
	my $ORF;
	my $minseq;
	my $minloc;
	my @o;
	my $count=0;
	#forward 123
	for ($i=0;$i<3;$i++){
		my $seq=substr($_[0],$i,);
		$o[$count]=&TranslateSeq($seq);
		$count++;
	}
	#reverse 123
	my $c = &Getcomplementseq($_[0]);
	for ($i=0;$i<3;$i++){
		my $seq=substr($c,$i,);
		&TranslateSeq($seq);
		$o[$count]=&TranslateSeq($seq);
		$count++;
	}
	return @o;
}


sub FindbestORF{
	#read in record,translate from 1,2,3,-1,-2,-3
	#count stop codon number,find best translation
	my $name='start';
	my $seq;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN =~/>/){
			if($name ne 'start'){
				print OUT $name;
				my $o=&GetbestORF($seq,$name);
				print OUT $o;
			}
			$name = $lineIN;
			$seq = '';
		}else{
			#remove -
			$lineIN =~tr/-//d;
			#print "$lineIN\n";
			$seq=$seq.$lineIN;
		}
	
	}
	print OUT $name;
	#&GetbestORF($seq,$name);
	my $o=&GetbestORF($seq,$name);
	print OUT $o;			
}
			
sub compairstringrow{
	#print "$_[0] $_[1] $_[1]\n";
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if ($_[2] eq "eq"){
			if ($tp[$_[1]] eq $tp[$_[3]] ){
				print OUT $lineIN."\n";
			}
		}
		if ($_[2] eq "ne"){
			if ($tp[$_[1]] ne $tp[$_[3]] ){
				print OUT $lineIN."\n";
			}
		}
	}			
}

	# compair int or float row
sub compairintrow{
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		my $add=0;
		$_=@_;
		for(my $i=3;$i<$_;$i++){
			if ($_[2] eq "eq"){
				if ($tp[$_[1]] == $tp[$_[$i]] ){
					$add=1;
				}				
			}
			if ($_[2] eq "gt"){
				if ($tp[$_[1]] > $tp[$_[$i]] ){
					$add=1;	
				}
		}
			if ($_[2] eq "lt"){
				if ($tp[$_[1]] < $tp[$_[$i]] ){
					print "$tp[$_[1]] $tp[$_[$i]]\n";
					$add=1;	
				}
			}
			if ($_[2] eq "ne"){
				if ($tp[$_[1]] != $tp[$_[$i]] ){
					$add=1;
				}				
			}			
		}
		if($add==1){
			print OUT $lineIN."\n";
		}
	}
}
		
sub GetbigSmallone{
	my $wh="big";
	if ($_[1]){
		$wh=$_[1];
	}
	print "$wh...\n";
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		if ($wh eq 'big'){
			
			my @sorted = sort{$b <=> $a } @tp;
			print OUT $sorted[0]."\t".$lineIN."\n";

		}
		elsif($wh eq 'small'){
			my @sorted = sort{$a <=> $b } @tp;
			print OUT $sorted[0]."\t".$lineIN."\n";
		}
	}			
}	


#SelectBuyAtleastOneRowValue
sub SelectBuyAtleastOneRowValue{
	#para
	$_=@_;
	my @row;
	my $count=0;
	my $dowhat="";
	my $cut;
	for(my $i=1;$i<$_;$i++){
		if($_[$i] eq "gt" || $_[$i] eq "lt" || $_[$i] eq "eq" || $_[$i] eq "ne" || $_[$i] eq "=" || $_[$i] eq "!="){
			$dowhat=$_[$i];
			$cut=$_[$i+1];
			last;
		}else{
			$row[$count]=$_[$i];
			$count++;
		}
	}
	print "$count dowhat=$dowhat cut=$cut\n";
	#read
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		my $add=0;
		for(my $i=0;$i<$count;$i++){
		
			if ($dowhat eq "gt" ){
				#print "a";
				if ($tp[$row[$i]] > $cut){	
					$add = 1; 
					last;				
				}
			}
			elsif ($dowhat eq "=" ){
				if ($tp[$row[$i]] == $cut){
					$add = 1; 
					last;						
				}
			}
			elsif ($dowhat eq "!=" ){
				if ($tp[$row[$i]] != $cut){
					$add = 1; 
					last;						
				}
			}
			elsif ($dowhat eq "lt" ){
				if ($tp[$row[$i]] < $cut){	
					$add = 1; 
					last;										
				}
			}
			elsif ($dowhat eq "eq" ){
				if ($tp[$row[$i]] eq $cut){
					$add = 1; 
					last;						
				}
			}
			elsif ($dowhat eq "ne" ){
				if ($tp[$row[$i]] ne $cut){
					$add = 1; 
					last;						
				}
			}
		}
		if($add ==1 ){
			print OUT $lineIN."\n";
		}
	}

}


#SelectBuyRowValue
sub SelectBuyRowValue{
	#para
	$_=@_;
	my @row;
	my $count=0;
	my $dowhat="";
	my $cut;
	for(my $i=1;$i<$_;$i++){
		if($_[$i] eq "gt" || $_[$i] eq "lt" || $_[$i] eq "eq" || $_[$i] eq "ne" || $_[$i] eq "=" || $_[$i] eq "!="){
			$dowhat=$_[$i];
			$cut=$_[$i+1];
			last;
		}else{
			$row[$count]=$_[$i];
			$count++;
		}
	}
	print "$count dowhat=$dowhat cut=$cut\n";
	#read
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		#my @tp = split(/\t/,$lineIN);
		my $add=1;
		for(my $i=0;$i<$count;$i++){
		
			if ($dowhat eq "gt" ){
				#print "a";
				if ($tp[$row[$i]] >= $cut){					
				}else{ $add = 0; }
			}
			elsif ($dowhat eq "=" ){
				if ($tp[$row[$i]] == $cut){
				}else{ $add = 0; }
			}
			elsif ($dowhat eq "!=" ){
				if ($tp[$row[$i]] != $cut){
				}else{ $add = 0; }
			}
			elsif ($dowhat eq "lt" ){
				if ($tp[$row[$i]] <= $cut){					
				}else{ $add = 0; }
			}
			elsif ($dowhat eq "eq" ){
				if ($tp[$row[$i]] eq $cut){
				}else{ $add = 0; }
			}
			elsif ($dowhat eq "ne" ){
				if ($tp[$row[$i]] ne $cut){
				}else{ $add = 0; }
			}
		}
		if($add ==1 ){
			print OUT $lineIN."\n";
		}
	}

}

sub SelectbyRowString{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if ($_[2] eq "ne" ){
			if ($tp[$_[1]] ne $_[3]){
				print OUT $lineIN."\n";
			}
		}	
		elsif ($_[2] eq "eq" ){
			if ($tp[$_[1]] eq $_[3]){
				print OUT $lineIN."\n";
			}
		}
	}
}

#select buy row String length
sub rowlengthcutoff{
	printf("ok\n");
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		$tplen = length($tp[$rowcutoff[1]]);
		if ($_[2] eq ">" ){
			if ( $tplen > $_[3]){
				print OUT $tplen."\t".$lineIN."\n";
			}
		}	
		elsif ($_[2] eq "=" ){
			if ($tplen == $_[3]){
				print OUT $tplen."\t".$lineIN."\n";
			}
		}
		elsif ($_[2] eq "<" ){
			if ($tplen < $_[3]){
				print OUT $tplen."\t".$lineIN."\n";
			}
		}
		elsif ($_[2] eq "<=" ){
			if ($tplen <= $_[3]){
				print OUT $tplen."\t".$lineIN."\n";
			}
		}	
		elsif ($_[2] eq ">=" ){
			if ($tplen >= $_[3]){
				print OUT $tplen."\t".$lineIN."\n";
			}
		}					
	}
}
	

sub Counttypecontinuationno{
	my $name="NA";
	my $count=0;
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		if($tp[$_[1]] eq $name){
			$count++;
		}else{
			if($name ne "NA"){
				print OUT "$name\t$count\n";
			}	
		$name=$tp[$_[1]];
			$count=1;
		}
	}
	print OUT "$name\t$count\n";
}


sub CountCombintypeValue{
	#input
	#name1	name2	value
	
	my @name;	#0-combin type 1-no 2-value
	my $count=0;
	
	while($lineIN=<IN>){
		chomp($lineIN);
		#@tp = &Getheadinfo($lineIN);
		my @tp=split(/\t/,$lineIN);
		my $tpname=$tp[$_[1]]."\t".$tp[$_[2]];
		my $add=1;
		for(my $i=0;$i<$count;$i++){
			if($name[$i][0] eq $tpname){
				$name[$i][1]++;
				$name[$i][2]+=$tp[$_[3]];
				$add=0;
				last;
			}
		}
		if($add==1){
			$name[$count][0] = $tpname;
			$name[$count][1] = 1;
			$name[$count][2] = $tp[$_[3]];
			$count++;
		}

	}
	#output
	for(my $i=0;$i<$count;$i++){
		print OUT "$name[$i][0]\t$name[$i][1]\t$name[$i][2]\n";
	}
	print "done\n";
}

sub CountCombintypeno{
	#input
	#name1	name2	name3
	#value1	value2	value3
	
	my @info;	#2ed
	my $infono=0;
	my @name;	#combin type
	my $nameno=0;
	my @type;
	#check combin
	$lineIN=<IN>;
	chomp($lineIN);
	my @tp = &Getheadinfo($lineIN);	
	my $tp=@tp;

	for(my $i=0;$i<$tp;$i++){
		for(my $j=$i+1;$j<$tp;$j++){
			$type[$infono] = "$tp[$i]=$tp[$j]";
			$infono++;
		}
	}
	print "infono=$infono\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);
		@tp = &Getheadinfo($lineIN);
		my $n=0;	
		for(my $i=0;$i<$tp;$i++){
			for(my $j=$i+1;$j<$tp;$j++){
				my $aa=$tp[$i].$tp[$j];
				my $addn=-1;
				for(my $k=0;$k<$nameno;$k++){
					if($name[$k] eq $aa){
						$info[$k][$n]++;
						$addn=1;
					}
				}
				if($addn == -1){
					#new
					$name[$nameno]=$aa;	
					#print "$aa\n";				
					for($l=0;$l<$infono;$l++){
						$info[$nameno][$l]=0;
					}
					$info[$nameno][$n]=1;
					$nameno++;
				}
				$n++;
			}
			
		}
	}
	#export result
	print OUT "Type";
	for(my $j=0;$j<$nameno;$j++){
		print OUT "\t$name[$j]";		
	}	
	print OUT "\n";
	for(my $i=0;$i<$infono;$i++){
		print OUT "$type[$i]";
		for(my $j=0;$j<$nameno;$j++){
			print OUT "\t$info[$j][$i]";	
		}			
		print OUT "\n";
	}
	print "done\n";
}


#count type no
sub Counttypeno{
	
	while($lineIN=<IN>){
		chomp($lineIN);
		#my @tp = split(/\t/,$lineDB);
		my @tp = &Getheadinfo($lineIN);	
		if (!$db{$tp[$_[1]]}){
			$db{$tp[$_[1]]} = 1;		
		}else{
			$db{$tp[$_[1]]}++;			
		}
	}
	#export result
	while(($type, $num) = each(%db)){
		print OUT $type."\t".$num."\n";
	}
}



sub CountlineLength{
	while($lineIN=<IN>){
		chomp($lineIN);
		my $len=length($lineIN);
		print OUT "$lineIN\t$len\n";
	}

}

sub CountlineTypeNo{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $tp=@tp;
		%db=();
		for (my $i=0; $i<$tp; $i++){
			if (!$db{$tp[$i]}){
				$db{$tp[$i]} = 1;		
			}else{
				$db{$tp[$i]}++;			
			}
		}
		#export result
		my @o=keys(%db);
		my $o=@o;
		print OUT "$tp\t$o\t@o\n";
	}

}

sub CountlineFamilyNo{
	my $sb = "_";

	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $tp=@tp;
		%db=();
		for (my $i=0; $i<$tp; $i++){
			my @f= &Getheadinfo($tp[$i], $sb);	
			if (!$db{$f[0]}){
				$db{$f[0]} = 1;		
			}else{
				$db{$f[0]}++;			
			}
		}
		#export result
		my @o=keys(%db);
		my $o=@o;
		
		print OUT "$tp\t$o\t@o\n";
	}

}


sub CounttypeTotalValue{
	#in
	#0-type 1-v1 2-v2 3-v3
	#@db save value 0- type 1-no 2-3-4-value
	my $n=0;
	my $vn=$_;
	my @para=@_;
	my $para=@para;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $add=-1;
		for(my $i=0;$i<$n;$i++){
			if($db[$i][0] eq $tp[$para[1]]){
				$add=$i;
				last;
			}
		}
		if($add==-1){   #add new type
			$db[$n][0]=$tp[$para[1]];
			$db[$n][1]=1;
			for(my $j=2;$j<$para;$j++){
				$db[$n][$j]=$tp[$para[$j]];
			}
			$n++;
			#print ".";
		}else{
			$db[$add][1]++;
			for(my $j=2;$j<$para;$j++){
				$db[$add][$j]+=$tp[$para[$j]];
			}
		}
	}
	print "$n type $para value\n";
	print OUT "Type\tNo.";
	for(my $j=2;$j<$para;$j++){
		if($i==0){
			print OUT "\tcol:$para[$j]";
		}
	}	
	print OUT "\n";
	
	for(my $i=0;$i<$n;$i++){
		print OUT "$db[$i][0]\t$db[$i][1]";
		for(my $j=2;$j<$para;$j++){
			print OUT "\t$db[$i][$j]";
		}
		print OUT "\n";
	}
}


sub Counttype_nonredundant_Length{
	#input must sorted
	my $st=0;
	my $ed=0;
	my $loop=0;
	my $count=0;
	my @data=();
	my $data=0;
	my $n=99;

	#read in data
	while($lineIN=<IN>){
		#my @tp = &Getheadinfo($lineIN);	
		my @tp= split(/\t/,$lineIN);
		$tp=@tp;
		$data[$count][0]=$tp[$_[1]];
		$data[$count][1]=$tp[$_[2]];
		$data[$count][2]=$tp[$_[3]];
		$count++;
	}
	print "Input $count data\n";

	#loop combin
	while($n>0){
		$n=0;
		$loop++;
		my @ndata=();
		my $ncount=0;
		$data=@data;
		#print "data=$data\n";

		print "Loop $loop\t";
		for(my $i=1;$i<$data;$i++){
			if($data[$i][0] eq $data[$i-1][0]){
				#new data
				if($data[$i][1] > $data[$i-1][2]){
					#print "+";
					$ndata[$ncount][0]=$data[$i-1][0];
					$ndata[$ncount][1]=$data[$i-1][1];
					$ndata[$ncount][2]=$data[$i-1][2];
					$ncount++;
				}else{
				#connect
					$n++;
					#print "=";
					$data[$i][1]=$data[$i-1][1];
					$data[$i][2]=$data[$i][2]>$data[$i-1][2]?$data[$i][2]:$data[$i-1][2];
				}
			}else{
			#new seq
				$ndata[$ncount][0]=$data[$i-1][0];
				$ndata[$ncount][1]=$data[$i-1][1];
				$ndata[$ncount][2]=$data[$i-1][2];
				#print "+";
				$ncount++;
			}
		}
		$ndata[$ncount][0]=$data[$data-1][0];
		$ndata[$ncount][1]=$data[$data-1][1];
		$ndata[$ncount][2]=$data[$data-1][2];
		#print "+";
		$ncount++;
		@data=@ndata;
		print "$n, Remain $ncount data\n";
	}

	#count total length
	$data=@data;
	#print "Remain $data data ";
	my $len=0;
	my $name="NA";
	my $totallen=0;
	for(my $i=0;$i<$data;$i++){	
		#print "$data[$i][0]\n";
		my $tplen=$data[$i][2]-$data[$i][1]+1;
		if($name eq "NA"){
			$name = $data[$i][0];
		}
		if($name ne $data[$i][0]){
			print ".";
			print OUT "$name\t$len\n";
			$name = $data[$i][0];
			$len=$tplen;
		}else{
			$len+=$tplen;
		}
		$totallen+=$tplen;
	}
	print OUT "$name\t$len\ntotal\t$totallen\n";
	print "\ndone\n";
}


sub CounttypeAverageValue{
	my $name;
	my $count=0;
	my @value;
	my $n=@_;
	$n=$n-2;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if($tp[$_[1]] eq $name){
			$count++;
			for(my $i=0;$i<$n;$i++){
				$value[$i]+=$tp[$_[$i+2]];
			}
		}else{
			if($count>0){
				print OUT "$name\t$count";
				for(my $i=0;$i<$n;$i++){
					my $e=$value[$i]/$count;
					print OUT "\t$e";
				}
				print OUT "\n";
			}
			$name=$tp[$_[1]];
			$count=1;
			for(my $i=0;$i<$n;$i++){
				$value[$i]=$tp[$_[$i+2]];
			}
		}	
	}
	print OUT "$name\t$count";
	for(my $i=0;$i<$n;$i++){
		my $e=$value[$i]/$count;
		print OUT "\t$e";
	}
	print OUT "\n";
}


sub ListsubtypeMember{
	#sorted
	my $type="NA";
	my $subtype;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if ($type ne $tp[$_[1]]){
			if($type ne "NA"){
				print OUT "$type\t$subtype\n";
			}
			$type=$tp[$_[1]];
			$subtype=$tp[$_[2]];
			$count=1;
		}else{
			$subtype=$subtype.",".$tp[$_[2]];
			$count++;
		}
	}
	print OUT "$type\t$subtype\n";
	
}
		
#count Type length by Type row & length row
sub counttypelength{
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		if (!$db{$tp[$_[1]]}){
		$db{$tp[$_[1]]} = $tp[$_[2]];		
		}else{
			$db{$tp[$_[1]]}+=$tp[$_[2]];			
		}
	}
	#export result
	while(($type, $num) = each(%db)){
		print OUT $type."\t".$num."\n";
	}		
}

sub getAllSeqAfterName{
	
	my $canadd=0;
	my $count=0;
	my $name=$_[1];
	
	my $canadd=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/$name/){
			$canadd=1;
		}
		if($canadd==1){
			print OUT $lineIN."\n";
		}
	}
}	


sub getseqwithoutsymble{
	
	my $sb=$_[1];
	my $canadd=0;
	my $canadd=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			$canadd=1;
			if(index($lineIN,$sb,)>0){
				$canadd=0;
			}
		}
		if($canadd==1){
			print OUT "$lineIN\n";
		}
			
	}
}

sub getseqfromsymble{
	
	my $sb=$_[1];
	my $canadd=0;
	my $canadd=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			$canadd=0;
			if(index($lineIN,$sb,)>0){
				$canadd=1;
			}
		}
		if($canadd==1){
			print OUT "$lineIN\n";
		}
			
	}
}	

#get seq have both side hit
#group 0 the same value, group 1 is L and R?
#group 3 the same value, group 4 is L and R?
sub getseqwithbothside{
	my @side1;
	my @side1;
	my $name1 = '';
	my $name2 = '';
	my $count = 0;	
	while($lineIN=<IN>){

		my @tp = &Getheadinfo($lineIN);
		if ($name1 eq $tp[0] && $name2 eq $tp[2]){
			$count++;
			@side1[$count] = $tp[1];	
			@side2[$count] = $tp[3];		
		}else{			
			if (($count == 1) && ($side1[0] ne $side1[1]) && ($side2[0] ne $side2[1])){
				print OUT $name1."\t".$name2."\n";		
			}
			$name1 = $tp[0];		
			$name2 = $tp[2];	
			$count = 0;	
			@side1[0] = $tp[1];	
			@side2[0] = $tp[3];
		}		
	}
	if (($count == 1) && ($side1[0] ne $side1[1]) && ($side2[0] ne $side2[1])){
		print OUT $name1."\t".$name2."\n";	
	}
}


sub GetRowtoLinesfromNumber{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		$_= @_;
		$tpout = $tp[$_[1]];
		for ($i = 2; $i < $_; $i++){
			$tpout = $tpout."\n".$tp[$_[$i]];		
		}
		print OUT $tpout."\n";		
	}		
}

sub GetfasNamefromNumber{
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			my @tp = &Getheadinfo($lineIN,'>');
			$_= @_;
			$tpout = $tp[$_[1]];
			for ($i = 2; $i < $_; $i++){
				$tpout = $tpout."\t".$tp[$_[$i]];		
			}
			print OUT ">$tpout\n";	
		}else{
			print OUT $lineIN."\n";
		}	
	}		
}


sub GetRowMatrixfromFilelist_sort_muti{
	
	#key = 0 
	#value local $_[2] $_[3] .....
	$n= @_;
	@loc=@_;
	
	#data saved in 
	my @matrix;
	my $count=0;
	my $add=1;
	my $maxlen=0;
	###
	$matrix[0][0]="ID";
	print "n=$n\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);
		for ($i = 1; $i < $n; $i++){
			$matrix[0][$count*($n-1)+$i]="$lineIN.$i";
		}		
		my $len=0;
		print "$lineIN\n";
		
		open(TPIN,'<',"$lineIN");		
		while($lineTPIN=<TPIN>){
			chomp($lineTPIN);
			my @tp = &Getheadinfo($lineTPIN);	
			#file len
			$len++;
			#add ID
			if($add==1){
				$matrix[$len][0]=$tp[0];
			}
			#check
			if($tp[0] eq $matrix[$len][0]){ 	
				$matrix[$len][0]=$tp[0];			
				for ($i = 1; $i < $n; $i++){
					if($tp[$loc[$i]]<0){
						print "loc=$loc[$i] line=$lineTPIN\n";
					}
					$matrix[$len][$count*($n-1)+$i]=$tp[$loc[$i]];						
				}					
			}else{
				print "ID name diff at $lineIN : $tp[0] ne $matrix[$len][0]\n";
			}
		}
		
		if($add==1){
				$maxlen=$len;
				$add=0;		
		}
		#check length
		if($len ne $maxlen){
			print "file length diff at $lineIN: $len ne $maxlen\n";
		}	

		close(TPIN);	
		$count++;
	}	
	
	#print OUT
	my $col=$count*($n-1)+1;
	print "col=$col line=$maxlen\n";
	
	for ($i = 0; $i < $maxlen; $i++){	
		print OUT "$matrix[$i][0]";
		for ($j = 1; $j < $col; $j++){	
			print OUT "\t$matrix[$i][$j]";
		}
		print OUT "\n";
	}
}

sub GetRowMatrixfromFilelist{
	
	#key = 0 
	$loc=1;
	if ($_[1]>0){
		$loc=$_[1];
	}
	
	#data saved in 
	my @matrix;  #2d
	my $count=0;
	my $len=0;
	my $add=1;
	###
	my $tittle= "ID";
	
	while($lineIN=<IN>){
		chomp($lineIN);
		$tittle=$tittle."\t$lineIN";
		$count++;
		open(TPIN,'<',"$lineIN");		
		while($lineTPIN=<TPIN>){
			chomp($lineTPIN);
			my @tp = &Getheadinfo($lineTPIN);	
			
			#add ID
			if($add==1){				
				$matrix[$len][0]=$tp[0];		
				$matrix[$len][1]=$tp[$loc];		
				$len++;
			}else{
				my $canadd=0;
				for ($i = 0; $i < $len; $i++){
					if($tp[0] eq $matrix[$i][0]){ 
						$canadd=1;	
						$matrix[$i][$count]=$tp[$loc];
						last;		
					}	
				}
				if($canadd==0){
					print "file length diff in $lineIN at $lineTPIN\n";
				}		
			}						
		}		
		if($add==1){
			print "file length = $len\n";
			$add=0;		
		}

		close(TPIN);	
	}	
	print "file no = $count and file length = $len\n";
	#print OUT
	print OUT "$tittle\n";
	
	for ($i = 0; $i < $len; $i++){	
		print OUT "$matrix[$i][0]";
		for ($j = 1; $j <=$count; $j++){	
			print OUT "\t$matrix[$i][$j]";
		}
		print OUT "\n";
	}
}



sub GetRowWithoutNumber{
	
	my %rd;
	$_= @_;
	print "remove";
	for (my $i = 1; $i < $_; $i++){
		$rd{$_[$i]} = 1;		
		print " $_[$i]";
	}	
	print "\n";

	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		#my @tp=split(/\t/,$lineIN);
		my $tp=@tp;
		$tpout = "";
		my $first=1;
		for (my $i = 0; $i < $tp; $i++){
			if($rd{$i} == 1){
				print ".";
			}else{
				if($first==1){
					$tpout =$tp[$i];
					$first=0;
				}else{
					$tpout = $tpout."\t".$tp[$i];	
				}	
			}
		}
		print OUT $tpout."\n";		
	}		
}


sub GetRowfromNumber_splite{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN,">",":","|","=");
		#print "@tp\n";
		#my @tp=split(/\t/,$lineIN);
		$_= @_;
		$tpout = $tp[$_[1]];
		for ($i = 2; $i < $_; $i++){
			$tpout = $tpout."\t".$tp[$_[$i]];		
		}
		print OUT $tpout."\n";		
	}		
}


sub GetRowfromNumber{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		#my @tp=split(/\t/,$lineIN);
		$_= @_;
		$tpout = $tp[$_[1]];
		for ($i = 2; $i < $_; $i++){
			$tpout = $tpout."\t".$tp[$_[$i]];		
		}
		print OUT $tpout."\n";		
	}		
}

sub TrimFastaHeadinfo{
	while($lineIN=<IN>){
		chomp($lineIN);
		if (index($lineIN,'>') >= 0){
			$lineIN=substr($lineIN,1,length($lineIN));
			my @tp = &Getheadinfo($lineIN);
			$_= @_;
			$tpout = $tp[$_[1]];
			for ($i = 2; $i < $_; $i++){
				$tpout = $tpout."\t".$tp[$_[$i]];		
			}
			print OUT ">".$tpout."\n";		
		}else{
			print OUT $lineIN."\n";
		}
	}		
}



sub countseqlength{
	my $tplen = 0;
	my $name = "";
	my $first=1;
	while($lineIN=<IN>){
		chomp($lineIN);		
		if (index($lineIN,'>') >= 0){
			my @tp = &Getheadinfo($lineIN);
			if ($first == 1){
				$first=-1;
			}else{
				print OUT "$name\t$tplen\n";
				#print "$name $tplen\n";
				$tplen = 0;
			}
			$name = substr($tp[0],1,length($tp[0]));
		}else{				
			$tplen += length($lineIN);
		}
	}
	print OUT "$name\t$tplen\n";
}

#batch simple compute ^t row ^t Do wat ^t constant 		
sub BatchSimpleCompute{
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp = &Getheadinfo($lineIN);
		my $r;
		if ($_[2] eq 'mult'){
			$r = $tp[$_[1]]* $_[3];
		}
		if ($_[2] eq '+'){
			$r = $tp[$_[1]] + $_[3];
		}
		if ($_[2] eq '-'){
			$r = $tp[$_[1]] - $_[3];
		}
		if ($_[2] eq 'div'){
			$r = $tp[$_[1]] / $_[3];
		}
		print OUT $lineIN."\t".$r."\n";
	}	
}


sub BatchComputeMedian{
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);
		my $tp=@tp;
		my $r;
		my @sorted = sort{$a <=> $b } @tp;
		if($tp%2 ==0){
			$r=($sorted[$tp/2]+$sorted[$tp/2-1])/2;
		}else{
			$r=$sorted[$tp/2-0.5];
		}
		print OUT $lineIN."\t$r\t$sorted[0]\t$sorted[$tp-1]\n";
	}	
}


sub BatchMutiRowCompute{
	#seperate
	#+ 2 3 4 5 6 7 : + 4 5 6
	my $n=@_;
	my @row;  #00-method1 01-no 02-n1 03-n2 
	my $row=0;
	$row[0][0]=$_[1];
	$row[0][1]=2;
	for(my $i=2;$i<$n;$i++){
		if($_[$i] eq ":"){
			$row++;
			$i++;
			$row[$row][0]=$_[$i];
			$row[$row][1]=2;
		}else{  #nurmic
			$row[$row][$row[$row][1]]=$_[$i];
			$row[$row][1]++;
		}
	}
	$row++;
	print "$row\n";
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		my @tp = &Getheadinfo($lineIN);	
		print OUT "$lineIN";	
		for(my $i=0;$i<$row;$i++){
			#print "$row[$i][0] $row[$i][1] $row[$i][2] $row[$i][3] $row[$i][4]\n";
			my $r=$tp[$row[$i][2]];
			for(my $j=3;$j<$row[$i][1];$j++){
				$r=&doSampleCompute($row[$i][0],$r,$tp[$row[$i][$j]]);
				#print "$r $row[$i][0] $tp[$row[$i][$j]]\n";
			}			
			print OUT "\t$r";
		}
		print OUT "\n";		
	}	
}

sub BatchTwoRowCompute{
	#seperate
	my $n=@_;
	my @dowhat;  #0-method 1=n 2-n1 3-n2
	my $d=0;
	my $tpn=2;
	for(my $i=1;$i<$n;$i++){
		if($_[$i] eq "+" || $_[$i] eq "-" || $_[$i] eq "mult" || $_[$i] eq "div"){			
			$dowhat[$d][0]=$_[$i];
			#print "$dowhat[$d][0]\n";
		}elsif($_[$i] eq ":" ){
			$dowhat[$d][1]=$tpn;			
			$d++;
			$tpn=2;
		}else{
			$dowhat[$d][$tpn]=$_[$i];
			#print "$dowhat[$d][$tpn]\n";
			$tpn++;
		}						
	}
	$dowhat[$d][1]=$tpn;
	$d++;
	print "$d pairs\n";
	for(my $i=0;$i<=$d;$i++){
		print "$dowhat[$i][0] $dowhat[$i][1]";
		for(my $j=2;$j<$dowhat[$i][1];$j++){
			print " $dowhat[$i][$j]";
		}
		print "\n";
	}
	
	while($lineIN=<IN>){
		chomp($lineIN);	
		print OUT $lineIN;
		my @tp = &Getheadinfo($lineIN);
		for(my $i=0;$i<$d;$i++){
			my $r=$tp[$dowhat[$i][2]];
			for(my $j=3;$j<$dowhat[$i][1];$j++){
				#print "$dowhat[$i][0],$r,$tp[$dowhat[$i][$j]]\n";
				$r=&doSampleCompute($dowhat[$i][0],$r,$tp[$dowhat[$i][$j]]);
			}
			#print "$r";
			print OUT "\t$r";
		}
		print OUT "\n";
	}	
	print "done\n";
}

sub doSampleCompute{
	#print "$_[0] $_[1] $_[2]\n";
	my $o=0;
	if ($_[0] eq 'mult'){
		$o = $_[1]* $_[2];
	}
	if ($_[0] eq '+'){
		$o = $_[1]+ $_[2];
	}
	if ($_[0] eq '-'){
		$o = $_[1]- $_[2];
	}
	if ($_[0] eq 'div'){
		$o = $_[1]/ $_[2];
	}	
	return $o;
}



sub BreakSeqbyNNN_loc{
	my $tps;

	my $count = 1;
	my $seq;

	my $outfname=$outf.".loc";
	open(OUTF,'>',"$outfname");	

	while($lineIN=<IN>){
		chomp($lineIN);
		if (index($lineIN,'>') >= 0){
			#print OUT
			if (length($seq)>0){
				my $st=1;
				my $ed=1;
				my $seqlen=length($seq);
				my @tp=split/[Nn]/,$seq;
				my $tp=@tp;
				for(my $i=0;$i<$tp;$i++){
					my $len=length($tp[$i]);
					if($len>0){
						print OUT $tps."=".$count."\n";
						&regulatedOUT($tp[$i]);
						$ed=$st+$len;
						print OUTF $tps."=".$count."\t$st\t$ed\t$seqlen\n";
						$count++;
						$st=$ed+1;
					}else{
						$st++;
					}
					
				}
			}
			
			#name
			my @tp=&Getheadinfo($lineIN);
			$tps=$tp[0];
			$count = 1;
			$seq='';
		}else{
			#if NNN new contig
			$seq=$seq.$lineIN;
		}
	}
	#print last one
	my @tp=split/[Nn]/,$seq;
	my $tp=@tp;
	my $st=1;
	my $ed=1;
	my $seqlen=length($seq);
	for(my $i=0;$i<$tp;$i++){
		my $len=length($tp[$i]);
		if($len>0){
			print OUT $tps."=".$count."\n";
			&regulatedOUT($tp[$i]);
			$ed=$st+$len;
			print OUTF $tps."=".$count."\t$st\t$ed\t$seqlen\n";
			$count++;		
			$st=$ed+1;
		}else{
			$st++;
		}			
	}	
	print "\ndone\n";
	close(OUTF);
}



#break sequence by NNN
sub BreakSeqbyNNN{
	my $tps;
	#my $cut=20;
	#if($_[1] >0){
	#	$cut=$_[1];
	#}
	my $count = 1;
	my $seq;
	while($lineIN=<IN>){
		chomp($lineIN);
		if (index($lineIN,'>') >= 0){
			#print OUT
			if (length($seq)>0){
				my @tp=split/[Nn]/,$seq;
				my $tp=@tp;
				for(my $i=0;$i<$tp;$i++){
					if(length($tp[$i])>0){
						print OUT $tps."=".$count."\n";
						&regulatedOUT($tp[$i]);
						$count++;		
					}					
				}
			}
			
			#name
			my @tp=&Getheadinfo($lineIN);
			$tps=$tp[0];
			$count = 1;
			$seq='';
		}else{
			#if NNN new contig
			$seq=$seq.$lineIN;
		}
	}
	#print last one
	my @tp=split/[Nn]/,$seq;
	my $tp=@tp;
	for(my $i=0;$i<$tp;$i++){
		if(length($tp[$i])>0){
			print OUT $tps."=".$count."\n";
			&regulatedOUT($tp[$i]);
			$count++;		
		}					
	}	
	print "\ndone\n";
}


sub GetPairLTRbyDirect
{
	#in
	#RLOth_654	Chr1	90.46	262	25	0	1	262	18007222	18007483	1.00E-86	321	210	+
	#list distancd for select
	#out
	#Family	Chrom	direct	s1	e1	s2	e2	distance
	my $d=13;
	if ($_[1]>0){
		$d=$_[1];
	}
	my @info;
	my $out=0;
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);
		if ($info[0] eq $tp[0] && $info[1] eq $tp[1] && $info[2] eq $tp[$d]){
			my @tpinfo=($tp[8],$tp[9]);
			@info=(@info,@tpinfo);
			$out++;
			print ".";
		}else{
			if($out==1){
				my $distence=$info[5]-$info[4];
				#print "$distence ";
				print OUT "@info\t$distence\n";
			}
			elsif($out>1){
				print OUT "@info\n";
			}
			@info=();
			$out=0;
			@info = ($tp[0],$tp[1],$tp[$d],$tp[8],$tp[9]);
		}		
	}
	print "\n";
}


sub GetFirstSubtype{
	my $type="";
	my $subtype="";
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		if ($type eq $tp[$_[1]]){
			if($subtype eq $tp[$_[2]]){
				print OUT $lineIN."\n";
			}
		}else{
			$type=$tp[$_[1]];
			$subtype=$tp[$_[2]];
			print OUT $lineIN."\n";
		}		
	}
}


sub GetTypeLongestfas{
	my $add=1;
	my $name="NA";
	my @seq;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN,"=","|");
		if($lineIN=~/>/){			
			if($tp[0] eq $name){
				$count++;	
				$seq[$count]="";			
			}else{
				if($name ne "NA"){
					#output
					my $max=0;
					my $maxi=0;
					for(my $i=0;$i<$count;$i++){
						my $len=length($seq[$count]);
						if($max>$len){
							$maxi=$i;
						}
					}
					print OUT "$seq[$maxi]\n";
				}
				$name =$tp[0];
				print OUT $name."\n"; 
				$count=0;
				@seq=();				
			}	
					
		}else{
			$seq[$count]=$seq[$count].$lineIN;
		}
	}
	
	#last one
	my $max=0;
	my $maxi=0;
	for(my $i=0;$i<$count;$i++){
		my $len=length($seq[$count]);
		if($max>$len){
			$maxi=$i;
		}
	}
	print OUT "$seq[$maxi]\n";	
}

sub Gettypefirstsequence{
	my $add=1;
	my $name=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		if($lineIN=~/>/){
			if($lineIN eq $name){
				$add=0;
				print ".";
			}else{
				$add=1;
			}
			$name = $lineIN;
		}
		if($add==1){
			print OUT $lineIN."\n";
		}
	}
}


sub GettypeToprows{
	my $tpinfo="NA";
	my $loc=$_[1];
	my $top=$_[2];
	my $count=0;
	print "Get Local $loc and Top $top\n";
	#
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $lineinfo=$tp[$loc[0]];	
		if ($tpinfo ne $lineinfo){
			print OUT $lineIN."\n";
			$tpinfo = $lineinfo;
			$count=0;
			#print "+";
		}else{
			$count++;
			if($count < $top){
				print OUT $lineIN."\n";
			}else{
				#print "-";
			}
		}
	}
}


sub Gettypefirstrow{
	my $tpinfo="NA";
	my @loc={0};
	my $infono=@_;
	$infono--;
	if($infono>=1){
		for(my $i=0;$i<$infono;$i++){
			$loc[$i]=$_[$i+1];
			#print "$_[$i+1]\n";
		}
	}
	#print "$infono $loc[0] $loc[1] $_[0] $_[1]\n";
	#
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $lineinfo=$tp[$loc[0]];
		for(my $i=1;$i<$infono;$i++){
			$lineinfo=$lineinfo.$tp[$loc[$i]];
		}		
		if ($tpinfo ne $lineinfo){
			print OUT $lineIN."\n";
			$tpinfo = $lineinfo;
			print "+";
		}else{
			print "-";
		}
	}
}

sub regulatestsfile{	
	while($lineIN=<IN>){
		chomp($lineIN);
		if (index($lineIN,'>') >= 0){
			if($tpinfo ne $lineIN){
				print OUT "\n".$lineIN;
				$tpinfo = $lineIN;
			}
		}else{
			print OUT "\t".$lineIN;
		}
	}
}	

sub addsubtype{	
	while($lineIN=<IN>){
		chomp($lineIN);
		my $l=@_;
		for (my $i=1; $i<$l; $i++){
			print OUT $lineIN."\t".$_[$i]."\n";
		}
	}
}

sub deletenestloc{
	@tploc;  # 0-loc, 1-start, 2-end, 3,4,5, etc
	$tpcount=0;
	my $tptype;
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);
		print $tp[0]." ".$tp[1]."\n";
		my @tpone=($tp[$_[2]],$tp[$_[3]],$tp[$_[4]]);
		if ($tp[$_[1]] ne $tptype){  #type
			print OUT $lineIN;	
			$tpcount=1;
			@tploc=@tpone;
			$tptype=$tp[$_[1]];
		}else{
			my $test=-1;
			for(my $i=0; $i<$tpcount; $i++){
				@searchone=($tploc[$i*3], $tploc[$i*3+1], $tploc[$i*3+2]);
				@test=Comparelocal(@tpone, @searchone);
				if ($test[0]>0){
					last;
				}
			}
			print $test."\n";
			if ($test<0){   #add new
				@tploc[$tpcount*3]=$searchone[0];	
				@tploc[$tpcount*3+1]=$searchone[1];	
				@tploc[$tpcount*3+2]=$searchone[2];
				$tploc=@tploc;
				print "$tploc\n";
				print OUT $lineIN;	
				$tpcount++;
			}
		}
	}
}
	
	
sub Floot_2_int_up{	
	my $n=0;
	if($_[1] >0 ){
		$n=$_[1];
	}
	
	#
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		#my @tp = split(/\t/,$lineDB);
		my @f = split(/\./,$tp[$n]);			
		if($f[1]>0){
			$tp[$n]=$f[0]+1;
		}else{
			$tp[$n]=$f[0];
		}		
		print OUT "@tp\n";
	}

}	

	
sub FindSameModule{	
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);
		my $tp=@tp;
		my $m0=$tp[0];
		my @samemodule;
		my $count=0;
		for (my $i=1;$i<$tp;$i+=$_[1]){
			#compair
			my $tt=&Modulenatch($tp[0], $tp[$i]);
			if ($tt > 0){
				$samemodule[$count]=$i;
				$count++;
			}
		}
		#find max
		my $nmax=0;
		my $moduleok;
		for (my $i=0;$i<$count;$i++){
			if ($tp[$samemodule[$i]+1]>$nmax){
				$nmax=$tp[$samemodule[$i]+1];
				$moduleok=$samemodule[$i];
			}
		}
		#result
		my $r=$tp[0];
		for (my $i=0; $i<$_[1];$i++){
			$r=$r."\t".$tp[$moduleok+$i];		
		}
		print OUT $r."\n";		
	}
}	

sub countstringlength{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $len=length($tp[$_[1]]);
		print OUT $lineIN."\t".$len."\n";
	}
}	

sub Typelowvaluemask{
	my $tptype="0";
	my $tpmaxno;
	my @valuelist;
	my	$count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		if ($tptype eq "0"){
			$tptype=$tp[0];
		}
		if ($tptype ne $tp[0]){
			#find max
			for (my $i=1; $i<$count; $i++){
				if ($valuelist[$tpmaxno] < $valuelist[$i]){
					$tpmaxno=$i;
				}
			}
			#output
			for (my $i=0; $i<$count; $i++){
				my $r=$tptype."\t".$valuelist[$i];
				if ($tpmaxno!=$i){
					$r=$r."\t"."-";
				}
				print OUT $r."\n";
			}
			#new data
			$tptype=$tp[0];
			$count=0;
			$tpmaxno=0;			
		}
		$valuelist[$count]=$tp[1];
		$count++;		
	}
	#find max
	for (my $i=1; $i<$count; $i++){
		if ($valuelist[$tpmaxno] < $valuelist[$i]){
			$tpmaxno=$i;
		}
	}
	#output
	for (my $i=0; $i<$count; $i++){
		my $r=$tptype."\t".$valuelist[$i];
		if ($tpmaxno!=$i){
			$r=$r."\t"."-";
		}
		print OUT $r."\n";
	}
}	

sub Combinclusterrecord{
	#ID as hash key,value as sorted integer
	#each pair, check record,
	#if new, add a new figure
	#if have, add same value
	#if both have,change v2 to v1,search other hash,change all v2 to v1
	$count=10;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $v0 = $db{$tp[0]};
		my $v1 = $db{$tp[1]};
		if ($v0<=0){
			if ($v1<=0){ #add
				$db{$tp[0]}=$count;	
				$db{$tp[1]}=$count;
				#print $tp[0]." ".$tp[1]." ".$db{$tp[0]}."\n";
				$count++;
			}else{	#$v0=$v1
				$db{$tp[0]}=$v1;
			}			
		}else{
			if ($v1<=0){ #$v1=$v0
				$db{$tp[1]}=$v0;	
			}else{	#all $v1=$v0
				#print "combin ".$v1."\n";
				while ((my $tpkey, my $tpvalue) = each(%db)) {
					if ($tpvalue==$v1){
						$db{$tpkey}=$v0;
					}
				}				
			}
		}
	}
	#regulate to other hash
	my %o;
	while ((my $dbkey, my $dbvalue) = each(%db)) {
		#print $dbkey." ".$dbvalue."\n";
		my $outvalue=$o{$dbvalue};
		if ($outvalue ne ""){
			$o{$dbvalue}=$outvalue."\t".$dbkey;
		}else{
			$o{$dbvalue}=$dbkey;
		}
	}
	
	#output
	my $outcount=1;
	while ((my $outkey, my $outvalue) = each(%o)) {	
		print OUT $outcount."\t".$outvalue."\n";
		$outcount++;
	}		
}

sub Listallelement{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $tp=@tp;
		for (my $i=0;$i<$tp;$i++){
			print OUT $tp[$i]."\n";
		}
	}
}


sub ListLineElementNo{
	
	my $count=0;
	my @element;

	print OUT "Type\t";
	$_=@_;
	for(my $i=1;$i<$_;$i++){
		$element[$count]=$_[$i];
		$count++;
		print OUT "$_[$i]\t";
	}
	print OUT "Total\n";	
	print "$count elements\n";	
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $tp=@tp;
		my @eno;
		my $total=0;
		for(my $i=0;$i<$count;$i++){
			$eno[$i]=0;
		}
		#check	
		for (my $i=0;$i<$tp;$i++){
			for(my $j=0;$j<$count;$j++){
				if($tp[$i] eq $element[$j]){
					$eno[$j]++;
					$total++;
					last;
				}
			}
		}
		#out
		print OUT "$tp[0]";
		for(my $i=0;$i<$count;$i++){
			print OUT "\t$eno[$i]";
		}
		print OUT "\t$total\n";
	}
}


sub Addlisttoinfile{
	my $count=$_[1];
	my $n=$_[2];
	while($lineIN=<IN>){
		chomp($lineIN);
		my $n=$_[2] - length($count);
     	my $tpadd=$count;
     	for(my $j=0;$j<$n;$j++){
     		$tpadd="0".$tpadd;
     	}
     	my $add=$_[3].$tpadd.$_[4];
		print OUT $add."\t".$lineIN."\n";	
		$count++;	
	}
}


sub AddlistandReplacetoDuplicateID{
	my $add=1;
	my @readin;
	my @date;	#0-check ID, 1,2,3,all
	my $count=0;
	my $n=0;
	my $loc=0;
	
	#read 
	while($lineIN=<IN>){
		chomp($lineIN);
		$readin[$count]=$lineIN;
		$count++;	
	}
	
	#sort
	@readin=sort(@readin);
	
	#print "@readin\n";
	
	#change
	
	for(my $i=0;$i<$count;$i++){
		my @tp = &Getheadinfo($readin[$i]);
		$n=@tp;
		$date[$i]=[@tp];
		$date[$i][$n]=$tp[$loc];
		#print "$date[$i][0]\n";
		if($i>0){
			if ($date[$i][$n] eq $date[$i-1][$n]){		
				print ".";		
				if($date[$i][$loc] eq $date[$i-1][$loc]){
					#chang last
					$date[$i-1][$loc]=$date[$i-1][$loc].".".$add;
					$add++;
				}
				#change current
				$date[$i][$loc]=$date[$i][$loc].".".$add;
				$add++;
							
			}else{
				#reset $add
				$add=1;
			}	
		}
	}	
	#print "@date\n";
	#output
	for(my $i=0;$i<$count;$i++){
		print OUT "$date[$i][0]";
		for(my $j=1;$j<$n;$j++){
			print OUT "\t$date[$i][$j]";
		}
		print OUT "\n";
	}
}


sub Addlisttoeachtype{
	my $count=$_[2];
	my $tptype='start';
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		if ($tp[$_[1]] eq $tptype){
			print OUT $lineIN."\t".$count."\n";
			$count++;			
		}else{
			$count=$_[2];
			$tptype=$tp[$_[1]];
			print OUT $lineIN."\t".$count."\n";
			$count++;
		}			
	}
}

sub CombinFastaSeqtoOne{
	my $add=0;
	my $no=200;
	my $sy='N';
	if ($_[1]>0){
		$no=$_[1];
	}
	if ($_[2] ne ''){
		$sy=$_[2];
	}
	while($lineIN=<IN>){
		if($add==0){
			print OUT ">Chr1\n";
			$add++;
		}else{
			chomp($lineIN);
			if($lineIN=~/>/){
				for(my $i=0; $i<$no;$i++){
					print OUT "$sy";
				} 	
				print OUT "\n";
			}else{
				print OUT "$lineIN\n";
			}
		}
	}
}		
		
sub Combinmutirow{
	my @c=();
	my $cno=@_;

	for(my $i=1;$i<$cno;$i++){
		$c[$i-1]=$_[$i];
	}
	$cno--;
	
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $tp=@tp;
		my $combin=$tp[$c[0]];
		for(my $i=1;$i<$cno;$i++){
			$combin=$combin."=".$tp[$c[$i]];
		}
		print OUT "$combin\t$lineIN\n";
	}
}	

sub Findnestloc{
	#input
	#name chr st ed direct
	#CUFF.31799.1	Chr9	22940001	22941221	- 

	my @ed;
	my @name;
	my $count=0;
	my $chrom;
	my $max=100000;
	if($_[1]>0){
		$max=$_[1];
	}
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		if($tp[2]-$ed[$count]<$_[1] && $tp[1] eq $chrom){
			for (my $i=0;$i<=$count;$i++){
				my $rage=$tp[2]-$ed[$count-$i];
				if($rage<$_[1]){
					my $o=$_[1]-$rage;
					print OUT $tp[0]."\t".$name[$count-$i]."\t".$rage."\t".$o."\n";
					print OUT $name[$count-$i]."\t".$tp[0]."\t".$rage."\t".$o."\n";
				}else{
					last;
				}				
			}
			$count++;									
		}else{
			$count=0;
		}

		$ed[$count]=$tp[3];
		$name[$count]=$tp[0];
		$chrom=$tp[1]
	}	
	#check
	
}

sub Deleteshortline{	
	
	####???????????????? 
	my @ed;
	my @name;
	my $count=0;
	my $chrom;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &_($lineIN);
		if($tp[2]-$ed[$count]<$Findnestloc[1] && $tp[1] eq $chrom){
			for (my $i=0;$i<=$count;$i++){
				my $rage=$tp[2]-$ed[$count-$i];
				if($rage<$Findnestloc[1]){
					my $o=$Findnestloc[1]-$rage;
					print OUT $tp[0]."\t".$name[$count-$i]."\t".$rage."\t".$o."\n";
					print OUT $name[$count-$i]."\t".$tp[0]."\t".$rage."\t".$o."\n";
				}else{
					last;
				}				
			}
			$count++;									
		}else{
			$count=0;
		}

		$ed[$count]=$tp[3];
		$name[$count]=$tp[0];
		$chrom=$tp[1]
	}	
}	

#0-name 1-chrom 
#record combin need same derict
#querry two record can't overlap
#exceed querry distance 2 times > target distance
sub CombinlBlastResult{	
	my @info;
	my @count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
		my $tp=@tp;
		if ($info[0][0]==$tp[0] && $info[0][1]==$tp[1]){  
			#add info and test combin
			my $canadd=1;
			for (my $i=0; $i<$count;$i++){
				#estimate querry
				if($tp[6]>$info[$i][6]){
					my $d=$tp[6]-$info[$i][7];
					if($d>0){
						#estimate target
						if($tp[8]<$tp[9] && $info[$i][8]<$info[$i][9]){
							my $t=$tp[8]-$info[$i][9];
							if($t>0 && $t<$d*4){
								#combin
								$info[$i][7]=$tp[7];
								$info[$i][9]=$tp[9];
								$info[$i][12]='=1';
								$info[$i][3]=$info[$i][7]-$info[$i][6]+1;
								$canadd=-1;
								last;
							}
						}
						if($tp[8]>$tp[9] && $info[$i][8]>$info[$i][9]){
							my $t=$info[$i][9]-$tp[8];
							if($t>0 && $t<$d*4){
								#combin
								$info[$i][7]=$tp[7];
								$info[$i][9]=$tp[9];
								$info[$i][12]='=2';
								$info[$i][3]=$info[$i][7]-$info[$i][6]+1;
								$canadd=-1;
								last;
							}						
						}
					}
				}else{
					my $d=$info[$i][6]-$tp[7];
					if($d>0){
						#estimate target
						if($tp[8]<$tp[9] && $info[$i][8]<$info[$i][9]){
							my $t=$info[$i][8]-$tp[9];
							if($t>0 && $t<$d*4){
								#combin
								$info[$i][6]=$tp[6];
								$info[$i][9]=$tp[9];
								$info[$i][12]='=3';
								$info[$i][3]=$info[$i][7]-$info[$i][6]+1;
								$canadd=-1;
								last;
							}
						}
						if($tp[8]>$tp[9] && $info[$i][8]>$info[$i][9]){
							my $t=$tp[9]-$info[$i][8];
							if($t>0 && $t<$d*4){
								#combin
								$info[$i][6]=$tp[6];
								$info[$i][9]=$tp[9];
								$info[$i][12]='=4';
								$info[$i][3]=$info[$i][7]-$info[$i][6]+1;
								$canadd=-1;
								last;
							}						
						}						
					}				
				}
			}			
			if ($canadd==1){
				$info[$count]=[@tp];
				$count++;
			}
		}else{
			if ($info[0][0]){
				#output
				for (my $i=0; $i<$count;$i++){
					my $r=$info[$i][0];
					for (my $j=1; $j<13; $j++){
						$r=$r."\t".$info[$i][$j];
					}	
					print OUT $r."\n";
				}
				$count=0;
			}
			$info[$count]=[@tp];
			$count++;	
			
		}	
	}
	
	for (my $i=0; $i<$count;$i++){
		my $r=$info[$i][0];
		for (my $j=1; $j<13; $j++){
			$r=$r."\t".$info[$i][$j];
		}	
		print OUT $r."\n";
	}	
}

#If some records was in same location, keep the record With have largest e value
#blast result have 12 rank,e.g.
#pa007_2_G	supercontig_69	78.69	413	80	4	2	406	605484	605072	6e-019	97.6
sub CombinltoBestQuerry{
	#use two-dimensional array
	#$b[$count]=[@tp];
	my @b;   #blast result
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $canadd=0;
		#check cross and reiteration
		for (my $i=0;$i<$count;$i++){
			my @a=&Comparelocal($tp[1],$tp[8],$tp[9],$b[$i][1],$b[$i][8],$b[$i][9]);
			if ($a[0]>=1){
				#reiteration,select by e value 
				if ($tp[10]<=$b[$i][10]){
					$b[$i]=[@tp];
				}else{
				}
				$canadd=1;
				last;	
				
			}
		}
		#if no reiteration, add record
		if ($canadd==0){
			$b[$count]=[@tp];
			$count++;
		}
	}	
	#output
	for (my $i=0;$i<$count;$i++){
		my $o=$b[$i][0];
		for (my $j=1;$j<12;$j++){
			$o=$o."\t".$b[$i][$j];
		}
		print OUT $o."\n";
	}
}

sub SortbyrownumberRepeat{
	my @r;
	my $count=0;
	while($lineIN=<IN>){
		chomp($lineIN);	
		$r[$count]=$lineIN;
		$count++;
	}
		my @sortr = sort{&memberlen($a) <=> &memberlen($b) } @r;
		
		#@array=sort {$a <=> $b} @array;
	for (my $i=$count -1; $i>=0;$i--){
		#print "$r[$i]\n";
		print OUT "$sortr[$i]\n";
	} 
}
sub memberlen{
	my @t=&Getheadinfo($_[0]);
	$t=@t;
	return $t;
}
		
sub Sortbyrownumber{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $k=$tp[$_[1]];
		#print $k."\n";
		#muti
		if ($db{$k}){
			if ($dba{$k}){
				$dba{$k}++;
				my $c=$k.$dba{$k};
				$db{$c}=$lineIN;
			}else{
				$dba{$k}=1;
				my $c=$k.$dba{$k};
				$db{$c}=$lineIN;
			}
			
		}else{
			$db{$k}=$lineIN;
		}		
	}
	
	my @dbkey = keys(%db);
	my @dbkey2=sort(@dbkey);
	my $count=@dbkey2;
	print $count."\n";
	if ($_[2] eq "+"){
		#A-Z	
		for(my $i=0;$i<$count;$i++){
			my $o=$db{@dbkey2[$i]};
			print OUT $o."\n";
		}
	}else{
		#Z-A
		for(my $i=$count-1;$i>=0;$i--){
			my $o=$db{@dbkey2[$i]};
			print OUT $o."\n";
		}		
	}
}	

sub Deletesameline{
	my $l='1';
	while($lineIN=<IN>){
		chomp($lineIN);
		if ($l eq "1"){
			$l=$lineIN;
			print OUT $lineIN."\n";
		}else{
			if ($l ne $lineIN){
				print OUT $lineIN."\n";
				$l=$lineIN;
			}
		}
	}
}	

sub FindPEmatch{
	#read in same name f and r, save to two-dimensional array,meat nest BES pire, analysis last group
	#analysis over, output and read in next group
	my @pe;
	my $name='first';
	my $count=0;
	my $infocount=$_[2];
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp=&Getheadinfo($lineIN);
		my $l=length($tp[9]);	
		my $tpname=substr($tp[9],0,$l-1);
		if($name ne $tpname){
			#analysis
			if($name ne 'first'){
				for(my $i=0; $i<$count;$i++){
					#use first direct to search second direct
					if($pe[0][9] ne $pe[$i][9]){
						last;
					}
					my $tpcount=1;
					for(my $j=$i+1;$j<$count;$j++){
						if($pe[$j][9] ne $pe[$i][9]){
							#compare locate
							my $d=$pe[$i][15]-$pe[$j][15];
							if($d<$_[1] && $d>(0-$_[1])){
								#output
								my $tpout=$d."\t".$tpcount;
								for (my $k=0;$k<=$infocount;$k++){
									$tpout=$tpout."\t".$pe[$i][$k];
								}
								for (my $k=0;$k<=$infocount;$k++){
									$tpout=$tpout."\t".$pe[$j][$k];
								}
								print OUT $tpout."\n";
								$tpcount++;
							}							
						}	
					}
				}	
			}
			#clean
			$count=0;
			$name=$tpname;
		}
		#add record
		$pe[$count]=[@tp];
		$count++;
	}
}	

sub add0behandnumber{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $tp=@tp;
		my $la=$_[2]-length($tp[$_[1]]);
		for (my $j=0;$j<$tp;$j++){
			if ($j==$_[1]){
				for (my $i=0;$i<$la;$i++){
					print OUT '0';
				}
			}
			print OUT $tp[$j];
			if ($j != $tp - 1){
				print OUT "\t";
			}else{
				print OUT "\n";
			}		
		}
	}
}

sub remove0behandnumber{
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		my $tp=@tp;
		#change
		$tp[$_[1]]+=0;

		#output
		for (my $j=0;$j<$tp;$j++){
			print OUT $tp[$j];
			if ($j != $tp - 1){
				print OUT "\t";
			}else{
				print OUT "\n";
			}		
		}
	}
}	
	
sub Countsubtypeno{
	my $count=1;
	my @subtype;
	my $type='start';
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		if ($type eq $tp[$_[1]]){
			my $canadd=0;
			for (my $i=0;$i<$count;$i++){
				if ($tp[$_[2]] eq $subtype[$i]){
					$canadd=1;
					print 
				}
			}
			if ($canadd==0){
				$subtype[$count]=$tp[$_[2]];
				$count++;
			}
		}else{
			if ( $type ne 'start'){
				print OUT $type."\t".$count."\n";
			}
			$type = $tp[$_[1]];
			$count=1;
			$subtype[0]=$tp[$_[2]];
		}			
	}
	print OUT $type."\t".$count."\n";					
}				

sub Countsubloc{
	#range=5 means in 5 lines same name combin to same loc
	my @subtype;
	my @range;
	my $count=0;
	my %type;
	while($lineIN=<IN>){
		my @tp = &Getheadinfo($lineIN);	
		my $canadd=0;
		for (my $i=1;$i<=$_[2];$i++){
			if ($tp[$_[1]] eq $range[$count-$i]){
				$canadd=1;
				last;
			}
		}
		if ($canadd==0){
			print "$tp[$_[1]]\n";
			if ($type{$tp[$_[1]]}){
				$type{$tp[$_[1]]}++;
			}else{
				$type{$tp[$_[1]]}=1;
			}
		}
		$range[$count]=$tp[$_[1]];
		$count++;	
	}
	my @typelist=keys(%type);
	my $typelist=@typelist;
	for(my $i=0;$i<$typelist;$i++){
		print OUT $typelist[$i]."\t".$type{$typelist[$i]}."\n";
	}			
}

sub Combinloc{
	#name start end 
	my $name='start';
	my $start;
	my $endl;
	while($lineIN=<IN>){
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);	
		if ($name eq $tp[0]){
			my $distance=$tp[1]-$endl;
			if($distance<=$_[1]){
				$endl=$tp[2];
			}else{
				print OUT $name."\t".$start."\t".$endl."\n";	
				$name=$tp[0];
				$start=$tp[1];
				$endl=$tp[2];	
			}
		}else{
			if($name ne 'start'){
				print OUT $name."\t".$start."\t".$endl."\n";
			}
			$name=$tp[0];
			$start=$tp[1];
			$endl=$tp[2];
		}	
	}
	print OUT $name."\t".$start."\t".$endl."\n";	
}


sub sortSequencebyLength{
	print "sortSequencebyLength\n";
	my $name;
	my %seq;
	my $seqtp="NA";
	my $add=1000000;
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		if($lineIN=~/>/){	
			if($seqtp ne "NA"){					
				my $len=$add+length($seqtp);				
				$tpname=$len."\t$name";
				$seq{$tpname}=$seqtp;
				#print "$tpname\n";
			}			
			$seqtp="";
			$name=$lineIN;
		}else{								
			$seqtp=$seqtp.$lineIN;			
		}
			
	}
	
	#sort
	my @o=keys(%seq);
	@o=sort(@o);
	my $o=@o;
	print "$o seq readed\n";
	print $seq{"1000362 >w10_157_5"}."\ntest\n";
	#printout
	for (my $i=0;$i<$o;$i++){
		my @tp = &Getheadinfo($o[$i]);	
		my $outseq=$seq{$o[$i]};
		#print "$o[$i] ";
		print OUT "$tp[1]\n$outseq\n";
	}
}

sub sortSequencebyName{
	print "sortSequencebyName\n";
	my @name;
	my %seq;
	my $count=-1;
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		if($lineIN=~/>/){
			$count++;
			$name[$count]=$lineIN;
		}else{
			$seq{$name[$count]}=$seq{$name[$count]}.$lineIN;
		}	
	}
	#sort
	@name=sort(@name);
	#printout
	for (my $i=0;$i<=$count;$i++){
		my $outseq=$seq{$name[$i]};
		print OUT "$name[$i]\n$outseq\n";
	}
}


sub sortbyChr{
	print "sortbyChr\n";
	my %row;  
	my $count=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		#my @tp = &Getheadinfo($lineIN);
		my @tp = split(/\t/,$lineIN);			
		if($row{$tp[$_[1]]}){
			$row{$tp[$_[1]]}=$row{$tp[$_1]}."$lineIN\n";
		}else{
			$row{$tp[$_[1]]}="$lineIN\n";
			$count++;	
			print "$tp[$_[1]]\n";
		}
		
	}
	print "read $count chr\n";
	
	my @ckeys=keys(%row);

	#sort
	@ckeys=sort(@ckeys);
	#printout
	for (my $i=0;$i<$count;$i++){
		print OUT $row{$ckeys[$i]};
	}

	print "\n";
}

sub sortbyStringrow{
	print "sortbyStringrow\n";
	my @row;
	my $count=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		my @tp = &Getheadinfo($lineIN);
		$row[$count]=$tp[$_[1]].' '.$lineIN ;
		$count++;	
	}
	#sort
	@row=sort(@row);
	#printout
	print "1.";
	if ($_[2] eq "+"){
		for (my $i=0;$i<$count;$i++){
			my $st=index($row[$i], ' ');
			my $out=substr($row[$i],$st+1,length($row[$i]));
			print OUT $out."\n";
		}
		print "2.";
	}else{
		for (my $i=$count-1;$i>=0;$i--){
			my $st=index($row[$i], ' ');
			my $out=substr($row[$i],$st+1,length($row[$i]));
			print OUT $out."\n";
		}	
		print "3.";
	}
	print "\n";
}

sub sortbyStrandIntrow{
	print "sortbyStrandIntrow $_[1] $_[2]\n";
	my @b;
	my @a1; #int
	my @a2; #str
	my $count=0;
	my $len=0;
	my $max1=0;  
	my $max2=0;
	my $loc1=0;	#chr
	my $loc2=0;	#int
	if ($_[1]>=0){
		$loc1=$_[1];
	}
	if ($_[2]>=0){
		$loc2=$_[2];
	}
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		my @tp = &Getheadinfo($lineIN);
		if($count==0){
			$len=@tp;
		}
		$a1[$count]=$tp[$loc1];
		$b[$count]=$lineIN;
		$a2[$count]=$tp[$loc2];
		if($max2<length($tp[$loc2])){
			$max2=length($tp[$loc2]);
		}
		if($max1<length($tp[$loc1])){
			$max1=length($tp[$loc1]);
		}
		$count++;	
	}
	#add key
	for(my $i=0;$i<$count;$i++){
		my $add='';
		for(my $j=length($a1[$i]);$j<$max1;$j++){
			$add=$add.' ';	
		}
		for(my $j=length($a2[$i]);$j<$max2;$j++){
			$add=$add.'0';	
		}		
		$b[$i]=$a1[$i].$add.$a2[$i]."\t".$b[$i];
		#print $b[$i]."\n";
	}
	#sort
	print "$count line $len row read\n";
	@b=sort(@b);
	print "Sort OK\nOutput...\n";
	#printout
	if ($_[3] eq "-"){
		for (my $i=$count-1;$i>=0;$i--){
			my $st=index($b[$i], "\t")+1;
			my $o=substr($b[$i],$st,length($b[$i]));
			print OUT $o."\n";			
		}	
		print "-.";
	}else{
		for (my $i=0;$i<$count;$i++){
			my $st=index($b[$i], "\t")+1;
			my $o=substr($b[$i],$st,length($b[$i]));
			print OUT $o."\n";	
		}
		print "+.";
	}
	print "\n";
}

sub sortbyIntOrStrrow{
	#read in sort list
	#int must be short than 10000000000000000000
	$_=@_;
	my $max=100000000000;
	my @list;	#0-type 1-value
	my $c=0;
	for(my $i=1;$i<$_;$i+=2){
		$list[$c][0]=$_[$i];
		$list[$c][1]=$_[$i+1];
		$c++;
	}
	
	#read in data
	my @a;
	my $count=0;
	my $len=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t$lineIN\n";
		my @tp = &Getheadinfo($lineIN);
		my $add="";
		my $l=0;
		for(my $i=0;$i<$c;$i++){
			my $addtp=$tp[$list[$i][1]];					
			if($list[$i][0] eq "int"){
				$addtp=$max+$addtp;
			}
			$add=$add.$addtp."_";
		}
		$a[$count]=$add."\t".$lineIN;
		#print "$a[$count]\n";
		$count++;	
	}

	#sort
	print "$count line $len row read\n";
	@a=sort(@a);
	print "Sort OK\nOutput...\n";
	for (my $i=0;$i<$count;$i++){
		my @tp = &Getheadinfo($a[$i]);
		my $tp=@tp;
		print OUT $tp[1];
		for(my $j=2;$j<$tp;$j++){
			print OUT "\t$tp[$j]";	
		}
		print OUT "\n";
	}
	print "\n";
}



sub sortbyIntandStrrow{
	print "sortbyIntandStrrow\n";
	my @b;
	my @a;
	my $count=0;
	my $len=0;
	my $max=0;
	my $loc1=0;
	my $loc2=1;
	if ($_[1]>0){
		$loc1=$_[1];
	}
	if ($_[2]>0){
		$loc2=$_[2];
	}
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		my @tp = &Getheadinfo($lineIN);
		if($count==0){
			$len=@tp;
		}
		$a[$count]=$tp[$loc1];
		$b[$count]=$tp[$loc2]." ".$lineIN;
		if($max<$tp[$loc1]){
			$max=$tp[$loc1];
		}
		$count++;	
	}
	#add key
	my $maxlen=length($max);
	for(my $i=0;$i<$count;$i++){
		my $add='';
		for(my $j=length($a[$i]);$j<$maxlen;$j++){
			$add=$add.'0';	
		}
		$b[$i]=$add.$a[$i].$b[$i];
		#print $b[$i]."\n";
	}
	#sort
	print "$count line $len row read\n";
	@b=sort(@b);
	print "Sort OK\nOutput...\n";
	#printout
	if ($_[3] eq "-"){
		for (my $i=$count-1;$i>=0;$i--){
			my $st=index($b[$i], ' ')+1;
			my $o=substr($b[$i],$st,length($b[$i]));
			print OUT $o."\n";			
		}	
		print "-.";
	}else{
		for (my $i=0;$i<$count;$i++){
			my $st=index($b[$i], ' ')+1;
			my $o=substr($b[$i],$st,length($b[$i]));
			print OUT $o."\n";	
		}
		print "+.";
	}
	print "\n";
}

sub sort2IntinEachRow{
	print "sort2IntinEachRow\n";
	my $a1=$_[1];
	my $a2=$_[2];

	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		my @tp = &Getheadinfo($lineIN);
		my $tp=@tp;
		my $i1=$tp[$_[1]];
		my $i2=$tp[$_[2]];
		
		if($i1>$i2){
			$tp[$_[1]]=$i2;
			$tp[$_[2]]=$i1;
		}
		#out
		print OUT $tp[0];
		for(my $i=1;$i<$tp;$i++){
			print OUT "\t$tp[$i]";
		}
		print OUT "\n";
	}

	print "\n";
}

sub sortbyIntrow{
	print "sortbyIntrow\n";
	my @b;
	my @a;
	my $count=0;
	my $len=0;
	my $max=0;
	my $loc=0;
	if ($_[1]>0){
		$loc=$_[1];
	}
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		my @tp = &Getheadinfo($lineIN);
		if($count==0){
			$len=@tp;
		}
		$a[$count]=$tp[$loc];
		$b[$count]=$lineIN;
		if($max<$tp[$loc]){
			$max=$tp[$loc];
		}
		$count++;	
	}
	#add key
	my $maxlen=length($max);
	for(my $i=0;$i<$count;$i++){
		my $add='';
		for(my $j=length($a[$i]);$j<$maxlen;$j++){
			$add=$add.'0';	
		}
		$b[$i]=$add.$a[$i]." ".$b[$i];
		#print $b[$i]."\n";
	}
	#sort
	print "$count line $len row read\n";
	@b=sort(@b);
	print "Sort OK\nOutput...\n";
	#printout
	$maxlen++;
	if ($_[2] eq "-"){
		for (my $i=$count-1;$i>=0;$i--){
			my $o=substr($b[$i],$maxlen,length($b[$i]));
			print OUT $o."\n";			
		}	
		print "-.";
	}else{
		for (my $i=0;$i<$count;$i++){
			my $o=substr($b[$i],$maxlen,length($b[$i]));
			print OUT $o."\n";	
		}
		print "+.";
	}
	print "\n";
}

	
sub sortby2Introw{
	print "sortby2Introw\n";
	#if indrect, change it
	my @b;
	my @row;
	my $count=0;
	while($lineIN=<IN>){	
		chomp($lineIN);
		#print "$count\t";
		my @tp = &Getheadinfo($lineIN);
		if($tp[$_[1]]>$tp[$_[2]]){
			my $t=$tp[$_[1]];
			$tp[$_[1]]=$tp[$_[2]];
			$tp[$_[2]]=$t;
		}
		$b[$count]=[@tp];
		$row[$count]=$tp[$_[1]];
		$count++;	
	}
	#sort
	my @sortrow=sort{$a <=> $b}(@row);
	#printout
	print "1.";
	if ($_[2] eq "+"){
		for (my $i=0;$i<$count;$i++){
			for (my $j=0;$j<$count;$j++){
				if ($row[$i] == $b[$j][$_[1]]){
					my $o=$b[$j][0];
					my $l=$b[$j];
					for (my $k=1;$k<$l;$k++){
						$o=$o."\t".$b[$j][$k];
					}
					print OUT $o."\n";
					$b[$j][$_[1]]=$b[$j][$_[1]]."removed";
				}
			}
		}
		print "2.";
	}else{
		for (my $i=$count-1;$i>=0;$i--){
			for (my $j=$count-1;$j>=0;$j--){
				if ($row[$i] == $b[$j][$_[1]]){
					my $o=$b[$j][0];
					my $l=$b[$j];
					for (my $k=1;$k<$l;$k++){
						$o=$o."\t".$b[$j][$k];
					}
					print OUT $o."\n";
					$b[$j][$_[1]]=$b[$j][$_[1]]."removed";
				}
			}
		}	
		print "3.";
	}
	print "\n";
}


sub RegulateblastQuerrytoloc{
	#If have parametre "+", amend the result
	#if have parametre "x", target is aa
	#defort: first hit is best
	
	
	my $tp_name="start";
	my $CountNo=1;
	my $SeqMax;
	
	my $addpare=0;
	#my $addpare=100;
	
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
        my $T_start = $tp[8];
        my $T_end = $tp[9];
        my $Q_start = $tp[6];
        my $Q_end = $tp[7];
        my $T_len = $tp[12];
		if ($_[1] eq "x"){
			$T_start=$T_start*3-2;	
			$T_end*=3;
			$T_len*=3;
		}
		
        if ($tp_name ne $tp[0]){
                $CountNo = 1;
                $tp_name = $tp[0];
                $SeqMax  = $T_len;
        }else{   
        	$CountNo++;
		}
		
		my $LeftMinus = 0;
        my $RightAdd = 0;
        #if ($_[2] eq "+")  #amend
        {
			if ($T_end < $T_start)   #reverse
        	{  
        		$LeftMinus = $T_end-1;
            	$RightAdd = $T_len-$T_start;
        	}else{         
            	$LeftMinus = $T_start-1;
            	$RightAdd = $T_len-$T_end;
        	}
        }
        
		#output
		#result BYO list
		my $o = "$tp[0]|$tp[1]\t$tp[0]\t";
        #if ($T_end < $T_start)   #reverse
        if ($Q_end < $Q_start)   #reverse
        {
        	$Q_end=$Q_end-$RightAdd;
        	$Q_start=$Q_start+$LeftMinus;
          	$o=$o."$Q_end\t$Q_start\t-";    
            #$Q_start=$Q_start-$LeftMinus-$addpare;
            #$Q_end=$Q_end+$RightAdd+$addpare;
            #$o=$o."complement(".$Q_start."..".$Q_end.")";         	        
        }else  #Forward
        {
            $Q_start=$Q_start-$LeftMinus-$addpare;
            $Q_end=$Q_end+$RightAdd+$addpare;
            $o=$o."$Q_start\t$Q_end\t+";
        }
        print OUT $o."\n";	
	}		
}

sub regulateblasttoloc{
	#If have parametre "+", amend the result
	#defort: first hit is best
	my $tp_name="start";
	my $CountNo=1;
	my $SeqMax;
	while($lineIN=<IN>){	
		chomp($lineIN);
		my @tp = &Getheadinfo($lineIN);
        my $T_start = $tp[8];
        my $T_end = $tp[9];
        my $Q_start = $tp[6];
        my $Q_end = $tp[7];
        my $T_len = $tp[12];

        if ($tp_name ne $tp[0])
        {
                $CountNo = 1;
                $tp_name = $tp[0];
                $SeqMax = $Q_end;
        }else{   
        	$CountNo++;
		}
		
		my $LeftMinus = 0;
        my $RightAdd = 0;
        if ($_[2] eq "+")  #amend
        {
            $LeftMinus = $Q_start-1;
            $RightAdd = $SeqMax-$Q_end;
        }
		#output
		#result BYO list
		#my $o = ">".$tp[0]."|".$tp[3]."|".$CountNo."|".$tp[1].":";
		#result no list   
		my $o = ">".$tp[0]."|".$tp[3]."|".$tp[1].":";
        if ($T_end < $T_start)   #reverse
        {
        	$T_end=$T_end-$RightAdd;
        	$T_start=$T_start+$LeftMinus;
            $o=$o."complement(".$T_end."..".$T_start.")";            
        }else  #Forward
        {
            $T_start=$T_start-$LeftMinus;
            $T_end=$T_end+$RightAdd;
            $o=$o.$T_start."..".$T_end;
        }
        print OUT $o."\n";	
	}		
}
				
close(IN);
close(OUT);
print "progrem over!!\n";

sub Usage{
    printf("$ARGV[0] can not find\n");
}

sub Comparelocal{	#para with two local array
	#return array, [0]-cross state, [1]-cross value
	#input chrom-1 st-1 ed-1 chrom-2 st-2 ed-2
	#print "aaa";
	my @a=@_;
	#judge order
	my @r=(-1,0);
	if ($_[1]>$_[2]){
		$a[1]=$_[2];
		$a[2]=$_[1];
	}
	if ($_[4]>$_[5]){
		$a[4]=$_[5];
		$a[5]=$_[4];
	}
	if ($_[0] ne $_[3]){
		return @r;  #nocross 
	}else{
		#1-2 before corse 1 2-1 before corse 2 3-2 in 1 4-1 in 2 5-same loc 6-no corse 1 before 2 7-no corse 2 before 1
		#same loc
		if ($a[2]==$a[5] && $a[1]==$a[4]){
			$r[0]=5;
			$r[1]=$a[2]-$a[1]+1;
		}
		#1 include 2
		elsif ($a[1]<=$a[4] && $a[2]>=$a[5]){
			$r[0]= 3;
			$r[1]=$a[5]-$a[4]+1;
		}
		#2 include 1
		elsif ($a[2]<=$a[5] && $a[1]>=$a[4]){
			$r[0]= 4;
			$r[1]=$a[2]-$a[1]+1;
		}
		#1 behand 2
		elsif ($a[1]>=$a[4] && $a[1]<=$a[5]){
			$r[0]= 1;
			$r[1]=$a[5]-$a[1]+1;
		}
		#1 before 2
		elsif ($a[2]>=$a[4] && $a[2]<=$a[5]){
			$r[0]= 2;
			$r[1]=$a[2]-$a[4]+1;
		}
		#no corse 1 before 2
		elsif($a[2]<$a[4]){
			$r[0]= 6;
			$r[1]=$a[4]-$a[2]+1;
		}
		#no corse 2 before 1
		elsif($a[1]>$a[5]){
			$r[0]= 7;
			$r[1]=$a[1]-$a[5]+1;
		}
		return @r;
	}	
}


sub GetbestORF{
	my @stopcount;
	my @transleted;
	my $min=9999;
	my $ORF;
	my $minseq;
	my $minloc;
	#forward 123
	for ($i=0;$i<3;$i++){
		my $seq=substr($_[0],$i,);
		$translated[$i]=&TranslateSeq($seq);
		my @sd=&CountCharNo($translated[$i]);		#O-stop codon
		$stopcount[$i]=$sd[0];
		if ($min>$stopcount[$i]){
			$min=$stopcount[$i];
			$minloc=$sd[1];
			$ORF=$i;
			$minseq=$translated[$i];
		}
	}
	#reverse 123
	my $c = &Getcomplementseq($_[0]);
	for ($i=0;$i<3;$i++){
		my $seq=substr($c,$i,);
		&TranslateSeq($seq);
		$translated[$i+3]=&TranslateSeq($seq);
		my @sd=&CountCharNo($translated[$i+3]);
		$stopcount[$i+3]=$sd[0];		
		if ($min>$stopcount[$i+3]){
			$min=$stopcount[$i+3];
			$minloc=$sd[1];
			$ORF=$i+3;
			$minseq=$translated[$i+3];
		}
	}
	if($min>0){
		#print "$_[1]\n";
	}
	my $out="\tORF=".$ORF."\tStopCodon=".$min."\t".$minloc."\n".$minseq."\n";
	return $out;
}

sub CountCharNo{
    my @seq = split('', $_[0]);
    my $seq=@seq;
    my $count=0;
    my $sl;
    for(my $i=0;$i<$seq;$i++){
		if ($seq[$i] eq "\*"){
			$count++;
			$sl=$sl.$i.",";
		}		
	}
	my @r;
	$r[0]=$count;
	$r[1]=$sl;
	
    return @r;	
}

sub TranslateSeq{
    #O-stop codon
	my %codon;
	$codon{"AAA"} = "K";
	$codon{"AAC"} = "N";
	$codon{"AAG"} = "K";
	$codon{"AAT"} = "N";
	$codon{"ACA"} = "T";
	$codon{"ACC"} = "T";
	$codon{"ACG"} = "T";
	$codon{"ACT"} = "T";
	$codon{"AGA"} = "R";
	$codon{"AGC"} = "S";
	$codon{"AGG"} = "R";
	$codon{"AGT"} = "S";
	$codon{"ATA"} = "I";
	$codon{"ATC"} = "I";
	$codon{"ATG"} = "M";
	$codon{"ATT"} = "I";
	$codon{"CAA"} = "Q";
	$codon{"CAC"} = "H";
	$codon{"CAG"} = "Q";
	$codon{"CAT"} = "H";
	$codon{"CCA"} = "P";
	$codon{"CCC"} = "P";
	$codon{"CCG"} = "P";
	$codon{"CCT"} = "P";
	$codon{"CGA"} = "R";
	$codon{"CGC"} = "R";
	$codon{"CGG"} = "R";
	$codon{"CGT"} = "R";
	$codon{"CTA"} = "L";
	$codon{"CTC"} = "L";
	$codon{"CTG"} = "L";
	$codon{"CTT"} = "L";
	$codon{"GAA"} = "E";
	$codon{"GAC"} = "D";
	$codon{"GAG"} = "E";
	$codon{"GAT"} = "D";
	$codon{"GCA"} = "A";
	$codon{"GCC"} = "A";
	$codon{"GCG"} = "A";
	$codon{"GCT"} = "A";
	$codon{"GGA"} = "G";
	$codon{"GGC"} = "G";
	$codon{"GGG"} = "G";
	$codon{"GGT"} = "G";
	$codon{"GTA"} = "V";
	$codon{"GTC"} = "V";
	$codon{"GTG"} = "V";
	$codon{"GTT"} = "V";
	$codon{"TAA"} = "\*";
	$codon{"TAC"} = "Y";
	$codon{"TAG"} = "\*";
	$codon{"TAT"} = "Y";
	$codon{"TCA"} = "S";
	$codon{"TCC"} = "S";
	$codon{"TCG"} = "S";
	$codon{"TCT"} = "S";
	$codon{"TGA"} = "\*";
	$codon{"TGC"} = "C";
	$codon{"TGG"} = "W";
	$codon{"TGT"} = "C";
	$codon{"TTA"} = "L";
	$codon{"TTC"} = "F";
	$codon{"TTG"} = "L";
	$codon{"TTT"} = "F";	
    
    my $l=length($_[0]);
    my $out='';
    for (my $i=0;$i<($l-3);$i+=3){
    	my $a='X';
    	$a=$codon{substr($_[0],$i,3)};
    	$out=$out.$a;
    }
    return $out;
}


sub Getlocinfo{
        $localinfo[0]=substr($_[0],0,index($_[0],":",,));     #chrom
        if ($_[0]=~/complement/)
        {
                $localinfo[1]=1;        #complement=1, else=0
                $localinfo[2]=substr($_[0],index($_[0],"(",,) + 1,index($_[0],"..",,) - index($_[0],"(",,) - 1);  #start
                $localinfo[3]=substr($_[0],index($_[0],"..",,) + 2,index($_[0],")",,) - index($_[0], "..") - 2);  #end
        }else
        {
                $localinfo[1]=0;
                $localinfo[2]=substr($_[0],index($_[0],":",,) + 1,index($_[0],"..",,) - index($_[0],":",,) - 1);  #start
                $localinfo[3]=substr($_[0],index($_[0],"..",,) + 2,length($_[0]));  #end
        }
    return @localinfo;
}

sub Getheadinfo{
	my $test = -1;
	my $li = 0;
	my @samble = (" ", "\t");
	$_=@_;
	my $sn=1;
	if($_>1){
		for (my $i=1;$i<$_;$i++){
			$samble[$i+$sn]="$_[$i]";
			$sn++;
			#print "$_[$i]\t";
		}
	}
	my $samble=@samble;
	#print "$samble\n";
	my @headinfo;
	my $tp;
	
	#chomp($_[0]);
	
	for (my $i = 0; $i < length($_[0]); $i++)
	{
		$tp = substr($_[0],$i,1);
		my $sp=0;
		
		for(my $i=0; $i<$samble;$i++){
			if($tp eq $samble[$i]){
					$sp=1;
			}
		}
		
		if ($sp==1)
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

sub regulatedOUT{
	for (my $i=0; $i<length($_[0]); $i+=80){
		print OUT substr($_[0],$i,80)."\n";	
	}
}

sub Modulenatch{
	#test two group module, same return 1, otherwise -1
	my $r=-1;
	#forward
	my $tp0=$_[0].$_[0];
	my $tp1=$_[1];
	if ($tp0 =~/$tp1/i){
		$r=1;
	} 
	#reverse
	my $tpr0=&Getcomplementseq($_[0]).&Getcomplementseq($_[0]);
	if ($tpr0 =~/$tp1/i){
		$r=1;
	} 
	return $r;	
}


sub Change_fas_to_complement{
	my $seq="";
	my $a=0;
	while($lineIN=<IN>){
		chomp($lineIN);
    if ($lineIN=~/>/){
    	if($a==1){
    		my $cseq=&Getcomplementseq($seq);
    		print OUT "$cseq\n";
    		$seq="";
    	}else{
    		$a=1;
    	}
    	print OUT "$lineIN.complement\n";
    	
    }else{
    	$seq=$seq."$lineIN";		
		}
	}
 	my $cseq=&Getcomplementseq($seq);
  print OUT "$cseq\n";	
	
}


sub Getcomplementseq{
        my @seq = split('', $_[0]);
        my $cseq;
        foreach $nucleotide(reverse(@seq)) 
				{
    	    			if($nucleotide =~ /a/i) {
      	          			$cseq = $cseq."T";
        				} elsif ($nucleotide =~ /t/i) {
                        $cseq = $cseq."A";
        				} elsif ($nucleotide =~ /g/i) {
                        $cseq = $cseq."C";
        				} elsif ($nucleotide =~ /c/i) {
                        $cseq = $cseq."G";
        				} elsif ($nucleotide =~ /n/i) {
                        $cseq = $cseq."N";
         				} elsif ($nucleotide =~ /-/g) {
                        $cseq = $cseq."-"; 
                }          			
    }
        return $cseq;
}

sub CountSymbleLocNumber{
	my @r;
	$r[0]=0;
	for(my $i=0; $i<length($_[1]);$i++){
		if(substr($_[1],$i,1) eq $_[0]){
			#print "$i\n";
			$r[0]++;
			$r[$r[0]]=$i;
		}
	}	
	#print "$r[0]\n";
	return @r;
}

sub ReadAln{
     	my @tpf= &Getheadinfo($_[0], '\/','\\');
     	my $tpf=@tpf;		
	#my $f=substr($tpf[$tpf-1],0,index($tpf[$tpf-1],'.',0));
	my $f=$_[0];
	print "$f\t";
	open(ALN,'<',"$_[0]");
	my $lineALN;
	my @tpname;
	my @tpseq;
	my $count;
	$lineALN=<ALN>;
	while($lineALN=<ALN>){
		$count++;
		chomp($lineALN);	
		my @tp= &Getheadinfo($lineALN);
		my $tp=@tp;
		if($tp == 2){
			if(substr($tp[0],0,1) ne '*'){
				$tpname=@tpname;
				my $add=1;
				for(my $i=0; $i<$tpname; $i++){
					if($tpname[$i] eq $tp[0]){
						$tpseq[$i]=$tpseq[$i].$tp[1];
						$add=0;
						print ".";
						last;
					}
				}
				if($add==1){
					print "$tp[0] ";
					$tpname[$tpname]=$tp[0];
					$tpseq[$tpname]=$tp[1];			
				}
			}
		}
	}
	close(ALN);
	my $r;
	$tpname=@tpname;
	for(my $i=0; $i<$tpname; $i++){
		$r=$r.">".$f."_".$tpname[$i]."\n".$tpseq[$i]."\n";
	}
	print OUT "\n";
	return $r;
}

sub FasttoArray{	
	my @headinfo;
	my @r;
	my $count=0;
	my $tpseq;
	my $tpname;
	my $canadd = 0;
	open(DB,'<',"$_[0]");
	my $lineDB;
	while($lineDB=<DB>){
		chomp($lineDB);
        if(index($lineDB,'>')>=0)
        {        	      	
        	if ($canadd > 0){
        		$r[$count][0] = $tpname;
                $r[$count][1] = $tpseq;
                $count++;
            }
			@headinfo = Getheadinfo($lineDB);
			$tpname=substr($headinfo[0],1,length($headinfo[0]));              
            $tpseq = '';
                
    	}else{
        	$tpseq = $tpseq.$lineDB;
        	$canadd = 1;
    	}
	}
    $r[$count][0] = $tpname;
    $r[$count][1] = $tpseq;
    return @r;
}
