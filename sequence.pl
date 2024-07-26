#!/home/liwei/SoftWare/perl/bin/perl
# count sequence infomation
# 2009-02-23 BGI xulizhi <xulz@genomics.org.cn>

use strict;
use warnings;

use FindBin qw($Bin $Script);
use lib "$Bin/lib";
use SEQUENCE qw(ReadFa);
use MATH qw(Max Sum NumN50 NumN90);

#get input information
use Getopt::Long;
my %opts = (length=>1000);
GetOptions(\%opts,"length:i","contig:s","scaff:s","gaplis:s","table","help");

if(@ARGV < 1 or exists $opts{help}){
    Usage(),exit(1);
}

# pretreatment
printf("%20s %9s %9s %9s %12s %12s %12s %12s %9s %12s %12s %12s %6s %9s %9s %9s %9s %12s\n","#SeqName","#TotalScf","#MinLen","#ScfNum","#ScfLen","#ScfN50","#ScfN90","#ScfMax","#CtgNum","#CtgLen","#CtgN50","#CtgMax","#GC(%)","#InGapNum","#InGapLen","#InGapMax","#mScfNum","#mScfLen") if(exists $opts{table});
if(exists $opts{scaff} and $opts{scaff} ne ""){
    open(SCF,">$opts{scaff}") or die "OpenError: Write \'$opts{scaff}\',$!";
}
if(exists $opts{contig} and $opts{contig} ne ""){
    open(CTG,">$opts{contig}") or die "OpenError: Write \'$opts{contig}\',$!";
}
if(exists $opts{contig} and exists $opts{gaplis} and $opts{gaplis} ne ""){
    open(GAP,">$opts{gaplis}") or die "OpenError: Write \'$opts{gaplis}\',$!";
}

foreach my $each (@ARGV){
    my (@list_seq,@scf_len,@ctg_len);
    my %hash_stat = (gap_max=>0,gap_len=>0,ctg_max=>0,ctg_num=>0,scf_max=>0,ctg_len=>0,scf_num=>0,gc=>0,gap_num=>0,scf_len=>0,less_num=>0,less_len=>0);
    ReadFa(\@list_seq,$each);
    foreach my $seq (@list_seq) {
        my $len = length $seq->[1];
        # count less length seq infomation
        if ($len < $opts{length}){
            $hash_stat{less_num} ++;
            $hash_stat{less_len} += Sum(map(length($_),split(/N+/,$seq->[1])));
            next;
        }
        if(exists $opts{scaff}){
            my $seq_tmp = $seq->[1];
            $seq_tmp =~ s/(\w{50})/$1\n/g;
            print SCF ">$seq->[0]\n$seq_tmp\n";
        }
        # count scaffold infomation
        $hash_stat{scf_num} ++;
        $hash_stat{scf_len} += $len;
        if(exists $hash_stat{scf_max}){
            $hash_stat{scf_max} = Max($hash_stat{scf_max},$len);
        }else{
            $hash_stat{scf_max} = $len;
        }
        push(@scf_len,$len);
        # count gap information
        my @gap = map(length($_),split(/[ATGC]+/,$seq->[1],));
        $hash_stat{gap_num} += grep($_ != 0,@gap);
        $hash_stat{gap_len} += Sum(@gap);
        if(exists $hash_stat{gap_max}){
            $hash_stat{gap_max} = Max(@gap,$hash_stat{gap_max});
        }else{
            $hash_stat{gap_max} = Max(@gap);
        }
        # count contigs information
        my @ctg=map(length($_),split(/N+/,$seq->[1]));
        if(exists $opts{contig}){
            my @ctg_tmp = split(/N+/,$seq->[1]);
            for(my $i=1;$i<=@ctg;$i++){
                $ctg_tmp[$i-1] =~ s/(\w{50})/$1\n/g;
                print CTG ">Contig".($hash_stat{ctg_num}+$i)."\n$ctg_tmp[$i-1]\n";
                print GAP "$seq->[0]\tContig".($hash_stat{ctg_num}+$i-1)."\tR\t$seq->[0]\tContig".($hash_stat{ctg_num}+$i)."\tL\t$gap[$i-1]\n" if (exists $opts{gaplis} and $i>1);
            }
        }
        $hash_stat{ctg_num} += @ctg;
        $hash_stat{ctg_len} += Sum(@ctg);
        push(@ctg_len,@ctg);
        if(exists $hash_stat{ctg_max}){
            $hash_stat{ctg_max} = Max($hash_stat{ctg_max},@ctg);
        }else{
            $hash_stat{ctg_max} = Max(@ctg);
        }
        # count gc
        $seq->[1] =~ s/[ATN]//g;
        $hash_stat{gc} += length($seq->[1]);
    }
    # print out
    if(exists $opts{table}){
        # "#SeqName","#MinLen","#mScfNum","#mScfLen","#ScfN50","#ScfMax","#CtgNum","#CtgLen","#CtgN50","#CtgMax","#GC(%)","#GapNum","#GapLen","#GapMax"
        printf("%20s %9d %9d %9d %12d %12d %12d %12d %9d %12d %12d %12d %3.2f %9d %9d %9d %9d %12d\n",
               $each,$hash_stat{scf_num}+$hash_stat{less_num},$opts{length},$hash_stat{scf_num},$hash_stat{scf_len},NumN50(@scf_len),NumN90(@scf_len),$hash_stat{scf_max},$hash_stat{ctg_num},
               $hash_stat{ctg_len},NumN50(@ctg_len),$hash_stat{ctg_max},$hash_stat{gc}/$hash_stat{ctg_len}*100,$hash_stat{gap_num},$hash_stat{gap_len},$hash_stat{gap_max},$hash_stat{less_num},$hash_stat{less_len});
    }else{
        printf("# $each Sequence information\n  Scaffold stat:\n");
        printf("  \t%9s %9s %9s %12s %12s %12s %12s %9s %12s\n","TotalNum","MinLen","ScfNum","ScfLen","ScfN50","ScfN90","ScfMax","mScfNum","mScfLen");
        printf("  \t%9d %9d %9d %12d %12d %12d %12d %9d %12d\n",$hash_stat{scf_num}+$hash_stat{less_num},$opts{length},$hash_stat{scf_num},$hash_stat{scf_len},NumN50(@scf_len),NumN90(@scf_len),$hash_stat{scf_max},$hash_stat{less_num},$hash_stat{less_len});
        printf("  Contigs stat:\n");
        printf("  \t%9s %12s %12s %12s %12s %6s\n","CtgNum","CtgLen","CtgN50","CtgN90","CtgMax","GC(%)");
        printf("  \t%9d %12d %12d %12d %12d %3.2f\n",$hash_stat{ctg_num},$hash_stat{ctg_len},NumN50(@ctg_len),NumN90(@ctg_len),$hash_stat{ctg_max},$hash_stat{gc}/$hash_stat{ctg_len}*100);
        printf("  Insert Gap stat:\n");
        printf("  \t%9s %12s %9s\n","GapNum","GapLen","MaxGap");
        printf("  \t%9d %12d %9d\n\n",$hash_stat{gap_num},$hash_stat{gap_len},$hash_stat{gap_max});
    }
}

close SCF if (exists $opts{scaff});
close CTG if (exists $opts{contig});

sub Usage{
    print qq(
        Usage:
            perl $0 [options] <seq.fa> ...
        Options:
            -contig <str>   Output Contigs sequence which longer than minlength
            -scaff  <str>   Output Scaffold sequence which longer than minlength
            -gaplis <str>   Output Scaffold insert gaplis which longer than minlength
            -length <int>   minimum length of count [default 1000]
            -table  < - >   print stat info as table format
            -help   < - >   show this help information
        Introductions:
            options '-contig' or '-scaff' only use one input sequence file.
            options '-gaplis' must use '-contig'.
    \n);
}
