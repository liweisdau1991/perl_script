#!/usr/bin/perl
use Bio::SeqIO; 
my ($file1,$file2)=@ARGV; 
my $seqin = Bio::SeqIO -> new (-format => 'fastq',-file => $file1); 
my $seqout = Bio::SeqIO -> new (-format => 'fasta',-file => ">$file2"); 

while (my $seq_obj = $seqin -> next_seq) 
{ 
   $seqout -> write_seq($seq_obj); 
}

