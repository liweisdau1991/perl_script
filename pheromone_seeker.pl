#!/usr/bin/perl -w
#Copyright 2016 Dr. Alija Bajro Mujic (amujic@gmail.com)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

# This script will search a multiple sequence nucleotide fasta file for open reading frames that contain CaaX motifs typical
# of small phermone precursor genes.  Both forward and reverse strands are searched.

use strict;
use Bio::SeqIO;
use Bio::DB::Fasta;

sub trim($);

my $usage="This script will search a multiple sequence nucleotide fasta file for open reading frames that contain CaaX motifs typical of small phermone precursor genes.  Both forward and reverse strands are searched\n\nusage: pheromone_seeker.pl inputfile";
#require file input from the command line
die("$usage\n") unless(@ARGV > 0); 
my $input=$ARGV[@ARGV-1];

#declare the pheromone motif we are looking for...
#pattern is structured as such: M|L protein-body C (V|L|I|M|){2} X
my $pheromone_pattern="ATG ((?!TAA|TAG|TGA)[TCAG]{3})+ TG[TC] (GT[TCAG] | CT[TACG] | AT[TCAG]){2} [TCAG]{3} (TAA|TAG|TGA)";

#declare the sequence handle and container stuff
my $active_seq;
my $seq_obj_revcom;
my $active_seq_revcom;
my $sequence_handle=Bio::SeqIO->new( -file=>$input, -format=>"Fasta");
my $db = Bio::DB::Fasta->new($input);
my $temp_header;
my @active_header;

#ok, now loop over all of the sequences in the file and output the appropriate stuff
while(my $seq_obj=$sequence_handle->next_seq){
	$temp_header=$db->header($seq_obj->id);
	chomp($temp_header);
	@active_header=split(/_/, $temp_header);
	$active_seq= $seq_obj->seq;
	$seq_obj_revcom=$seq_obj->revcom;
	$active_seq_revcom=$seq_obj_revcom->seq;
	while($active_seq =~ /$pheromone_pattern/gix){
		print ">". $seq_obj->id . "_". (length($`) + 1) . "_" . (length($`) + 1 + length($&)) . "_" . "+" ."\n" . $& . "\n";
	}
	while($active_seq_revcom =~ /$pheromone_pattern/gix){
		print ">" .  $seq_obj->id . "_" . (length($`) + 1) . "_" . (length($`) + 1 + length($&)) . "_" . "-"."\n" . $& . "\n";
	}
}#end while


# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

