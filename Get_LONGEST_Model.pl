#!/usr/bin/perl

my $file = $ARGV[0]; # protein.fasta

my($gid,$pid,%SEQ,%LEN,%LONGEST,%HIT_ID,%HIT_SEQ);

open IN,'<',$file;
while(<IN>){
  chomp;
  if(/^>/){
  s/^>//;
  ($pid) = (split/\s+/,$_)[0];
   $SEQ{$pid} = '';
   $LEN{$pid} = 0;
  }
  else{
   s/\*//g;
   $SEQ{$pid}.=$_;
   $LEN{$pid}+=length($_);
  }
}
close IN;

foreach (keys %LEN){
  chomp;
  ($gid) = (split/\./,$_)[0];
  if($LONGEST{$gid} < $LEN{$_}){
   $LONGEST{$gid} = $LEN{$_};
   $HIT_ID{$gid} = $_;
   $HIT_SEQ{$gid} = $SEQ{$_};
  }
}

foreach (keys %LONGEST){
  chomp;
  print '>'.$HIT_ID{$_}."\n".$HIT_SEQ{$_}."\n";
}
