#!/usr/bin/perl
my $filename = $ARGV[0]; 

open(my $fh, '<', $filename)
  or die "Could not open file '$filename' $!";

my $line_num  = 0;
while (my $line= <$fh>) {
    if ($line_num % 4 == 0) {
        print(">". substr($line, 1));
        $line_num = 0;
    }
    if ($line_num  == 1) {
        print($line);
    }

    $line_num++;
}

