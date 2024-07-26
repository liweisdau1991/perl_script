#!/usr/bin/perl
use Getopt::Std;
getopts("d:");

my $dir=$Getopt::Std::opt_d;

if (! $dir ) {
print "\nUsage: perl $0 -d dir\n";
exit;
}
opendir(DIR,$dir);

my @files=grep(/(\S+)\.fa\.tmp/,readdir(DIR));

system("mkdir out");
foreach my $file (sort @files) {
	chomp($file);
	my($id)=$file=~/(\d+)/;
	#system("mafft --thread 10 $dir\/$file > \.\/out/$id.out");
	print "mafft --thread 10 $dir\/$file > \.\/out/$id.out\n";
	}
closedir(DIR);
#print "\n-------------Finished!\n";

