#!/usr/bin/perl

use strict;

my $base = "/Users/Rohandinho/Desktop/Evolutionary\ rate\ project/data/";
my $mut_prefix = "mutator_40K_diffs/";
my $nonmut_prefix = "non-mutator_40K_diffs/";
my $mutdir = $base.$mut_prefix;
my $nonmutdir = $base.$nonmut_prefix;

count_mobs($mutdir);
count_mobs($nonmutdir);

##input: a directory containing gb files.
## the total number of mob diff lines in the line.
sub count_mobs {
    my $diffdir = shift;
    opendir(my $diff_dirh, $diffdir) or die $!;
    while (my $curfile = readdir($diff_dirh)) {
	
	next unless ($curfile =~ m/.gd$/);
	## now open the current diff file.
	my $curpath = $diffdir.$curfile; 
	open(my $curfile_fh, $curpath);
	print "#",$curfile, "\n";
	my $total_mobs = 0;
	while (<$curfile_fh>) {
	    ##skip everything except mobile element annotations.
	    my $curline = $_;
	    next unless $curline =~ m/^MOB/;
	    print $curline;
	    $total_mobs = $total_mobs+1;
	}
	print "##TOTAL MOBS: ", $total_mobs, "\n";
	print "//\n";
    }

}
