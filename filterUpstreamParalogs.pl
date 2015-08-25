#!/usr/bin/perl -w


use strict;

#Check that the number of arguments is correct.
if(scalar @ARGV != 2) {
     die "$0 <.filterParalogs file> <upstream paralogy cut-off>\n";
}

#INFILE:
#1) File with clusters that have been filtered for paralogy of the downstream CDS (.filterParalogs)
my $infile = $ARGV[0];
#2) Cut-off for upstream paralogy (0.0-1.0) This variable is taken from infile paramaters.txt
my $cut_off = $ARGV[1];


#OUTFILES:
#1) STDOUT with clusters that have been filtered for paralogy of the upstream CDS.
#The STDOUT is redirected to a file (.filterUpParalogs) by the pipeline.



open INFILE, $infile or die "cannot open $infile\n";

while (<INFILE>)
{
	chomp $_;
	my @line = split(/\t/, $_);
        my $up_score = $line[10];

        if ($up_score <= $cut_off)
        {
        	print "$_\n";

        }
}
close (INFILE);