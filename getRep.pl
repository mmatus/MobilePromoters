#!/usr/bin/perl -w

use strict;

#Check that number of arguments is correct
if(scalar @ARGV != 2) {
     die "$0 <.filterUpParalogs file> <.Database file>\n";
}

#INPUT:
#1) File with clusters that have been filtered for paralogy of upstream CDSs (.filterUpParalogs)
my $infile = $ARGV[0];
#2) Database with all sequences (.Database)
my $database = $ARGV[1];

#OUTPUT:
#1) STDOUT which is redirected to two files, one with representatives (.Repr1) and one with all members of a cluster (.Members1)


my @fileName = split (/\./, $infile);
my $run_name = $fileName[0];
my $outfile_Rep = $run_name.".Repr1";
my $outfile_Mem = $run_name.".Members1";
my @pids_Rep;
my @pids_Mem;




open INFILE, $infile or die "Cannot open $infile\n";

while (<INFILE>)
{
	chomp $_;
        my @line = split (/\t/, $_);

        push (@pids_Rep, $line[11]);

        my @members = split(" ", $line[12]);
        push (@pids_Mem, @members);

}

close(INFILE);


foreach my $p (@pids_Rep)
{
	print "Looking for $p\n";
        `grep ">$p" $database -A 1 -m 1 >> $outfile_Rep`;

}

foreach my $p (@pids_Mem)
{
	print "Looking for $p\n";
        `grep ">$p" $database -A 1 -m 1 >> $outfile_Mem`;
}