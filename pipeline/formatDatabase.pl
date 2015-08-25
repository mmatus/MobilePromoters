#!/usr/bin/perl -w

use strict;


if(scalar @ARGV != 2) {
     die "$0 <.filterUpParalogs file> <.ReprEXT1 or .MembersEXT1 file>\n";
}

my $infile1 = $ARGV[0];
my $infile2 = $ARGV[1];
my %clusters;
my $cluster_id = 0;
open CLUSTERS, $infile1 or die "File $infile1 cannot be opened\n";
my %plus;
my %minus;

while(<CLUSTERS>)
{
	chomp $_;
        my @line1 = split (/\t/, $_);
        $cluster_id++;
        $clusters{$line1[12]} = $line1[0]." Cl-".$cluster_id;
	 #print "$line1[12] $clusters{$line1[12]}\n";


}

close(CLUSTERS);



open MEMBERS, $infile2 or die "File $infile2 cannot be opened\n";

while(<MEMBERS>)
{
	chomp $_;
        if ($_ =~ m/>/)
        {
        	$_ =~ s/>//;
        	my @line = split(/\t/, $_);
                $line[5] =~ s/\W/_/g;
                chop $line[5];
                
                my $pid = $line[0];
		 
                foreach my $k (keys %clusters)
                {
                	if ($k =~ m/\b$pid\b/)
                        {
                        	print ">@line ";
                                print $clusters{$k};

                        }
                }
        }  else
        {
        	print "\n$_\n";
        }


}
close(MEMBERS);

