#!/usr/bin/perl -w
use strict;

#Check that user introduces necessary data.
if(scalar @ARGV != 2) {
     die "$0 <.Clusters file> <paralogs cut-off>\n";
}

#INFILE:
#1) Clusters file (.Clusters)
my $infile = $ARGV[0];

#OUTFILES:
my @fileName = split (/\./, $infile);
my $run_name = $fileName[0];
#1) File that keeps track of how many clusters we have of only promoters, only CDS, and mixed (for analysis of numbers).
my $sumFile = $run_name.".Sum";
#2) File that prints all clusters with paralogy information (as control, and if there's interest).
my $paralogsFile = $run_name.".ClustersParalogs";
#3) Print to STDOUT the clusters that pass the target CDS paralogy cut-off

# Target CDS paralogy cut-off.
my $cutoff = $ARGV[1];

#Variables declaration.
my %promClusters;
my %cdsClusters;
my %paralogs;
my %upParalogs;
my %otherInfo;
my %numMembers;
my %representative;


#Open files
#INFILE 1
open CLUSTERFILE, $infile or die "cannot open $infile\n";
#OUTFILE 1
open SUMA, ">$sumFile" or die "cannot open $sumFile\n";
#OUTFILE 2
open PARALOGS, ">$paralogsFile" or die "cannot open $paralogsFile\n";

while (<CLUSTERFILE>)
{
	chomp $_;
        my @line = split(/\t/,$_);

	if ($_ =~ /_prom/ && $_ =~ /CDS/) #Skip the mixed clusters
        {
        	print SUMA "MIXED\n";
                $cdsClusters{$line[7]} = $line[0]; #I need to store it to filter properly paralogy

        } elsif ($_ =~ /_prom/)#Store the clusters that have only promoters
        {

        	$promClusters{$line[7]} = $line[0];
                $numMembers{$line[7]} = $line[1];
        	$otherInfo{$line[7]} = "$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]";
                $representative{$line[7]} = $line[6];
                print SUMA "PROM\n";

        } else #Store the clusters that have only CDS
        {
         	$cdsClusters{$line[7]} = $line[0];
                print SUMA "CDS\n";

        }
}

close (CLUSTERFILE);
close (SUMA);

my @proms = keys(%promClusters);
my %clustersNumber;
my %upClustersNumber;

foreach my $i (@proms)
{
	$paralogs{$i} = 0;
        $upParalogs{$i} = 0;
        my @members = split(" ", $i);
        my $clusterId;
        $clustersNumber{$i} = 0;
        $upClustersNumber{$i} = 0;
        my %count;
        my %upCount;

        foreach my $j (@members)
        {
        	$j =~ m/(\d+)_\w+/;
                my $pid = $1;
                my $cdsHit = $pid."_CDS";
                my $upHit = $pid."_upCDS";
                my $grep_pid = `grep $pid $infile | cut -f 8`; #Reduce the search space
                my @grep_pid = split(/\n/, $grep_pid);

        	foreach my $c(@grep_pid)
                {
                	chomp $c;
                        $clusterId = $cdsClusters{$c};

                        if ($c =~ m/\b$cdsHit\b/)
                        {
                               	$count{$clusterId}++;

                        }

                        if ($c =~ m/\b$upHit\b/)
                        {
                               	$upCount{$clusterId}++;

                        }
                }


	}

        #CDS paralogs filtering

        foreach my $k (keys %count)
        {
        	if ($count{$k} > 1)
                {
                	$paralogs{$i} += $count{$k};
                        $clustersNumber{$i} ++;
                }
        }



        my $paralogsPerc = $paralogs{$i}/$numMembers{$i};

        #Upstream CDS paralogs filtering
        foreach my $k (keys %upCount)
        {
		if ($upCount{$k} > 1)
                {
                	$upParalogs{$i} += $upCount{$k};
                        $upClustersNumber{$i} ++;

                }
        }


        my $upParalogsPerc = $upParalogs{$i}/$numMembers{$i};

        #Print paralogy information for all clusters as a check-up

        print PARALOGS "$otherInfo{$i}\t";
        print PARALOGS "$paralogs{$i}\t$clustersNumber{$i}\t";
        print PARALOGS sprintf ("%.2f",$paralogsPerc);
        print PARALOGS "\t$upParalogs{$i}\t$upClustersNumber{$i}\t";
        print PARALOGS sprintf ("%.2f",$upParalogsPerc);
        print PARALOGS "\t$representative{$i}\t$i\n";


        if ($paralogsPerc <= $cutoff)  #Print only the filtered ones to continue with the pipeline
        {
                print "$otherInfo{$i}\t";
		print "$paralogs{$i}\t$clustersNumber{$i}\t";
                printf ("%.2f",$paralogsPerc);
                print "\t$upParalogs{$i}\t$upClustersNumber{$i}\t";
        	printf ("%.2f",$upParalogsPerc);
                print "\t$representative{$i}\t$i\n";
        }


}