#!/usr/bin/perl -w

use strict;


if(scalar @ARGV != 5) {
     die "$0 <blastfile> <score cut-off> <min length> <min identity> <run name>\n";
}

my $blastFile = $ARGV[0];
my $cutoff = $ARGV[1];
my $MINLENGTH = $ARGV[2];
my $MINIDENT = $ARGV[3];
my $RUN_NAME = $ARGV[4];
my @fileName = split(/\./, $blastFile);
my $edgeFile = $fileName[0].".Edges";

#INPUT:
#1) Open BLAST file to parse results.
open BLASTFILE, $blastFile or die "cannot open $blastFile\n";


#OUTPUT:
#1) Create the edges file for each genome (to be used by NetClust).
open EDGELIST, ">$edgeFile" or die "Cannot open $edgeFile\n";
#2) Print to STDOUT the Clusters obtained with Netclust and some more information about them.
#3) Print to STDERR all promoter hits information: strands, length, identity, score (for filterOrientation program and downstream analysis of results)


my $query = "";
my $subject = "";
my $rawscore = 0;
my $len;
my $id;
my $strandQ;
my $strandS;
my $flag;
my %goodPairs;

while(<BLASTFILE>)
{
    if(/Query= (\S+)/)
    {
        $query = $1;
    } elsif (/>(\S+)/)
    {
        $subject = $1;
        $flag = 0; #Flag that indicates that a new subject has been found.

    } elsif (/Score = *([.0-9]+) bits \((\d+)\),/)
    {
        $rawscore = $2;

    } elsif (m#Identities\s=\s(\d+)/(\d+)\s\((\d+)%\)#)
    {
    	$len = $1;
        $id = $3;

    } elsif (m#Strand\s=\s(\w+)\s/\s(\w+)#)
    {
    	$strandQ = $1;
        $strandS = $2;

    	if ($query ne $subject) #Skip self-hits.
        {
        	$flag++;

                #Store relevant hits.

                if ($flag == 1) #Store one hit.

                {
        		if ($query =~ m/prom/)  #Promoter query. Only hits with min length and identity are stored.
        		{
                        	#Print all promoter hits in a file, later to be used by other programs (just the first hit)
                        	print STDERR "$query\t$strandQ\t$subject\t$strandS\t$len\t$id\t$rawscore\n";

                                #Store only hits with min length and identity
           			if ($len >= $MINLENGTH && $id >= $MINIDENT)
                		{
                       			$goodPairs{$query}{$subject} = $rawscore;
                		}

        		} else   #CDS query. Don´t store self-hits of CDS and upCDS.
        		{
                        	#Possible redundant entries are +/- one PID number
                                $query =~ /(\d+)_/;
                                my $pid = $1;
                                my $query_plus = $pid+1;
                                my $query_less = $pid-1;

                        	if ($subject =~ m/$query_plus/ or $subject =~ /$query_less/)
                                {
                                	#Redundant CDS/upCDS pairs
                                } else
                                {
                                	$goodPairs{$query}{$subject} = $rawscore;
                                }
        		}

                        #Print relevant hits in the edges and orientation files.
                        if (defined($goodPairs{$query}{$subject}))
        		{
        			print EDGELIST "$query $subject $goodPairs{$query}{$subject}\n"; #Print the raw scores of all hits with min length and min identity per genome.
        		}
                }
        }
    }

}

close (BLASTFILE);
close (EDGELIST);

#Use netclust
if (-s $edgeFile) #Need to check since some genomes did not have any paralogs and therefore do not produce an edges file.
{

`netindex $edgeFile`;
my @clusterLines = `netclust $edgeFile O2 S $cutoff`;

foreach my $clusterLine (@clusterLines) {
    chomp $clusterLine;
    my ($clusterID, $clusterCount, $clusterMemberField) = split(/\t/, $clusterLine, 3);
    my @clusterMembers = split(/ /, $clusterMemberField);
    if(scalar @clusterMembers != $clusterCount)
    {
        die "No valid cluster line: $clusterLine\n";
    }
    my ($cEdges, $maxDegree, $representative) = &countEdges(@clusterMembers);
    my $maxEdges = $clusterCount * ($clusterCount - 1) / 2;
    my $filter = $cEdges/$maxEdges;

    print "$clusterID\t$clusterCount\t$maxDegree\t$cEdges\t$maxEdges\t";
    printf("%.2f", $filter);
    print "\t$representative\t@clusterMembers\n";

}

}

#Subroutine used to count the number of edges in the clusters and to get a representative that has the most connections in the cluster.
sub countEdges {
    my $cEdges = 0;
    my $maxDegree = 0;
    my $representative;

    for(my $i = 0; $i < scalar (@_); $i++) {
        my $degree = 0;
        for(my $j = 0; $j < $i; $j++) {
             if ( (defined($goodPairs{$_[$i]}{$_[$j]})) || (defined($goodPairs{$_[$j]}{$_[$i]})) )
             {
                $degree++;
             }
        }

        if ($degree > $maxDegree) {
            $maxDegree = $degree;
            $representative = $_[$i];
        }

        $cEdges += $degree;
    }

    return ($cEdges, $maxDegree, $representative);
}