#!/usr/bin/perl -w
use strict;

if(scalar @ARGV != 2) {
    die "$0 <.CDHIT file> <.Repr2 or .Members2 file>\n";
}


#INPUT:
#1) File of non-redundant representatives, which is the output of CD-HIT (.CDHIT)
my $nr = $ARGV[0];

#OUTPUT:
#1) STDOUT which is redirected to a file (.ReprFinal or .MembersFinal)


my @clusterID;

open INFILE, $nr or die "Could not open infile $nr\n";

while(<INFILE>)
{
	chomp $_;
        if ($_ =~ /^>/)
        {
        	my @line = split(/\*/, $_);
        	push (@clusterID, $line[7]);
        }
}

close(INFILE);

#Infile $members is file with all clusters members previous to redundancy filter (from pipeline)
my $members = $ARGV[1];
my $flag = 0;

open INFILE, $members or die "Could not open infile $members\n";

while(<INFILE>)
{

	if ($_ =~ m/^>/)
        {

        	foreach my $i (@clusterID)
        	{
        		if ($_ =~ m/\b$i\b/)
                	{
                		print $_;
                                $flag = 1;
                                last;
               	 	}
                }
        } else {
        	if ($flag == 1)
                {
                	print $_;
                        $flag = 0;
                }
        }

}