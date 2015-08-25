#!/usr/bin/perl -w

use strict;

#Check that number of arguments is correct
if(scalar @ARGV != 1) {
    die "$0 <.Repr2>\n";
}

#INPUT:
#1) File with representatives (.Repr2)
my $infile = $ARGV[0];

#OUTPUT:
#1) STDOUT that is redirected to a file (.formatcdhit) by the pipeline

open INFILE, $infile or die "File $infile cannot be opened\n";

while(<INFILE>)
{
	chomp $_;
        if ($_ =~ />/)
        {
        	my @line = split(" ", $_);

                foreach my $l (@line)
                {
                	print $l."*";
                }

        } else
        {

        	print "\n$_\n";
        }
}
close(INFILE);