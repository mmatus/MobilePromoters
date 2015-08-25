#!/usr/bin/perl -w

use strict;

#Check that the number of arguments is correct
if(scalar @ARGV != 7) {
     die "$0 <chromosome.fna> <min promoter size> <max promoter size> <min cds size> <max cds size> <run name> <set-off>\n";
}

#INPUT:
#1) Name of file of a chromosome in FASTA format (.fna). This variable is defined by the pipeline starting from the directory name in parameters.txt.
#2) Min/max promoter size; min/max CDS size; run name and set-off. All these variables are taken from the infile parameters.txt.

#OUTPUT:
#1) STDOUT is redirected to a file (.fas) that contains all promoters, CDS and upCDS extracted sequences (1 per genome).
#2) STDERR is redirected to a global database (.Database) that contains all sequences from all genomes.


#Declare of global variables
my $CHROMOSOME = $ARGV[0];
my $MIN_PROMOTER = $ARGV[1];
my $MAX_PROMOTER = $ARGV[2];
my $MIN_CDS_LENGTH = $ARGV[3];
my $MAX_CDS_LENGTH = $ARGV[4];
my $RUN_NAME = $ARGV[5];
my $SET_OFF = $ARGV[6];
my @data = split (/\//, $CHROMOSOME);
my $genome = $data[1];
my $chrom = $data[2];
my $sequence = "";
my $header;
my $typeDNA;

#*******************************SUBROUTINES************************************
# Sub to reverse-complement the promoter and CDS sequences if the gene is in the minus strand.
# Necessary when performing the BLAST search only in the plus (or minus) strand later on.
sub RevComp($)
{
        my $Sequence = shift;

        $Sequence = reverse($Sequence);
	$Sequence =~ tr/ACGTacgt/TGCAtgca/;

        return ($Sequence);
}



#Sub to check the promoter size of the first gene of a DNA molecule when it is in the plus strand (so there's no previous gene).
sub checkFirstProm ($)
{
	my $start = shift;
        my $promoterSequence;
        my $upstreamRegion = $start; #$start - 0 which is $start
        print LENGTHS "$upstreamRegion\n"; #modified

        if ($start < $MIN_PROMOTER)
        {
        } elsif ($start < $MAX_PROMOTER)
        {
        	$promoterSequence = substr($sequence, 0, $start); #check examples !!!
        } else
        {
         	$promoterSequence = substr($sequence, $start - $MAX_PROMOTER, $MAX_PROMOTER); #as original, seems to work.
        }
        return ($promoterSequence);
}

#Sub to check the promoter size of the last gene of a DNA molecule when it is in the minus strand (so there's no previous gene).
sub checkLastProm ($)
{
	my $end = shift;
        my $promoterSequence;
        my $sequenceLength = length($sequence);
        my $upstreamRegion = $sequenceLength - $end;
        print LENGTHS "$upstreamRegion\n"; #modified

        if ($sequenceLength - $end < $MIN_PROMOTER)
        {
        } elsif ($sequenceLength - $end < $MAX_PROMOTER)
        {
        	$promoterSequence = substr($sequence, $end+1, $sequenceLength - $end);#Modified. Seems to work.
        } else
        {
        	$promoterSequence = substr($sequence, $end+1, $MAX_PROMOTER);  #Modified. Seems to work.
        }

	return($promoterSequence);
}

#Subroutine to check promoter sequences of middle genes in the plus strand.
sub checkPlusProm ($$)
{
	my $start = $_[0];
	my $prevGeneData = $_[1];
        my ($prevPID, $prevStart, $prevEnd, $prevStrand, $prevCOG, $prevProduct) = @{$prevGeneData};
        my $promoterSequence;
        my $upstreamRegion = $start - $prevEnd-1;#modified!
        print LENGTHS "$upstreamRegion\n"; #Store all promoter region sizes

        if($prevStrand eq "-") #Case 1. "Previous" gene is in the minus strand. So the two genes are divergent.
        {
		if($upstreamRegion < (2*$MIN_PROMOTER + $SET_OFF)) #Skip all regions smaller than 2*MIN (200) + set-off of previous gene (50 nt)
                {
                } elsif ($upstreamRegion < (2*$MAX_PROMOTER + $SET_OFF)) #Take promoters between 2* (MIN and MAX)
                {
                	$upstreamRegion = int ($upstreamRegion/2); #Divide the available promoter region into 2. Take only integer numbers as result.
                        $promoterSequence = substr($sequence, $start - $upstreamRegion, $upstreamRegion);  #as original. Seems to work.
                } else  #Else take the MAX
                {
                	$promoterSequence = substr($sequence, $start - $MAX_PROMOTER, $MAX_PROMOTER); #as original. Seems to work.
                }
	} else #Case 2. "Previous" gene is in the plus strand so both genes are in the same orientation.
        {
               if($upstreamRegion < $MIN_PROMOTER) #You only need the minimum (100 nt)
               {
               } elsif ($upstreamRegion < $MAX_PROMOTER) #Take promoters between MIN and MAX
               {
                        $promoterSequence = substr($sequence, $start - $upstreamRegion, $upstreamRegion);#as original.
	       } else #Else take the MAX
               {
               		$promoterSequence = substr($sequence, $start - $MAX_PROMOTER, $MAX_PROMOTER); #as original. Seems to work.
               }
	}

        return ($promoterSequence);

}

#Subroutine to check promoter sequences of middle genes in the minus strand.
sub checkMinusProm ($$)
{
	my $end = $_[0];
        my $prevGeneData = $_[1];
        my ($prevPID, $prevStart, $prevEnd, $prevStrand, $prevCOG, $prevProduct) = @{$prevGeneData};
        my $promoterSequence;
        my $upstreamRegion = $prevStart - $end -1; #modified!
        print LENGTHS "$upstreamRegion\n"; #Store all promoter region sizes

        if($prevStrand eq "-") #Case 1. "Previous" gene is in the minus strand so both genes are in the same orientation.
        {

        	if($upstreamRegion < $MIN_PROMOTER)#Skip all regions smaller than MIN.
               {
               } elsif ($upstreamRegion < $MAX_PROMOTER)#Take promoters between MIN and MAX.
               {
                        $promoterSequence = substr($sequence, $end+1, $upstreamRegion);#MODIFIED!! Seems to work.

               } else   #Else take the MAX
               {
                	$promoterSequence = substr($sequence, $end+1, $MAX_PROMOTER); #MODIFIED. Seems to work.
               }

        } else #Case 2. "Previous" gene is in the plus strand. So the two genes are divergent.
        {
                if($upstreamRegion < (2*$MIN_PROMOTER + $SET_OFF))#Skip all regions smaller than 2*MIN.
               {
               } elsif ($upstreamRegion < (2*$MAX_PROMOTER + $SET_OFF))#Take promoters between 2* (MIN and MAX).
               {
               		$upstreamRegion = int ($upstreamRegion/2); #Divide the available promoter region into 2. Take only integer numbers as result.
                        $promoterSequence = substr($sequence, $end+1, $upstreamRegion); #modified. Seems to work

               } else   #Else take the MAX.
               {
                	$promoterSequence = substr($sequence, $end+1, $MAX_PROMOTER); #MODIFIED. Seems to work.
               }

        }
        return ($promoterSequence);
}


#Subroutine to get upstream CDS sequences
sub getUpCDS($)
{
        my $prevGeneData = shift;
        my ($prevPID, $prevStart, $prevEnd, $prevStrand, $prevCOG, $prevProduct) = @{$prevGeneData};
        my $upSequence;
        my $upRegion = $prevEnd - $prevStart;

        #Take upstream CDS region between. No minimum cut-off
       	my $endUp = ($upRegion < $MAX_CDS_LENGTH) ? $upRegion + 1 : $MAX_CDS_LENGTH;

        if ($prevStrand eq "+")
        {
        	$prevStart = $prevStart+$SET_OFF;
        	$upSequence = substr($sequence, $prevStart, $endUp);
        } else
        {
        	$prevStart = $prevEnd-$endUp-$SET_OFF+1;
                $upSequence = substr($sequence, $prevStart, $endUp);
                $upSequence = RevComp($upSequence);
        }

        return ($upSequence);

}

#############################################################################################

#Main sub

#Read the FNAFILE all at once
open FNAFILE, "$CHROMOSOME.fna" or die "Could not open fna file $CHROMOSOME.fna\n";
while (<FNAFILE>) {
	chomp $_;
	if(/^>/) {
		$header = $_;
                $header =~ m/(plasmid|chromosome|genome)/i;
                if ($1)
                {
                	$typeDNA = $1;
                } else
                {
                	$typeDNA = "-";
                }
	} else {
		$sequence .= $_;
	}
}
close FNAFILE;

# Parse ptt file. Example:
#4653 proteins
#Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
#914..2173       -       419     295131885       -       ZPR_0001        -       COG0738G        Multiple antibiotic resistance (MarC)-related protein
#2629..4425      +       598     295131886       -       ZPR_0002        -       COG0481M        GTP-binding protein LepA
#4701..5459      +       252     295131887       -       ZPR_0003        -       COG0217S        hypothetical protein
#5571..7448      -       625     295131888       -       ZPR_0004        -       COG0025P        sodium/hydrogen exchanger family protein

open PTTFILE, "$CHROMOSOME.ptt" or die "Could not open ptt file $CHROMOSOME.ptt\n";
my $titleLine = <PTTFILE>;
my $cProteins = <PTTFILE>;
my $headerLine = <PTTFILE>;
my @genes;

while(<PTTFILE>) {
	chomp;
	my @fields = split /\t/, $_;
	if (scalar @fields != 9) {
		die "Unexpected number of fields in this line $_\n";
	}
	my $location = $fields[0];
	my $strand = $fields[1];
	my $PID = $fields[3];
        my $COG = $fields[7];
        my $product = $fields[8];

	my @locationFields = split /\.\./, $location;
	if (scalar @locationFields != 2) {
		die "Unexpected number of fields in this location $location\n";
	}

	my $start = $locationFields[0] - 1;
	my $end = $locationFields[1] - 1;
	my @geneData = ($PID, $start, $end, $strand, $COG, $product);
	push (@genes, \@geneData);
}
close PTTFILE;

my $lengths_file = $RUN_NAME.".PromLengths";
open (LENGTHS, ">>$lengths_file") or die "Could not create file $lengths_file\n";#Append the results for all genomes.

for(my $iGene = 0; $iGene < scalar @genes; $iGene++) {
	my $geneData = $genes[$iGene];
        my ($PID, $start, $end, $strand, $COG, $product) = @{$geneData};
        my $cdsSequence;
        my $promoterSequence;
        my $upSequence;
        my $startCDS;

        #1) Filter genes by MIN size of CDS.
        my $cdsRegion = $end - $start;
        if ($cdsRegion < $MIN_CDS_LENGTH) #Skip genes shorter then MIN size
        {
               	next;
        }

        #Take CDS region between MIN and MAX
       	my $endCDS = ($cdsRegion < $MAX_CDS_LENGTH) ? $cdsRegion + 1 : $MAX_CDS_LENGTH;

        #Depending if the CDS is in the plus or minus strand is the location from which we take the sequence.
        if ($strand eq "+"){
        	$startCDS = $start + $SET_OFF;
        } else {
        	$startCDS = $end  - $endCDS - $SET_OFF + 1;
        }

        $cdsSequence = substr($sequence, $startCDS, $endCDS);

        #2) Filter genes by MIN size of promoter.
        if ($iGene == 0 )  #First gene of a linear DNA molecule. It only has one adjacent gene (either previous or next).
        {

        	if ($strand eq "+")#The gene is in the plus strand. There's no previous gene.
                {
                	$promoterSequence = checkFirstProm ($start-$SET_OFF);  #MODIFIED

                } else  #The gene is in the minus. There's a gene in the 5' direction.
                {
                        my $prevGeneData = $genes[$iGene + 1]; #Get previous gene info
			$promoterSequence = checkMinusProm($end+$SET_OFF, $prevGeneData); #MODIFIED
                        $upSequence = getUpCDS($prevGeneData);
                } #CHECK WHEN YOU DONT HAVE ANYTHING BEFORE OR AFTER

        } elsif ($iGene + 1 == scalar @genes)  #Last gene of a linear DNA molecule. It only has one adjacent gene (either previous or next).
        {
        	if ($strand eq "+") #There's a previous gene (gene in the 5' direction).
                {
                        my $prevGeneData = $genes[$iGene - 1]; #Get previous gene info
                	$promoterSequence = checkPlusProm($start-$SET_OFF, $prevGeneData);
                        $upSequence = getUpCDS($prevGeneData);

                } else  #The gene is in the minus strand. There's no previous gene.
                {
                	$promoterSequence = checkLastProm ($end+$SET_OFF);
                }


        } else #Any middle gene.
        {
                if ($strand eq "+") #Gene is in plus strand.
                {
                	my $prevGeneData = $genes[$iGene - 1]; #Get previous gene info
                	$promoterSequence = checkPlusProm($start-$SET_OFF, $prevGeneData);
                        $upSequence = getUpCDS($prevGeneData);

                } else  #Gene is in minus strand.
                {
                	my $prevGeneData = $genes[$iGene + 1]; #Get previous gene info
			$promoterSequence = checkMinusProm($end+$SET_OFF, $prevGeneData);
                        $upSequence = getUpCDS($prevGeneData);
                }

        }

        #3) Reverse-complement promoter and CDS when gene is in the minus strand.
        if ($strand eq "-")
        {
        	$promoterSequence = RevComp ($promoterSequence);
                $cdsSequence = RevComp ($cdsSequence);
        }

        #4) Print the promoter/gene pairs that pass the filter criteria.
        if ($promoterSequence)
        {
       		print ">".$PID."_prom\n".$promoterSequence."\n";  #Print promoter sequence of PID
                print ">".$PID."_CDS\n".$cdsSequence."\n"; #Print CDS sequence of PID
                print STDERR ">".$PID."_prom\t$genome\t$chrom\t$typeDNA\t$COG\t$product\n".$promoterSequence."\n";
                print STDERR ">".$PID."_CDS\t$genome\t$chrom\t$typeDNA\t$COG\t$product\n".$cdsSequence."\n";



	        #5)Print upstream CDS sequence if it exists

	        if ($upSequence)
	        {
	                print ">".$PID."_upCDS\n".$upSequence."\n"; #Print upstream CDS sequence of PID
	                print STDERR ">".$PID."_upCDS\t$genome\t$chrom\t$typeDNA\t$COG\t$product\n".$upSequence."\n";
	        }
        }
}

close(LENGTHS);