#!/usr/bin/perl -w

use strict;

#Check number of arguments.
if(scalar @ARGV != 1) {
    die "$0 <parameters file>\n";
}

#INPUT:
#1) Direct input is file called parameters.txt
my $paramFile = $ARGV[0];
#......................

my $RUN_PATH = ($1) ? $1 : ".";
my %parameters = ();
open (PARAMETERS, $paramFile) or die "File $paramFile cannot be opened.\n";
while (<PARAMETERS>)
{
	my $line = $_;
	$line =~ s/\s+//g; #Remove spaces and new line characters
        if ($line =~ m/^(\\\\)/) #Ignore comments
        {
        } elsif ($line =~ m/^\w/) #Avoid taking blank lines
        {
        	$line =~ m/(\w+)=(.+)/;
                $parameters{$1}=$2;
        }
}
close(PARAMETERS);

#Remove file results that have the same run name. This prevents problems because we are concatenating results.

	if (-f $parameters{run_name}.".Database")
        {
        	print "Removing previous files $parameters{run_name}\n";
        	`rm $parameters{run_name}.*`;   #Add condition. Only do it if the files exist
        }

#Parse through all genomes to obtain sequences and do BLAST

my $genome;
my @dirstack;
my $dirName;
my $dirHandle;
push @dirstack, $parameters{dir_genomes};

while ($dirName = pop @dirstack) {
        opendir($dirHandle, $dirName) or die "cannot open $dirName\n";

                while($genome = readdir($dirHandle)) {
                        if ($genome !~ /^\.\.?$/ )  {
                                my $path = "$dirName/$genome"; #$genome has the name of each species(GENOME)
                                print "Processing $path ...\n";

                                #Even if BLAST is not run, a database with all promoter sequences must be formed for future steps (.Database file)
                                my $all_chrom = `find $path -name "*.ptt"`; #Store the names of all DNA molecules of a genome (store pathway)
                                my @all_chrom = split(/\n/, $all_chrom);

                                foreach my $chrom (@all_chrom)  #Obtain sequences of each DNA molecule
                                {
                                $chrom =~ s/.ptt$//;
                                #Store results in two files: .seqs file within each subdirectory and .Database as a single file in the running path
                                `perl $RUN_PATH/getCDSandPromRegions.pl $chrom $parameters{min_prom} $parameters{max_prom} $parameters{min_cds} $parameters{max_cds} $parameters{run_name} $parameters{set_off} > $chrom.seqs 2>> $parameters{run_name}.Database`;  #Retrieve promoter/CDS sequences for each gene

                                }

                                if ($parameters{run_blast} eq "T") #Run blast for each GENOME
                                {
                                        #Concatenate the sequences of all DNA molecules of the same GENOME
                                        `cat $path/*.seqs > $path/$genome.fas`;
                                        #Format the retrieved sequences as a BLAST database
                                        `formatdb -i $path/$genome.fas -p F`;
                                        #Do the BLAST search
                                        `blastall -a $parameters{cpus} -W $parameters{word_size} -p blastn -i $path/$genome.fas -d $path/$genome.fas -F F -e $parameters{max_eval} -o $path/$genome.bls`;
                                }

                                #Obtain clusters of promoters that have min lenght and min identity but CDS that have any hit
                                `perl $RUN_PATH/blastClust.pl $path/$genome.bls $parameters{score} $parameters{min_length} $parameters{min_identity} $parameters{run_name} >> $parameters{run_name}.Clusters 2>> $parameters{run_name}.AllProm`;

                           }
                }


}

#Filter downstream CDS paralogy
print "\nFiltering for paralogy of downstream CDSs...\t";
`perl $RUN_PATH/filterParalogs.pl $parameters{run_name}.Clusters $parameters{max_paralogs} > $parameters{run_name}.filterParalogs`;
print "Done.\n";

#Filter upstream CDS paralogy
print "\nFiltering for paralogy of upstream CDSs...\t";
`perl $RUN_PATH/filterUpstreamParalogs.pl $parameters{run_name}.filterParalogs $parameters{max_upparalogs} > $parameters{run_name}.filterUpParalogs`;
print "Done.\n";

#Obtain a representative from each cluster (oufile.Repr1) and all members (outfile.Members1) from .Database file printed before
print "\nGetting representatives of each cluster of PMPs...\t";
`perl $RUN_PATH/getRep.pl $parameters{run_name}.filterUpParalogs $parameters{run_name}.Database`;
print "Done.\n";

####
#Add some extra information to each cluster  (number of members per cluster and cluster id)
`perl $RUN_PATH/formatDatabase.pl $parameters{run_name}.filterUpParalogs $parameters{run_name}.Repr1 > $parameters{run_name}.Repr2`;
`perl $RUN_PATH/formatDatabase.pl $parameters{run_name}.filterUpParalogs $parameters{run_name}.Members1 > $parameters{run_name}.Members2`;

#Format representatives file for CD-HIT: substitute white spaces for * so that full names appear in the clusters files
`perl $RUN_PATH/formatcdhit.pl $parameters{run_name}.Repr2 > $parameters{run_name}.formatCDHIT`;

#Use CD-HIT to remove redundancy in the representatives.
# c= defines the level of identity (0.0 to 1.0)
# aL = defines the minimum fraction of the representative that must be covered in the aln (0.0 to 1.0)
print "\nDoing inter-genome clustering with CD-HIT...\t";
`$parameters{dir_CDHIT}/./cd-hit-est -i $parameters{run_name}.formatCDHIT -o $parameters{run_name}.CDHIT -d 0 -c $parameters{c} -n 5 -r 0 -g 1 -b 20 -G 0 -aL $parameters{aL} -s 0.0 -aS 0.0`;
print "Done.\n";

#Format the representatives and members file for final visualization (replace * for white spaces)
`perl $RUN_PATH/finalFormat.pl $parameters{run_name}.CDHIT $parameters{run_name}.Repr2 > $parameters{run_name}.ReprFinal`;
`perl $RUN_PATH/finalFormat.pl $parameters{run_name}.CDHIT $parameters{run_name}.Members2 > $parameters{run_name}.MembersFinal`;

#Search the RFAM database of regulatory RNAs. Use .CDHIT file to search the representatives because it has the complete annotation of each cluster.
print "\nSearching the RFAM database for regulatory RNAs...\t";
`mpd &`;
`/usr/bin/mpirun -np 10 /geninf/prog64/bin/cmsearch --mpi --tc $parameters{dir_RFAM}/Rfam.cm $parameters{run_name}.CDHIT > $parameters{run_name}.ReprCMSEARCH`;
`grep 'prom' $parameters{run_name}.ReprCMSEARCH -B 1 -A 5 > $parameters{run_name}.ReprRFAMResults`;
`/usr/bin/mpirun -np 10 /geninf/prog64/bin/cmsearch --mpi --tc $parameters{dir_RFAM}/Rfam.cm $parameters{run_name}.MembersFinal > $parameters{run_name}.MembersCMSEARCH`;
`grep 'prom' $parameters{run_name}.MembersCMSEARCH > $parameters{run_name}.MembersRFAMResults`;
print "Done.\n";

print "\nPIPELINE FINISHED. HAVE A GOOD DAY!\n\n";
#Searth the ISFinder database through web server. Use file $parameters{run_name}.CDHIT