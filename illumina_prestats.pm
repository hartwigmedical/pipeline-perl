#!/usr/bin/perl -w

###################################
### illumina_prestats.pm
### - Run fastqc for each fastq
###
### Author: S.W.Boymans
###################################

package illumina_prestats;

use strict;
use POSIX qw(tmpnam);

sub runPreStats {
    ###
    # Run fastqc
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $jobIds = {};
    
    my $mainJobID = "$opt{OUTPUT_DIR}/jobs/PreStatsMainJob_".get_job_id().".sh";

    open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
    print "Creating FASTQC report for the following fastq.gz files:\n";

    foreach my $input (keys %{$opt{FASTQ}}){
	my $coreName = undef;
	$coreName = (split("/", $input))[-1];
	$coreName =~ s/\.fastq.gz//;
	my ($sampleName) =  split("_", $coreName);
	print "\t$input\n"; #print fastq filename

	if(! -e "$opt{OUTPUT_DIR}/$sampleName/logs/PreStats_$sampleName.done"){

	    my $preStatsJobId = "PreStat_$coreName\_".get_job_id();
	    push(@{$jobIds->{$sampleName}}, $preStatsJobId);
	    open PS,">$opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh";
	    print PS "\#!/bin/sh\n\n";
	    print PS "cd $opt{OUTPUT_DIR}/$sampleName\n\n";
	    print PS "echo \"Start PreStats\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    print PS "$opt{FASTQC_PATH}/fastqc $input -o QCStats --noextract -t $opt{PRESTATS_THREADS}\n";
	    print PS "touch logs/PreStats_$sampleName.done\n";
	    print PS "echo \"End PreStats\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    close PS;

	    print QSUB "qsub -pe threaded $opt{PRESTATS_THREADS} -m a -M $opt{MAIL} -q $opt{PRESTATS_QUEUE} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $opt{OUTPUT_DIR}/$sampleName/logs/PreStat_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/PreStats_$coreName.err -N $preStatsJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh\n";
	} else {
	    print "\t WARNING: FASTQC report for $input already exists, skipping.\n";
	}

    }

    close QSUB;

    system("sh $mainJobID");
}

############
sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}
############
1;