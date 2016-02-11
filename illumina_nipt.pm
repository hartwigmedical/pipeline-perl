#!/usr/bin/perl -w

#######################################################
### illumina_nipt.pm
### - Run NIPT statistics
###
### Authors: R.F.Ernst
###
#######################################################

package illumina_nipt;

use strict;
use POSIX qw(tmpnam);
use FindBin;

sub runNipt {
    ###
    # Run NIPT statistics
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @runningJobs; #internal job array
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my $jobID = "NIPT_".get_job_id();
    my $jobIDCheck = "NIPT_Check_".get_job_id();

    ## Get running sample jobs
    foreach my $sample (@{$opt{SAMPLES}}){
	if (@{$opt{RUNNING_JOBS}->{$sample}}) {
	    push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	}
    }

    ### Run Chromate
    if(! -e "$opt{OUTPUT_DIR}/logs/NIPT.done"){
	my $command = "python $opt{CHROMATE_PATH} -f ";
	$command .= "-d $opt{NIPT_REFERENCESET} ";
	$command .= "-x $opt{OUTPUT_DIR}/ ";
	
	my $bashFile = $opt{OUTPUT_DIR}."/jobs/".$jobID.".sh";
	my $logDir = $opt{OUTPUT_DIR}."/logs";
        
	open NIPT_SH, ">$bashFile" or die "cannot open file $bashFile\n";
	print NIPT_SH "#!/bin/bash\n\n";
	print NIPT_SH "cd $opt{OUTPUT_DIR}\n";
	print NIPT_SH "echo \"Start NIPT\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
	print NIPT_SH "$command\n";
	close NIPT_SH;
	
	if (@runningJobs){
	    system "qsub -q $opt{NIPT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{NIPT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/NIPT_$runName.out -e $logDir/NIPT_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
	} else {
	    system "qsub -q $opt{NIPT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{NIPT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/NIPT_$runName.out -e $logDir/NIPT_$runName.err -N $jobID $bashFile";
	}
	
	### Check Chromate result
	my $bashFileCheck = $opt{OUTPUT_DIR}."/jobs/".$jobIDCheck.".sh";
	open NIPTCHECK_SH, ">$bashFileCheck" or die "cannot open file $bashFileCheck\n";
	print NIPTCHECK_SH "cd $opt{OUTPUT_DIR}\n";
	print NIPTCHECK_SH "if [ -s $runName.pdf ]\nthen\n";
	print NIPTCHECK_SH "\ttouch logs/NIPT.done \n";
	print NIPTCHECK_SH "fi\n";
	print NIPTCHECK_SH "echo \"Finished NIPT\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
	close NIPTCHECK_SH;

	system "qsub -q $opt{NIPT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{NIPT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/NIPT_$runName.out -e $logDir/NIPT_$runName.err -N $jobIDCheck -hold_jid bamMetrics_report_".$runName.",$jobID $bashFileCheck";
	return $jobIDCheck;

    } else {
	print "WARNING: $opt{OUTPUT_DIR}/logs/NIPT.done exists, skipping\n";
    }
}

############
sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}

############ 

1;