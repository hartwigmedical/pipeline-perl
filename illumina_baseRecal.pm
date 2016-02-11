#!/usr/bin/perl -w

####################################
### illumina_baseRecal.pm
### - Runs gatk base recalibration
###
### Author: R.F.Ernst
####################################

package illumina_baseRecal;

use strict;
use POSIX qw(tmpnam);

sub runBaseRecalibration {
    ###
    # Run base recalibration
    ###
    my $configuration = shift;
    my %opt = %{$configuration};

    print "Running base recalibration for the following BAM-files:\n";
    
    foreach my $sample (@{$opt{SAMPLES}}){
	### Set in and output bam file
	my $inBam = $opt{BAM_FILES}->{$sample};
	(my $inFlagstat = $inBam) =~ s/bam/flagstat/;
	(my $outBam = $inBam) =~ s/bam/recalibrated.bam/;
	(my $outBai = $inBam) =~ s/bam/recalibrated.bai/;
	(my $outFlagstat = $inBam) =~ s/bam/recalibrated.flagstat/;
	
	print "\t$opt{OUTPUT_DIR}/$sample/mapping/$inBam\n";
	$opt{BAM_FILES}->{$sample} = $outBam;
	
	### Skip base recalibration if done file present
	if (-e "$opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration_$sample.done"){
	    print "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration_$sample.done exists, skipping \n";
	    next;
	}
	
	### Build Queue command
	my $command = "java -Xmx".$opt{BASERECALIBRATION_MASTERMEM}."G -jar $opt{QUEUE_PATH}/Queue.jar ";
	# cluster options
	$command .= "-jobQueue $opt{BASERECALIBRATION_QUEUE} -jobNative \"-pe threaded $opt{BASERECALIBRATION_THREADS} -P $opt{CLUSTER_PROJECT}\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration.jobReport.txt "; #Queue options
	# baseRecalibration options
	$command .= "-S $opt{BASERECALIBRATION_SCALA} -R $opt{GENOME} -I $opt{OUTPUT_DIR}/$sample/mapping/$inBam -mem $opt{BASERECALIBRATION_MEM} -nct $opt{BASERECALIBRATION_THREADS} -nsc $opt{BASERECALIBRATION_SCATTER} ";
	
	### Parsing known files and add them to $command.
	my @knownFiles;
	if($opt{BASERECALIBRATION_KNOWN}) {
	    @knownFiles = split('\t', $opt{BASERECALIBRATION_KNOWN});
	    foreach my $knownFile (@knownFiles) {
		if(! -e $knownFile){ die"ERROR: $knownFile does not exist\n" }
		else { $command .= "-knownSites $knownFile " }
	    }
	}
	### retry option
	if($opt{QUEUE_RETRY} eq 'yes'){
	    $command  .= "-retry 1 ";
	}
	$command .= "-run";

	### Create bash script
	my $jobID = "BR_".$sample."_".get_job_id();
	my $bashFile = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobID.".sh";
	my $logDir = $opt{OUTPUT_DIR}."/".$sample."/logs";

	open BASERECAL_SH, ">$bashFile" or die "cannot open file $bashFile \n";
	print BASERECAL_SH "#!/bin/bash\n\n";
	print BASERECAL_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
	print BASERECAL_SH "cd $opt{OUTPUT_DIR}/$sample/tmp/\n";
	print BASERECAL_SH "echo \"Start base recalibration\t\" `date` \"\t$inBam\t\" `uname -n` >> ../logs/$sample.log\n\n";
	
	print BASERECAL_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$inBam ]\n";
	print BASERECAL_SH "then\n";
	print BASERECAL_SH "\t$command\n";
	print BASERECAL_SH "else\n";
	print BASERECAL_SH "\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$inBam does not exist.\" >&2\n";
	print BASERECAL_SH "fi\n";
	close BASERECAL_SH;
	
	### Submit baserecal bash script
	if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
	    system "qsub -q $opt{BASERECALIBRATION_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{BASERECALIBRATION_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/BaseRecalibration_$sample.out -e $logDir/BaseRecalibration_$sample.err -N $jobID -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." $bashFile";
	} else {
	    system "qsub -q $opt{BASERECALIBRATION_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{BASERECALIBRATION_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/BaseRecalibration_$sample.out -e $logDir/BaseRecalibration_$sample.err -N $jobID $bashFile";
	}
	
	### Create flagstat bash script
	my $jobIDFS = "BRFS_".$sample."_".get_job_id();
	my $bashFileFS = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobIDFS.".sh";
	open BASERECALFS_SH, ">$bashFileFS" or die "cannot open file $bashFileFS \n";
	### Generate FlagStats if gatk .done file present
	print BASERECALFS_SH "cd $opt{OUTPUT_DIR}/$sample/tmp/\n";
	print BASERECALFS_SH "if [ -f $opt{OUTPUT_DIR}/$sample/tmp/.$outBam.done ]\n";
	print BASERECALFS_SH "then\n";
	print BASERECALFS_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{FLAGSTAT_THREADS} $opt{OUTPUT_DIR}/$sample/tmp/$outBam > $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat\n";
	print BASERECALFS_SH "fi\n\n";

	### Check FlagStats and move files if correct else print error
	print BASERECALFS_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat ]\n";
	print BASERECALFS_SH "then\n";
	print BASERECALFS_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BASERECALFS_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BASERECALFS_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	print BASERECALFS_SH "\tthen\n";
	print BASERECALFS_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/$outBam $opt{OUTPUT_DIR}/$sample/mapping/\n";
	print BASERECALFS_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/$outBai $opt{OUTPUT_DIR}/$sample/mapping/\n";
	print BASERECALFS_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/*_baseRecalibration.pdf $opt{OUTPUT_DIR}/$sample/logs/\n";
	print BASERECALFS_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/*_post_recal_data.table $opt{OUTPUT_DIR}/$sample/logs/\n";
	print BASERECALFS_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/*_recal_data.table $opt{OUTPUT_DIR}/$sample/logs/\n";

	print BASERECALFS_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration_$sample.done\n";
	print BASERECALFS_SH "\telse\n";
	print BASERECALFS_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat and $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat do not have the same read counts\" >>../logs/BaseRecalibration_$sample.err\n";
	print BASERECALFS_SH "\tfi\n";
	print BASERECALFS_SH "else\n";
	print BASERECALFS_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat or $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat is empty.\" >> ../logs/BaseRecalibration_$sample.err\n";
	print BASERECALFS_SH "fi\n\n";
	print BASERECALFS_SH "echo \"End base recalibration\t\" `date` \"\t$inBam\t\" `uname -n` >> ../logs/$sample.log\n";
	close BASERECALFS_SH;
	
	### Submit flagstat bash script
	system "qsub -q $opt{FLAGSTAT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{FLAGSTAT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/BaseRecalibrationFS_$sample.out -e $logDir/BaseRecalibrationFS_$sample.err -N $jobIDFS -hold_jid $jobID $bashFileFS";

	push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobID);
	push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobIDFS);
    }
    return \%opt;
}

############
sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}
############ 

1;