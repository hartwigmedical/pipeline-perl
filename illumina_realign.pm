#!/usr/bin/perl -w

###################################################
### illumina_realign.pm
### - Realign bam files using gatk indel realigner
###
### Author: S.W.Boymans & R.F.Ernst
####
###################################################

package illumina_realign;

use strict;
use POSIX qw(tmpnam);

sub runRealignment {
    ###
    # Submit indel realign jobs.
    # Two modes: single and multi sample.
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $realignJobs = {};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];

    print "Running $opt{REALIGNMENT_MODE} sample indel realignment for the following BAM-files:\n";
    
    ### Parsing known indel files
    my @knownIndelFiles;
    if($opt{REALIGNMENT_KNOWN}) {
	@knownIndelFiles = split('\t', $opt{REALIGNMENT_KNOWN});
    }
    
    
    ### Multi sample realignment
    if($opt{REALIGNMENT_MODE} eq 'multi'){
	my $mainJobID = "$opt{OUTPUT_DIR}/jobs/RealignMainJob_".get_job_id().".sh";
	open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
	print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";

	my $jobId = "RE_".get_job_id();
	my $cleanupJobId = "REALIGN_CLEANUP\_".get_job_id();
	my $mergeJobs = "";
	my @waitFor = ();
	
	open REALIGN_SH,">$opt{OUTPUT_DIR}/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/jobs/$jobId.sh\n";
	print REALIGN_SH "\#!/bin/sh\n\n";
	print REALIGN_SH ". $opt{CLUSTER_PATH}/settings.sh\n\n";
	print REALIGN_SH "cd $opt{OUTPUT_DIR}/tmp\n";
	print REALIGN_SH "uname -n > ../logs/$jobId.host\n";
	print REALIGN_SH "echo \"Start indel realignment\t\" `date` >> ../logs/$runName.log\n";
	print REALIGN_SH "java -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp/ -Xmx".$opt{REALIGNMENT_MASTERMEM}."G -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -jobQueue $opt{REALIGNMENT_QUEUE} -nt $opt{REALIGNMENT_THREADS} -mem $opt{REALIGNMENT_MEM} -nsc $opt{REALIGNMENT_SCATTER} -mode $opt{REALIGNMENT_MODE} -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS} -P $opt{CLUSTER_PROJECT}\" -run ";
	
	if($opt{REALIGNMENT_KNOWN}) {
	    foreach my $knownIndelFile (@knownIndelFiles) {
		print REALIGN_SH "-known $knownIndelFile ";
	    }
	}
	
	if($opt{QUEUE_RETRY} eq 'yes'){
	    print REALIGN_SH "-retry 1 ";
	}
	
	open CLEAN_SH, ">$opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh\n";
	print CLEAN_SH "\#!/bin/sh\n\n";
	print CLEAN_SH "uname -n > $opt{OUTPUT_DIR}/logs/$cleanupJobId.host\n";
	print CLEAN_SH "PASS=0\n";
	
	foreach my $sample (@{$opt{SAMPLES}}){
	    my $bam = $opt{BAM_FILES}->{$sample};
	    (my $flagstat = $bam) =~ s/\.bam/\.flagstat/;
	    (my $realignedBam = $bam) =~ s/\.bam/\.realigned\.bam/;
	    (my $realignedBai = $bam) =~ s/\.bam/\.realigned\.bai/;
	    (my $realignedFlagstat = $bam) =~ s/\.bam/\.realigned\.flagstat/;
	    $opt{BAM_FILES}->{$sample} = $realignedBam;
	    
	    print "\t$opt{OUTPUT_DIR}/$sample/mapping/$bam\n";
	    
	    ## Check for realigned bam file, skip sample if realigned bam file already exist.
	    if (-e "$opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done"){
		print "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done exists, skipping\n";
		next;
	    }

	    push(@waitFor, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	    print REALIGN_SH "-I $opt{OUTPUT_DIR}/$sample/mapping/$bam";

	    my $mergeJobId = "REALIGN_MERGE_$sample\_".get_job_id();

	    open MERGE_SH, ">$opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
	    print MERGE_SH "\#!/bin/sh\n\n";
    	    print MERGE_SH "cd $opt{OUTPUT_DIR}/tmp/.queue/\n";
	    print MERGE_SH "CHUNKS=`find \$PWD -name '*$realignedBam' | sort | xargs`\n";
	    print MERGE_SH "if [ ! -f $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done ]\n";
	    print MERGE_SH "then\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{REALIGNMENT_MERGETHREADS} $opt{OUTPUT_DIR}/$sample/mapping/$realignedBam \`echo \$CHUNKS\` 1>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.log 2>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{REALIGNMENT_MERGETHREADS} $opt{OUTPUT_DIR}/$sample/mapping/$realignedBam $opt{OUTPUT_DIR}/$sample/mapping/$realignedBai\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{REALIGNMENT_MERGETHREADS} $opt{OUTPUT_DIR}/$sample/mapping/$realignedBam > $opt{OUTPUT_DIR}/$sample/mapping/$$realignedFlagstat\n\n";
	    print MERGE_SH "fi\n\n";
	    print MERGE_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$flagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat ]\n";
	    print MERGE_SH "then\n";
	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print MERGE_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print MERGE_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	    print MERGE_SH "\tthen\n";
	    print MERGE_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done\n";
	    print MERGE_SH "\telse\n";
	    print MERGE_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$flagstat and $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat do not have the same read counts\" >>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "\tfi\n";
	    print MERGE_SH "else\n";
	    print MERGE_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$flagstat or $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat is empty.\" >> $opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "fi\n";
	    close MERGE_SH;
	    
	    print CLEAN_SH "if [ ! -f $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done ]\n";
	    print CLEAN_SH "then\n";
	    print CLEAN_SH "\tPASS=1\n";
	    print CLEAN_SH "else\n";
	    print CLEAN_SH "\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$realignedBam didn't finish properly.\" >> $opt{OUTPUT_DIR}/logs/realn_cleanup.err\n";
	    print CLEAN_SH "fi\n\n";
	    
	    $mergeJobs .= "qsub -q $opt{REALIGNMENT_QUEUE} -p 100 -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MERGETHREADS} -P $opt{CLUSTER_PROJECT} -o $opt{OUTPUT_DIR}/$sample/logs -e $opt{OUTPUT_DIR}/$sample/logs -N $mergeJobId -hold_jid $jobId $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
	    push(@{$opt{RUNNING_JOBS}->{$sample}}, $mergeJobId);
	}
	
	print REALIGN_SH "-jobRunner GridEngine 1>>$opt{OUTPUT_DIR}/logs/$jobId.host 2>>$opt{OUTPUT_DIR}/logs/$jobId.host\n";
	close REALIGN_SH;
	
	print CLEAN_SH "if [ \$PASS -eq 0 ]\n";
	print CLEAN_SH "then\n";
	print CLEAN_SH "echo \"Finished indel realignment\t\" `date` >> ../logs/$runName.log\n";
	print CLEAN_SH "\tmv $opt{OUTPUT_DIR}/tmp/IndelRealigner.jobreport.txt $opt{OUTPUT_DIR}/logs/IndelRealigner.jobreport.txt\n";
	print CLEAN_SH "\tmv $opt{OUTPUT_DIR}/tmp/IndelRealigner.jobreport.pdf $opt{OUTPUT_DIR}/logs/IndelRealigner.jobreport.pdf\n";
	print CLEAN_SH "fi\n";
	close CLEAN_SH;
	
	print QSUB "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $opt{OUTPUT_DIR}/logs -e $opt{OUTPUT_DIR}/logs -N $jobId -hold_jid ".join(",", @waitFor)." $opt{OUTPUT_DIR}/jobs/$jobId.sh\n";
	print QSUB "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $opt{OUTPUT_DIR}/logs -e $opt{OUTPUT_DIR}/logs -N $cleanupJobId -hold_jid $jobId $opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh\n";
	print QSUB $mergeJobs."\n";
	
	system("sh $mainJobID");
    }
    
    ### Single sample indel realignment
    elsif($opt{REALIGNMENT_MODE} eq 'single'){
	foreach my $sample (@{$opt{SAMPLES}}){
	    my $bam = $opt{BAM_FILES}->{$sample};
	    (my $flagstat = $bam) =~ s/\.bam/.flagstat/;
	    (my $realignedBam = $bam) =~ s/\.bam/\.realigned\.bam/;
	    (my $realignedBai = $bam) =~ s/\.bam/\.realigned\.bai/;
	    (my $realignedFlagstat = $bam) =~ s/\.bam/\.realigned\.flagstat/;
	    $opt{BAM_FILES}->{$sample} = $realignedBam;

	    print "\t$opt{OUTPUT_DIR}/$sample/mapping/$bam\n";

	    ## Check for realigned bam file, skip sample if realigned bam file already exist.
	    if (-e "$opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done"){
		print "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done exists, skipping\n";
		next;
	    }
	    
	    ### Create realign bash script
	    my $logDir = $opt{OUTPUT_DIR}."/".$sample."/logs";
	    my $jobID = "Realign_".$sample."_".get_job_id();
	    my $bashFile = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobID.".sh";

	    open REALIGN_SH,">$bashFile" or die "Couldn't create $bashFile\n";

	    print REALIGN_SH "\#!/bin/bash\n\n";
	    print REALIGN_SH ". $opt{CLUSTER_PATH}/settings.sh\n\n";
	    print REALIGN_SH "cd $opt{OUTPUT_DIR}/$sample/tmp \n\n";
	    print REALIGN_SH "echo \"Start indel realignment\t\" `date` \"\t$bam\t\" `uname -n` >> $logDir/$sample.log\n\n";
	    
	    print REALIGN_SH "if [ -f $opt{OUTPUT_DIR}/$sample/mapping/$bam ]\n";
	    print REALIGN_SH "then\n";
	    print REALIGN_SH "\tjava -Djava.io.tmpdir=$opt{OUTPUT_DIR}/$sample/tmp -Xmx".$opt{REALIGNMENT_MASTERMEM}."G -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -jobQueue $opt{REALIGNMENT_QUEUE} -nt $opt{REALIGNMENT_THREADS} -mem $opt{REALIGNMENT_MEM} -nsc $opt{REALIGNMENT_SCATTER} -mode $opt{REALIGNMENT_MODE} -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS} -P $opt{CLUSTER_PROJECT}\" ";
	    
	    if($opt{REALIGNMENT_KNOWN}) {
		foreach my $knownIndelFile (@knownIndelFiles) {
		    if(! -e $knownIndelFile){ die"ERROR: $knownIndelFile does not exist\n" }
		    else { print REALIGN_SH "-known $knownIndelFile " }
		}
	    }

	    if($opt{QUEUE_RETRY} eq 'yes'){
		print REALIGN_SH "-retry 1 ";
	    }

	    print REALIGN_SH "-run -I $opt{OUTPUT_DIR}/$sample/mapping/$bam -jobRunner GridEngine\n";
	    print REALIGN_SH "else\n";
	    print REALIGN_SH "echo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$bam does not exist.\" >&2\n";
	    print REALIGN_SH "fi\n\n";
	    
	    close REALIGN_SH;

	    ### Submit realign bash script
	    if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
		system "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/Realignment_$sample.out -e $logDir/Realignment_$sample.err -N $jobID -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." $bashFile";
	    } else {
		system "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/Realignment_$sample.out -e $logDir/Realignment_$sample.err -N $jobID $bashFile";
	    }
	    
	    ### Create flagstat bash script
	    my $jobIDFS = "RealignFS_".$sample."_".get_job_id();
	    my $bashFileFS = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobIDFS.".sh";

	    open REALIGNFS_SH, ">$bashFileFS" or die "cannot open file $bashFileFS \n";
	    print REALIGNFS_SH "cd $opt{OUTPUT_DIR}/$sample/tmp\n";
	    
	    print REALIGNFS_SH "if [ -s $opt{OUTPUT_DIR}/$sample/tmp/$realignedBam ]\n";
	    print REALIGNFS_SH "then\n";
	    print REALIGNFS_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/$sample/tmp/$realignedBam > $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat\n";
	    print REALIGNFS_SH "\tmv $opt{OUTPUT_DIR}/$sample/tmp/$realignedBam $opt{OUTPUT_DIR}/$sample/mapping/$realignedBam\n";
	    print REALIGNFS_SH "\tmv $opt{OUTPUT_DIR}/$sample/tmp/$realignedBai $opt{OUTPUT_DIR}/$sample/mapping/$realignedBai\n";
	    print REALIGNFS_SH "fi\n\n";
	    
	    print REALIGNFS_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$flagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat ]\n";
	    print REALIGNFS_SH "then\n";
	    print REALIGNFS_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print REALIGNFS_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print REALIGNFS_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	    print REALIGNFS_SH "\tthen\n";
	    print REALIGNFS_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done\n";
	    print REALIGNFS_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/IndelRealigner.jobreport.txt $opt{OUTPUT_DIR}/$sample/logs/IndelRealigner.jobreport.txt\n";
	    print REALIGNFS_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/IndelRealigner.jobreport.pdf $opt{OUTPUT_DIR}/$sample/logs/IndelRealigner.jobreport.pdf\n";
	    print REALIGNFS_SH "\telse\n";
	    print REALIGNFS_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$flagstat and $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat do not have the same read counts\" >>../logs/Realignment_$sample.err\n";
	    print REALIGNFS_SH "\tfi\n";
	    print REALIGNFS_SH "else\n";
	    print REALIGNFS_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$flagstat or $opt{OUTPUT_DIR}/$sample/mapping/$realignedFlagstat is empty.\" >> ../logs/Realignment_$sample.err\n";
	    print REALIGNFS_SH "fi\n\n";
	    
	    print REALIGNFS_SH "echo \"End indel realignment\t\" `date` \"\t$sample\_dedup.bam\t\" `uname -n` >> $logDir/$sample.log\n"; 
	    close REALIGNFS_SH;
	    
	    ### Submit flagstat bash script
	    system "qsub -q $opt{FLAGSTAT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{FLAGSTAT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/RealignmentFS_$sample.out -e $logDir/RealignmentFS_$sample.err -N $jobIDFS -hold_jid $jobID $bashFileFS";
	    
	    push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobID);
	    push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobIDFS);
	}
	
    }else{
	die "ERROR: Invalid REALIGNMENT_MODE $opt{REALIGNMENT_MODE} , use 'single' or 'multi'\n";
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