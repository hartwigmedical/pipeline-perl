#!/usr/bin/perl -w

########################################################################
### illumina_check.pm
### - Check the result of the pipeline based on the selected modules.
### - Remove tmp files if pipeline completed successfully.
### Author: R.F.Ernst
########################################################################

package illumina_check;

use strict;
use POSIX qw(tmpnam);
use FindBin;

sub runCheck {
    ### 
    # Run checks and email result
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my $doneFile;
    my @runningJobs;

    ### Create bash file
    my $jobID = get_job_id();
    my $bashFile = "$opt{OUTPUT_DIR}/jobs/check_".get_job_id().".sh";
    open (BASH,">$bashFile") or die "ERROR: Couldn't create $bashFile\n";
    print BASH "\#!/bin/sh\n . $opt{CLUSTER_PATH}/settings.sh\n\n";

    ### Log file
    my $logFile = "$opt{OUTPUT_DIR}/logs/PipelineCheck.log";
    print BASH "failed=false \n";
    print BASH "rm $logFile \n";
    print BASH "echo \"Check and cleanup for run: $runName \" >>$logFile\n";

    ### pipeline version
    my $version = `git --git-dir $FindBin::Bin/.git describe --tags`;
    print BASH "echo \"Pipeline version: $version \" >>$logFile\n\n";
    print BASH "echo \"\">>$logFile\n\n"; ## empty line between samples

    ### Check sample steps
    foreach my $sample (@{$opt{SAMPLES}}){
	if(! $opt{VCF}) {
	    print BASH "echo \"Sample: $sample\" >>$logFile\n";
	    if($opt{PRESTATS} eq "yes" && ! $opt{BAM}){
		$doneFile = $opt{OUTPUT_DIR}."/$sample/logs/PreStats_$sample.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"\t PreStats: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"\t PreStats: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
	    }
	    if($opt{MAPPING} eq "yes" && ! $opt{BAM}){
		$doneFile = $opt{OUTPUT_DIR}."/$sample/logs/Mapping_$sample.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"\t Mapping: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"\t Mapping: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
	    }
	    if($opt{INDELREALIGNMENT} eq "yes"){
		$doneFile = $opt{OUTPUT_DIR}."/$sample/logs/Realignment_$sample.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"\t Indel realignment: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"\t Indel realignment: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
	    }
	    if($opt{BASEQUALITYRECAL} eq "yes"){
		$doneFile = $opt{OUTPUT_DIR}."/$sample/logs/BaseRecalibration_$sample.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"\t Base recalibration: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"\t Base recalibration: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
	    }
	    print BASH "echo \"\">>$logFile\n\n"; ## empty line between samples
	}
	## Running jobs
	if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
	}
    }

    ### Check run steps
    if($opt{POSTSTATS} eq "yes" && ! $opt{VCF}){
	$doneFile = $opt{OUTPUT_DIR}."/logs/PostStats.done";
	print BASH "if [ -f $doneFile ]; then\n";
	print BASH "\techo \"PostStats: done \" >>$logFile\n";
	print BASH "else\n";
	print BASH "\techo \"PostStats: failed \">>$logFile\n";
	print BASH "\tfailed=true\n";
	print BASH "fi\n";
	if ( $opt{RUNNING_JOBS}->{'postStats'} ){
	    push( @runningJobs, $opt{RUNNING_JOBS}->{'postStats'} );
	}
    }
    if($opt{NIPT} eq "yes" && ! $opt{VCF}){
	$doneFile = $opt{OUTPUT_DIR}."/logs/NIPT.done";
	print BASH "if [ -f $doneFile ]; then\n";
	print BASH "\techo \"NIPT: done \" >>$logFile\n";
	print BASH "else\n";
	print BASH "\techo \"NIPT: failed \">>$logFile\n";
	print BASH "\tfailed=true\n";
	print BASH "fi\n";
	if ( $opt{RUNNING_JOBS}->{'nipt'} ){
	    push( @runningJobs, $opt{RUNNING_JOBS}->{'nipt'} );
	}
    }
    if($opt{VARIANT_CALLING} eq "yes" && ! $opt{VCF}){
	$doneFile = $opt{OUTPUT_DIR}."/logs/VariantCaller.done";
	print BASH "if [ -f $doneFile ]; then\n";
	print BASH "\techo \"Variant caller: done \" >>$logFile\n";
	print BASH "else\n";
	print BASH "\techo \"Variant caller: failed \">>$logFile\n";
	print BASH "\tfailed=true\n";
	print BASH "fi\n";
    }
    if($opt{SOMATIC_VARIANTS} eq "yes"){
	print BASH "echo \"Somatic variants:\" >>$logFile\n";
	foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){
	    foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){
		foreach my $sample_ref (@{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}){
		    my $sample_tumor_name = "$sample_ref\_$sample_tumor";
		    my $done_file = "$opt{OUTPUT_DIR}/somaticVariants/$sample_tumor_name/logs/$sample_tumor_name.done";
		    print BASH "if [ -f $done_file ]; then\n";
		    print BASH "\techo \"\t $sample_tumor_name: done \" >>$logFile\n";
		    print BASH "else\n";
		    print BASH "\techo \"\t $sample_tumor_name: failed \">>$logFile\n";
		    print BASH "\tfailed=true\n";
		    print BASH "fi\n";
		}
	    }
	}
	if ( $opt{RUNNING_JOBS}->{'somVar'} ){
	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{'somVar'}} );
	}
    }
    if($opt{COPY_NUMBER} eq "yes"){
	print BASH "echo \"Copy number analysis:\" >>$logFile\n";
	if($opt{CNV_MODE} eq "sample_control"){
	    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){
		foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){
		    foreach my $sample_ref (@{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}){
			my $sample_tumor_name = "$sample_ref\_$sample_tumor";
			my $done_file = "$opt{OUTPUT_DIR}/copyNumber/$sample_tumor_name/logs/$sample_tumor_name.done";
			print BASH "if [ -f $done_file ]; then\n";
			print BASH "\techo \"\t $sample_tumor_name: done \" >>$logFile\n";
			print BASH "else\n";
			print BASH "\techo \"\t $sample_tumor_name: failed \">>$logFile\n";
			print BASH "\tfailed=true\n";
			print BASH "fi\n";
		    }
		}
	    }
	} elsif($opt{CNV_MODE} eq "sample"){
	    foreach my $sample (@{$opt{SAMPLES}}){
		my $done_file = "$opt{OUTPUT_DIR}/copyNumber/$sample/logs/$sample.done";
		print BASH "if [ -f $done_file ]; then\n";
		print BASH "\techo \"\t $sample: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"\t $sample: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
	    }
	}
	if ( $opt{RUNNING_JOBS}->{'CNV'} ){
	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{'CNV'}} );
	}
    }
    if($opt{SV_CALLING} eq "yes"){
	print BASH "echo \"SV calling:\" >>$logFile\n";
	# per sv type done file check
	my @svTypes = split/\t/, $opt{DELLY_SVTYPE};
	foreach my $type (@svTypes){
	    my $done_file = "$opt{OUTPUT_DIR}/DELLY/logs/DELLY_$type.done"; 
	    print BASH "if [ -f $done_file ]; then\n";
	    print BASH "\techo \"\t $type: done \" >>$logFile\n";
	    print BASH "else\n";
	    print BASH "\techo \"\t $type: failed \">>$logFile\n";
	    print BASH "\tfailed=true\n";
	    print BASH "fi\n";
	}
	if ( $opt{RUNNING_JOBS}->{'sv'} ){
	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{'sv'}} );
	}
    }
    if($opt{FILTER_VARIANTS} eq "yes"){
	$doneFile = $opt{OUTPUT_DIR}."/logs/VariantFilter.done";
	print BASH "if [ -f $doneFile ]; then\n";
	print BASH "\techo \"Variant filter: done \" >>$logFile\n";
	print BASH "else\n";
	print BASH "\techo \"Variant filter: failed \">>$logFile\n";
	print BASH "\tfailed=true\n";
	print BASH "fi\n";
    }
    if($opt{ANNOTATE_VARIANTS} eq "yes"){
	$doneFile = $opt{OUTPUT_DIR}."/logs/VariantAnnotation.done";
	print BASH "if [ -f $doneFile ]; then\n";
	print BASH "\techo \"Variant annotation: done \" >>$logFile\n";
	print BASH "else\n";
	print BASH "\techo \"Variant annotation: failed \">>$logFile\n";
	print BASH "\tfailed=true\n";
	print BASH "fi\n";
    }
    if($opt{VCF_UTILS} eq "yes"){
	$doneFile = $opt{OUTPUT_DIR}."/logs/VCF_UTILS.done";
	print BASH "if [ -f $doneFile ]; then\n";
	print BASH "\techo \"VCF Utils: done \" >>$logFile\n";
	print BASH "else\n";
	print BASH "\techo \"VCF Utils: failed \">>$logFile\n";
	print BASH "\tfailed=true\n";
	print BASH "fi\n";
	if ( $opt{RUNNING_JOBS}->{'VCF_UTILS'} ){
	    push( @runningJobs, $opt{RUNNING_JOBS}->{'VCF_UTILS'} );
	}
    }

    ### Check failed variable and mail report
    print BASH "echo \"\">>$logFile\n\n"; ## empty line after stats

    ### Pipeline failed
    print BASH "if [ \"\$failed\" = true  ]\n";
    print BASH "then\n";
    print BASH "\techo \"One or multiple step(s) of the pipeline failed. \" >>$logFile\n";
    print BASH "\tmail -s \"IAP FAILED $runName\" \"$opt{MAIL}\" < $logFile\n";

    ### Pipeline done
    print BASH "else\n";
    print BASH "\techo \"The pipeline completed successfully. The md5sum file will be created.\">>$logFile\n";
    print BASH "\tmail -s \"IAP DONE $runName\" \"$opt{MAIL}\" < $logFile\n";
    
    # Remove all tmp folders and empty logs except .done files if pipeline completed successfully
    print BASH "\trm -r $opt{OUTPUT_DIR}/tmp\n";
    print BASH "\trm -r $opt{OUTPUT_DIR}/*/tmp\n";
    print BASH "\tfind $opt{OUTPUT_DIR}/logs -size 0 -not -name \"*.done\" -delete\n";
    print BASH "\tfind $opt{OUTPUT_DIR}/*/logs -size 0 -not -name \"*.done\" -delete\n";
    print BASH "\tfind $opt{OUTPUT_DIR}/somaticVariants/*/logs -size 0 -not -name \"*.done\" -delete\n";
    if($opt{INDELREALIGNMENT} eq "yes"){
	foreach my $sample (@{$opt{SAMPLES}}){
	    if($opt{MAPPING_MARKDUP} eq "sample" || $opt{MAPPING_MARKDUP} eq "lane"){
		print BASH "\trm $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.ba*\n";
	    } else {
		print BASH "\trm $opt{OUTPUT_DIR}/$sample/mapping/$sample.ba*\n";
	    }
	}
    }
    # Create md5sum.txt
    print BASH "\n\tcd $opt{OUTPUT_DIR}\n";
    print BASH "\tfind . -type f \\( ! -iname \"md5sum.txt\" \\) -exec md5sum \"{}\" \\; > md5sum.txt\n";

    print BASH "fi\n";
    
    #Sleep to ensure that email is send from cluster.
    print BASH "sleep 5s \n";

    #Start main bash script
    if (@runningJobs){
	system "qsub -q $opt{CHECKING_QUEUE} -m as -M $opt{MAIL} -pe threaded $opt{CHECKING_THREADS} -P $opt{CLUSTER_PROJECT} -o /dev/null -e /dev/null -N check_$jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{CHECKING_QUEUE} -m as -M $opt{MAIL} -pe threaded $opt{CHECKING_THREADS} -P $opt{CLUSTER_PROJECT} -o /dev/null -e /dev/null -N check_$jobID $bashFile";
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
