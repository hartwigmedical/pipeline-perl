#!/usr/bin/perl -w

package illumina_finalize;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin";
use FindBin;
use illumina_sge;

sub runFinalize {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my $doneFile;
    my @runningJobs;

    my $jobID = $runName."_".get_job_id();
    my $bashFile = "$opt{OUTPUT_DIR}/jobs/Finalize_".$jobID.".sh";
    open (BASH,">$bashFile") or die "ERROR: Couldn't create $bashFile\n";
    print BASH "\#!/bin/sh\n . $opt{CLUSTER_PATH}/settings.sh\n\n";

    my $logFile = "$opt{OUTPUT_DIR}/logs/PipelineCheck.log";
    print BASH "failed=false \n";
    print BASH "rm $logFile \n";
    print BASH "echo \"Check and cleanup for run: $runName \" >>$logFile\n";

    my $version = `git --git-dir $opt{PIPELINE_PATH}/.git describe --tags`;
    print BASH "echo \"Pipeline version: $version \" >>$logFile\n\n";
    print BASH "echo \"\">>$logFile\n\n";

    foreach my $sample (@{$opt{SAMPLES}}){
	    print BASH "echo \"Sample: $sample\" >>$logFile\n";
	    if($opt{PRESTATS} eq "yes"){
			$doneFile = $opt{OUTPUT_DIR}."/$sample/logs/PreStats_$sample.done";
			print BASH "if [ -f $doneFile ]; then\n";
			print BASH "\techo \"\t PreStats: done \" >>$logFile\n";
			print BASH "else\n";
			print BASH "\techo \"\t PreStats: failed \">>$logFile\n";
			print BASH "\tfailed=true\n";
			print BASH "fi\n";
	    }

	    if($opt{MAPPING} eq "yes"){
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

	    if($opt{BAF} eq "yes"){
			$doneFile = $opt{OUTPUT_DIR}."/$sample/logs/BAF_$sample.done";
			print BASH "if [ -f $doneFile ]; then\n";
			print BASH "\techo \"\t BAF analysis: done \" >>$logFile\n";
			print BASH "else\n";
			print BASH "\techo \"\t BAF analysis: failed \">>$logFile\n";
			print BASH "\tfailed=true\n";
			print BASH "fi\n";
			if ( $opt{RUNNING_JOBS}->{'baf'} ){
				push( @runningJobs, @{$opt{RUNNING_JOBS}->{'baf'}} );
			}
	    }
	    print BASH "echo \"\">>$logFile\n\n";

		if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
			push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
		}
    }

    if($opt{POSTSTATS} eq "yes"){
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

    if($opt{VARIANT_CALLING} eq "yes"){
		$doneFile = $opt{OUTPUT_DIR}."/logs/GermlineCaller.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"Germline caller: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"Germline caller: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
    }

	if ($opt{FILTER_VARIANTS} eq "yes") {
		$doneFile = $opt{OUTPUT_DIR}."/logs/GermlineFilter.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"Germline filter: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"Germline filter: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
	}

	if ($opt{ANNOTATE_VARIANTS} eq "yes") {
		$doneFile = $opt{OUTPUT_DIR}."/logs/GermlineAnnotation.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"Germline annotation: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"Germline annotation: failed \">>$logFile\n";
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
		if ($opt{CNV_MODE} eq "sample_control") {
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
		} elsif ($opt{CNV_MODE} eq "sample") {
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

    if ($opt{KINSHIP} eq "yes") {
		$doneFile = $opt{OUTPUT_DIR}."/logs/Kinship.done";
		print BASH "if [ -f $doneFile ]; then\n";
		print BASH "\techo \"Kinship: done \" >>$logFile\n";
		print BASH "else\n";
		print BASH "\techo \"Kinship: failed \">>$logFile\n";
		print BASH "\tfailed=true\n";
		print BASH "fi\n";
		if ( $opt{RUNNING_JOBS}->{'Kinship'} ){
			push( @runningJobs, $opt{RUNNING_JOBS}->{'Kinship'} );
		}
    }

    print BASH "echo \"\">>$logFile\n\n";

	print BASH "if [ \"\$failed\" = true  ]\n";
    print BASH "then\n";
    print BASH "\techo \"One or multiple step(s) of the pipeline failed. \" >>$logFile\n";
    print BASH "\tmail -s \"Pipeline FAILED $runName\" \"$opt{MAIL}\" < $logFile\n";

    print BASH "else\n";
    print BASH "\techo \"The pipeline completed successfully.\">>$logFile\n";

    print BASH "\trm -r $opt{OUTPUT_DIR}/tmp\n";
    print BASH "\trm -r $opt{OUTPUT_DIR}/*/tmp\n";
    print BASH "\tfind $opt{OUTPUT_DIR}/logs -size 0 -not -name \"*.done\" -delete\n";
    print BASH "\tfind $opt{OUTPUT_DIR}/*/logs -size 0 -not -name \"*.done\" -delete\n";
    print BASH "\tfind $opt{OUTPUT_DIR}/somaticVariants/*/logs -size 0 -not -name \"*.done\" -delete\n";

	if($opt{INDELREALIGNMENT} eq "yes"){
		foreach my $sample (@{$opt{SAMPLES}}) {
			print BASH "\trm $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.ba*\n";
		}
    }

    if($opt{SOMATIC_VARIANTS} eq "yes" && $opt{SOMVAR_VARSCAN} eq "yes" && $opt{FINALIZE_KEEP_PILEUP} eq "no"){
		foreach my $sample (@{$opt{SAMPLES}}){
			print BASH "\trm $opt{OUTPUT_DIR}/$sample/mapping/$sample*.pileup.gz\n";
			print BASH "\trm $opt{OUTPUT_DIR}/$sample/mapping/$sample*.pileup.gz.tbi\n";
		}
    }

    $doneFile = "$opt{OUTPUT_DIR}/logs/PipelineCheck.done";
    print BASH "\tmail -s \"Pipeline DONE $runName\" \"$opt{MAIL}\" < $logFile\n";
    print BASH "\ttouch $doneFile\n";

    print BASH "fi\n";

    #Sleep to ensure that email is send from cluster.
    print BASH "sleep 5s \n";

    my $qsub = &qsubTemplate(\%opt, "FINALIZE");
    if (@runningJobs){
		system "$qsub -o /dev/null -e /dev/null -N Finalize_$jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
		system "$qsub -o /dev/null -e /dev/null -N Finalize_$jobID $bashFile";
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
