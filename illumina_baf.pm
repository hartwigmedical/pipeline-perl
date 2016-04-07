#!/usr/bin/perl -w

##################################################################
### illumina_baf.pm
### - Run gatk Unified Genotyper and make baf plots.
###
### Authors: R.F.Ernst
##################################################################

package illumina_baf;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_sge;

sub runBAF {
    ###
    # Run BAF analysis per sample
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @baf_jobs;

    foreach my $sample (@{$opt{SAMPLES}}){
	###
	# Setup sample variables
	###
	my $sample_bam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
	my $log_dir = $opt{OUTPUT_DIR}."/".$sample."/logs/";
	my $tmp_dir = $opt{OUTPUT_DIR}."/".$sample."/tmp/";
	my $job_dir = $opt{OUTPUT_DIR}."/".$sample."/jobs/";
	my $output_dir = $opt{OUTPUT_DIR}."/".$sample."/";
	my $command;
	my @running_jobs;

	if (-e "$log_dir/BAF_$sample.done"){
	    print "WARNING: $log_dir/BAF_$sample.done exists, skipping BAF analysis for $sample \n";
	} else {
	    ## Setup baf sh script
	    my $jobID = "BAF_$sample\_".get_job_id();
	    my $bashFile = $job_dir.$jobID.".sh";
	    my $output_vcf = $sample."_BAF_SNPS.vcf";
	    my $output_baf = $sample."_BAF.txt";
	    my $output_bafplot = $sample."_BAF.pdf";

	    ## Running jobs
	    if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
		push( @running_jobs, @{$opt{RUNNING_JOBS}->{$sample}} );
	    }

	    ###
	    # Run Unified Genotyper
	    ###
	    ### Skip if .done file exist
	    my $ug_ok = 0;
	    if (-e "$log_dir/BAF_UG_$sample.done"){
		print "WARNING: $log_dir/BAF_UG_$sample.done exists, skipping Unified Genotyper for $sample \n";
	    } else {
		$ug_ok = 1;
	    }

	    ###
	    # Make BAF file
	    ###
	    ### Skip if .done file exist
	    my $baf_file = 0;
	    if (-e "$log_dir/BAF_FILE_$sample.done"){
		print "WARNING: $log_dir/BAF_FILE_$sample.done exists, skipping BAF file for $sample \n";
	    } else {
		$baf_file = 1;
	    }
	    ###
	    # Run BAF plots
	    ###
	    my $baf_plots = 0;
	    if (-e "$log_dir/BAF_PLOT_$sample.done"){
		print "WARNING: $log_dir/BAF_PLOT._$sample done exists, skipping BAF plot for $sample \n";
	    } else {
		$baf_plots = 1;
	    }

            from_template("BAF_Job.sh.tt", $bashFile, tmp_dir => $tmp_dir, ug_ok => $ug_ok, log_dir => $log_dir, sample => $sample,
		sample_bam => $sample_bam, output_vcf => $output_vcf, output_dir => $output_dir, baf_file => $baf_file, output_baf => $output_baf,
		output_bafplot => $output_bafplot, baf_plots => $baf_plots, opt => \%opt);

	    ###
	    # Submit BAF JOB
	    ###
	    my $qsub = &qsubJava(\%opt,"BAF");
	    if (@running_jobs){
		system "$qsub -o $log_dir/BAF_$sample.out -e $log_dir/BAF_$sample.err -N $jobID -hold_jid ".join(",",@running_jobs)." $bashFile";
	    } else {
		system "$qsub -o $log_dir/BAF_$sample.out -e $log_dir/BAF_$sample.err -N $jobID $bashFile";
	    }
	    push(@baf_jobs, $jobID);
	}
    }
    return \@baf_jobs;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
