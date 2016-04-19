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

	    open BAF_SH, ">$bashFile" or die "cannot open file $bashFile \n";
	    print BAF_SH "#!/bin/bash\n\n";
	    print BAF_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
	    print BAF_SH "cd $tmp_dir\n\n";

	    ## Running jobs
	    if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
		push( @running_jobs, @{$opt{RUNNING_JOBS}->{$sample}} );
	    }

	    ###
	    # Run Unified Genotyper
	    ###
	    ### Skip if .done file exist
	    if (-e "$log_dir/BAF_UG_$sample.done"){
		print "WARNING: $log_dir/BAF_UG_$sample.done exists, skipping Unified Genotyper for $sample \n";
	    } else {
		### Build gatk command
		$command = "java -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp/ -Xmx".$opt{BAF_MEM}."G -jar $opt{QUEUE_PATH}/GenomeAnalysisTK.jar ";
		$command .= "-T UnifiedGenotyper ";
		$command .= "-R $opt{GENOME} ";
		$command .= "-L $opt{BAF_SNPS} ";
		$command .= "-I $sample_bam ";
		$command .= "-o $output_vcf ";
		$command .= "--output_mode EMIT_ALL_SITES ";
		## AD NT / NCT settings?

		#Create UG bash script
		print BAF_SH "echo \"Start Unified Genotyper\t\" `date` \"\t\" `uname -n` >> $log_dir/BAF_$sample.log\n";

		print BAF_SH "if [ -s $sample_bam ]\n";
		print BAF_SH "then\n";
		print BAF_SH "\t$command\n";
		print BAF_SH "else\n";
		print BAF_SH "\techo \"ERROR: Sample bam file do not exist.\" >&2\n";
		print BAF_SH "fi\n\n";

		print BAF_SH "if [ \"\$(tail -n 1 $output_vcf | cut -f 1,2)\" = \"\$(tail -n 1 $opt{BAF_SNPS} | cut -f 1,3)\" ]\n";
		print BAF_SH "then\n";
		print BAF_SH "\tmv $output_vcf $output_dir\n";
		print BAF_SH "\tmv $output_vcf.idx $output_dir\n";
		print BAF_SH "\ttouch $log_dir/BAF_UG_$sample.done\n";
		print BAF_SH "fi\n";
		print BAF_SH "echo \"Finished Unified Genotyper\t\" `date` \"\t\" `uname -n` >> $log_dir/BAF_$sample.log\n\n";
	    }

	    ###
	    # Make BAF file
	    ###
	    ### Skip if .done file exist
	    if (-e "$log_dir/BAF_FILE_$sample.done"){
		print "WARNING: $log_dir/BAF_FILE_$sample.done exists, skipping BAF file for $sample \n";
	    } else {
		$command = "cat $output_dir/$output_vcf | ";
		$command .= "$opt{BIOVCF_PATH}/bio-vcf --num-threads $opt{BAF_THREADS} -i ";
		$command .= "--sfilter '!s.empty? and s.dp>=20' ";
		$command .= "--eval '[r.chrom,r.pos,r.ref+\">\"+r.alt[0]]' ";
		$command .= "--seval 'tot=s.ad.reduce(:+) ; ((tot-s.ad[0].to_f)/tot).round(2)' ";
		$command .= "> $output_baf ";

		print BAF_SH "echo \"Start Make BAF file\t\" `date` \"\t\" `uname -n` >> $log_dir/BAF_$sample.log\n";
		print BAF_SH "if [ -s $output_dir/$output_vcf -a -e $log_dir/BAF_UG_$sample.done ]\n";
		print BAF_SH "then\n";
		print BAF_SH "\t$command\n";
		print BAF_SH "else\n";
		print BAF_SH "\techo \"ERROR: Sample BAF vcf and UG done file do not exist.\" >&2\n";
		print BAF_SH "fi\n\n";

		print BAF_SH "if [ -s $output_baf ]\n";
		print BAF_SH "then\n";
		print BAF_SH "\tmv $output_baf $output_dir\n";
		print BAF_SH "\ttouch $log_dir/BAF_FILE_$sample.done\n";
		print BAF_SH "fi\n";
		print BAF_SH "echo \"Finished Make BAF file\t\" `date` \"\t\" `uname -n` >> $log_dir/BAF_$sample.log\n\n";
	    }
	    ###
	    # Run BAF plots
	    ###
	    if (-e "$log_dir/BAF_PLOT_$sample.done"){
		print "WARNING: $log_dir/BAF_PLOT._$sample done exists, skipping BAF plot for $sample \n";
	    } else {
		$command = "Rscript $opt{BAF_PLOTSCRIPT} $tmp_dir $output_dir/$output_baf ";
		print BAF_SH "echo \"Start BAF plotting\t\" `date` \"\t\" `uname -n` >> $log_dir/BAF_$sample.log\n";
		print BAF_SH "if [ -s $output_dir/$output_baf -a -e $log_dir/BAF_FILE_$sample.done ]\n";
		print BAF_SH "then\n";
		print BAF_SH "\t$command\n";
		print BAF_SH "else\n";
		print BAF_SH "\techo \"ERROR: Sample BAF file and baf file done file do not exist.\" >&2\n";
		print BAF_SH "fi\n\n";

		print BAF_SH "if [ -s $output_bafplot ]\n";;
		print BAF_SH "then\n";
		print BAF_SH "\tmv $output_bafplot $output_dir\n";
		print BAF_SH "\ttouch $log_dir/BAF_PLOT_$sample.done\n";
		print BAF_SH "fi\n";
		print BAF_SH "echo \"Finished Make BAF plot\t\" `date` \"\t\" `uname -n` >> $log_dir/BAF_$sample.log\n\n";
	    }

	    ###
	    # Check all output files
	    ###
	    print BAF_SH "if [ -e $log_dir/BAF_UG_$sample.done -a -e $log_dir/BAF_FILE_$sample.done -a -e $log_dir/BAF_PLOT_$sample.done ]\n";
	    print BAF_SH "then\n";
	    print BAF_SH "\ttouch $log_dir/BAF_$sample.done\n";
	    print BAF_SH "fi\n";

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
