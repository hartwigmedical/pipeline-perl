#!/usr/bin/perl -w

##################################################################################################################################################
### illumina_copyNumber.pm
### - Run copy number tools
###   - Contra and freec
###   - Tow modes: sample_control (CPCT) & sample (WGS only)
###
### Author: R.F.Ernst
##################################################################################################################################################

package illumina_copyNumber;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);

sub parseSamples {
    ###
    # Parse samples names used in sample_control mode
    # Expects CPCT samples (CPCT........T/R)
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my %somatic_samples;

    foreach my $sample (@{$opt{SAMPLES}}){
	# Parse cpct samples based on expected naming
	my ($cpct_name,$origin) = ($sample =~ /$opt{CNV_REGEX}/);
	if ( (! $cpct_name) || (! $origin) ){
	    print "WARNING: $sample is not passing copy number samplename parsing, skipping \n\n";
	    next;
	}
	
	# Reference sample
	if ($origin =~ m/R.*/){
	    push(@{$somatic_samples{$cpct_name}{"ref"}},$sample);
	}

	# Tumor samples
	elsif ($origin =~ m/T.*/){
	    push(@{$somatic_samples{$cpct_name}{"tumor"}},$sample);
	}
    }

    $opt{SOMATIC_SAMPLES} = {%somatic_samples};
    return \%opt;
}

sub runCopyNumberTools {
    ### 
    # Run copy number tools and check completion
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @check_cnv_jobs;
    
    #### Sample Control mode with control or ref sample
    if($opt{CNV_MODE} eq "sample_control"){
	### Loop over somatic samples
	foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){

	    # Check correct sample ref
	    if (! $opt{SOMATIC_SAMPLES}{$sample}{'ref'}){
		print "WARNING: No ref sample for $sample, skipping \n";
		next;
	    }
	    ### Loop over tumor samples for each somatic sample (multiple tumor samples possible)
	    foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){
		foreach my $sample_ref (@{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}){
		    my @cnv_jobs;
		    ## Create output, log and job directories
		    my $sample_tumor_name = "$sample_ref\_$sample_tumor";
		    my $sample_tumor_out_dir = "$opt{OUTPUT_DIR}/copyNumber/$sample_tumor_name";
		    my $sample_tumor_log_dir = "$sample_tumor_out_dir/logs/";
		    my $sample_tumor_job_dir = "$sample_tumor_out_dir/jobs/";

		    if(! -e $sample_tumor_out_dir){
			make_path($sample_tumor_out_dir) or die "Couldn't create directory:  $sample_tumor_out_dir\n";
		    }
		    if(! -e $sample_tumor_job_dir){
			make_path($sample_tumor_job_dir) or die "Couldn't create directory: $sample_tumor_job_dir\n";
		    }
		    if(! -e $sample_tumor_log_dir){
			make_path($sample_tumor_log_dir) or die "Couldn't create directory: $sample_tumor_log_dir\n";
		    }

		    ## Lookup running jobs and bams
		    my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
		    my @running_jobs;
		    if ( @{$opt{RUNNING_JOBS}->{$sample_tumor}} ){
			push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
		    }
		    my $sample_ref_bam = "$opt{OUTPUT_DIR}/$sample_ref/mapping/$opt{BAM_FILES}->{$sample_ref}";
		    if ( @{$opt{RUNNING_JOBS}->{$sample_ref}} ){
			push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_ref}});
		    }

		    ## Print sample and bam info
		    print "\n$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

		    ## Skip Copy number tools if .done file exist
		    if (-e "$sample_tumor_log_dir/$sample_tumor_name.done"){
			print "WARNING: $sample_tumor_log_dir/$sample_tumor_name.done exists, skipping \n";
			next;
		    }

		    ## Run CNV callers
		    if($opt{CNV_CONTRA} eq "yes"){
			print "\n###SCHEDULING CONTRA####\n";
			my $contra_job = runContra($sample_tumor, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
			if($contra_job){push(@cnv_jobs, $contra_job)};
		
			my $contravis_job = runContraVisualization($sample_tumor, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir,$contra_job, \%opt);
			if($contravis_job){push(@cnv_jobs, $contravis_job)};
		    }
		    if($opt{CNV_FREEC} eq "yes"){
			print "\n###SCHEDULING FREEC####\n";
			my $freec_job = runFreec($sample_tumor, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
			if($freec_job){push(@cnv_jobs, $freec_job)};
		    }
		    ## Check copy number analysis
		    my $job_id = "CHECK_".$sample_tumor."_".get_job_id();
		    my $bash_file = $sample_tumor_job_dir."/".$job_id.".sh";

		    open CHECK_SH, ">$bash_file" or die "cannot open file $bash_file \n";
		    print CHECK_SH "#!/bin/bash\n\n";
		    print CHECK_SH "echo \"Start Check\t\" `date` `uname -n` >> $sample_tumor_log_dir/check.log\n\n";
		    print CHECK_SH "if [[ ";
		    if($opt{CNV_CONTRA} eq "yes"){print CHECK_SH "-f $sample_tumor_log_dir/contra.done && -f $sample_tumor_log_dir/contra_visualization.done && "}
		    if($opt{CNV_FREEC} eq "yes"){print CHECK_SH "-f $sample_tumor_log_dir/freec.done"}
		    print CHECK_SH " ]]\n";
		    print CHECK_SH "then\n";
		    print CHECK_SH "\ttouch $sample_tumor_log_dir/$sample_tumor_name.done\n";
		    print CHECK_SH "fi\n\n";
		    print CHECK_SH "echo \"End Check\t\" `date` `uname -n` >> $sample_tumor_log_dir/check.log\n";
		    close CHECK_SH;
		    if ( @cnv_jobs ){
			system "qsub -q $opt{CNVCHECK_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CNVCHECK_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id -hold_jid ".join(",",@cnv_jobs)." $bash_file";
		    } else {
			system "qsub -q $opt{CNVCHECK_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CNVCHECK_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id $bash_file";
		    }
		    push(@check_cnv_jobs, $job_id);
		}
	    }
	}
    ### Sample mode without control or ref sample
    } elsif($opt{CNV_MODE} eq "sample"){
	### Loop over samples
	foreach my $sample (@{$opt{SAMPLES}}){
	    my @cnv_jobs;
	    ## Create output, log and job directories
	    my $sample_out_dir = "$opt{OUTPUT_DIR}/copyNumber/$sample";
	    my $sample_log_dir = "$sample_out_dir/logs/";
	    my $sample_job_dir = "$sample_out_dir/jobs/";

		if(! -e $sample_out_dir){
		    make_path($sample_out_dir) or die "Couldn't create directory:  $sample_out_dir\n";
		}
		if(! -e $sample_job_dir){
		    make_path($sample_job_dir) or die "Couldn't create directory: $sample_job_dir\n";
		}
		if(! -e $sample_log_dir){
		    make_path($sample_log_dir) or die "Couldn't create directory: $sample_log_dir\n";
		}

		## Lookup running jobs and bams
		my $sample_bam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
		my @running_jobs;
		if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
		    push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample}});
		}
		## Print sample and bam info
		print "\n$sample \t $sample_bam \n";

		## Skip Copy number tools if .done file exist
		if (-e "$sample_log_dir/$sample.done"){
		    print "WARNING: $sample_log_dir/$sample.done exists, skipping \n";
		    next;
		}

		## Run CNV callers only Freec for wgs data
		if($opt{CNV_FREEC} eq "yes"){
		    print "\n###SCHEDULING FREEC####\n";
		    my $freec_job = runFreec($sample, $sample_out_dir, $sample_job_dir, $sample_log_dir, $sample_bam, "", \@running_jobs, \%opt);
		    if($freec_job){push(@cnv_jobs, $freec_job)};
		}
		## Check copy number analysis
		my $job_id = "CHECK_".$sample."_".get_job_id();
		my $bash_file = $sample_job_dir."/".$job_id.".sh";

		open CHECK_SH, ">$bash_file" or die "cannot open file $bash_file \n";
		print CHECK_SH "#!/bin/bash\n\n";
		print CHECK_SH "echo \"Start Check\t\" `date` `uname -n` >> $sample_log_dir/check.log\n\n";
		print CHECK_SH "if [[ -f $sample_log_dir/freec.done ]]\n";
		print CHECK_SH "then\n";
		print CHECK_SH "\ttouch $sample_log_dir/$sample.done\n";
		print CHECK_SH "fi\n\n";
		print CHECK_SH "echo \"End Check\t\" `date` `uname -n` >> $sample_log_dir/check.log\n";
		close CHECK_SH;
	    
		if ( @cnv_jobs ){
		    system "qsub -q $opt{CNVCHECK_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CNVCHECK_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_log_dir -e $sample_log_dir -N $job_id -hold_jid ".join(",",@cnv_jobs)." $bash_file";
		} else {
		    system "qsub -q $opt{CNVCHECK_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CNVCHECK_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_log_dir -e $sample_log_dir -N $job_id $bash_file";
		}
		push(@check_cnv_jobs, $job_id);
	}
    }
    return \@check_cnv_jobs;
}

### Copy number analysis tools
sub runFreec {
    ###
    # Run freec and plot result
    ###
    my ($sample_name, $out_dir, $job_dir, $log_dir, $sample_bam, $control_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    
    ## Skip Contra if .done file exist
    if (-e "$log_dir/freec.done"){
	print "WARNING: $log_dir/freec.done exists, skipping \n";
	return;
    }

    ## Create FREEC output directory
    my $freec_out_dir = "$out_dir/freec";
    if(! -e $freec_out_dir){
	make_path($freec_out_dir) or die "Couldn't create directory: $freec_out_dir\n";
    }
    
    ## Create freec config
    my $freec_config = $freec_out_dir."/freec_config.txt";
    open FREEC_CONFIG, ">$freec_config" or die "cannot open file $freec_config \n";

    print FREEC_CONFIG "[general]\n";
    print FREEC_CONFIG "chrLenFile= $opt{FREEC_CHRLENFILE}\n";
    print FREEC_CONFIG "ploidy=2\n";
    print FREEC_CONFIG "samtools=$opt{SAMTOOLS_PATH}/samtools\n";
    print FREEC_CONFIG "chrFiles= $opt{FREEC_CHRFILES}\n";
    print FREEC_CONFIG "window=$opt{FREEC_WINDOW}\n";
    print FREEC_CONFIG "maxThreads=$opt{FREEC_THREADS}\n";
    print FREEC_CONFIG "telocentromeric=$opt{FREEC_TELOCENTROMERIC}\n";
    print FREEC_CONFIG "BedGraphOutput=TRUE\n";
    print FREEC_CONFIG "outputDir=$freec_out_dir\n";
    print FREEC_CONFIG "gemMappabilityFile=/hpc/local/CentOS6/cog_bioinf/freec/hg19_mappability_tracks/out100m1_hg19.gem\n";
    print FREEC_CONFIG "gemMappabilityFile=/hpc/local/CentOS6/cog_bioinf/freec/hg19_mappability_tracks/out100m2_hg19.gem\n";
    
    print FREEC_CONFIG "[sample]\n";
    print FREEC_CONFIG "mateFile=$sample_bam\n";
    print FREEC_CONFIG "inputFormat=BAM\n";
    print FREEC_CONFIG "mateOrientation=FR\n";
    if($control_bam){
	print FREEC_CONFIG "[control]\n";
	print FREEC_CONFIG "mateFile=$control_bam\n";
	print FREEC_CONFIG "inputFormat=BAM\n";
	print FREEC_CONFIG "mateOrientation=FR\n";
    }

    if ($opt{CNV_TARGETS}){
	print FREEC_CONFIG "[target]\n";
	print FREEC_CONFIG "captureRegions=$opt{CNV_TARGETS}\n";
    }

    close FREEC_CONFIG;

    ## Create freec bash script
    my $job_id = "FREEC_".$sample_name."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";
    my $sample_bam_name = (split('/',$sample_bam))[-1];
    my $control_bam_name = (split('/',$control_bam))[-1];

    open FREEC_SH, ">$bash_file" or die "cannot open file $bash_file \n";

    print FREEC_SH "#!/bin/bash\n\n";
    if($control_bam){
	print FREEC_SH "if [ -s $sample_bam -a -s $control_bam ]\n";
    } else {
	print FREEC_SH "if [ -s $sample_bam ]\n";
    }
    print FREEC_SH "then\n";
    print FREEC_SH "\techo \"Start FREEC\t\" `date` \"\t $sample_bam \t $control_bam\t\" `uname -n` >> $log_dir/freec.log\n\n";

    print FREEC_SH "\t$opt{FREEC_PATH}/freec -conf $freec_config\n";
    print FREEC_SH "cd $freec_out_dir\n";
    print FREEC_SH "\tcat $opt{FREEC_PATH}/assess_significance.R | R --slave --args ".$sample_bam_name."_CNVs ".$sample_bam_name."_ratio.txt\n";
    print FREEC_SH "\tcat $opt{FREEC_PATH}/makeGraph.R | R --slave --args 2 ".$sample_bam_name."_ratio.txt\n";
    print FREEC_SH "touch $log_dir/freec.done\n";
    print FREEC_SH "\techo \"End FREEC\t\" `date` \"\t $sample_bam \t $control_bam\t\" `uname -n` >> $log_dir/freec.log\n\n";
    print FREEC_SH "else\n";
    print FREEC_SH "\techo \"ERROR: $sample_bam or $control_bam does not exist.\" >> $log_dir/freec.log\n";
    print FREEC_SH "fi\n";
    
    close FREEC_SH;
    
    ## Run job
    if ( @running_jobs ){
	system "qsub -q $opt{FREEC_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{FREEC_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{FREEC_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{FREEC_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
}

sub runContra {
    ###
    # Run contra
    ###
    my ($sample_tumor, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $contra_out_dir = "$out_dir/contra";

    ## Skip Contra if .done file exist
    if (-e "$log_dir/contra.done"){
	print "WARNING: $log_dir/contra.done exists, skipping \n";
	return;
    }

    ## Create contra bash script
    my $job_id = "CNTR_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    open CONTRA_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print CONTRA_SH "#!/bin/bash\n\n";
    print CONTRA_SH "if [ -s $sample_tumor_bam -a -s $sample_ref_bam ]\n";
    print CONTRA_SH "then\n";
    print CONTRA_SH "\techo \"Start Contra\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/contra.log\n\n";

    # Run CONTRA
    print CONTRA_SH "\t$opt{CONTRA_PATH}/contra.py -s $sample_tumor_bam -c $sample_ref_bam -f $opt{GENOME} -t $opt{CNV_TARGETS} -o $contra_out_dir/ --sampleName $sample_tumor $opt{CONTRA_FLAGS} \n\n";
    print CONTRA_SH "\tcd $contra_out_dir\n";
    # Check contra completed
    print CONTRA_SH "\tif [ -s $contra_out_dir/table/$sample_tumor*.vcf ]\n";
    print CONTRA_SH "\tthen\n";
    print CONTRA_SH "\t\ttouch $log_dir/contra.done\n";
    print CONTRA_SH "\tfi\n\n";

    print CONTRA_SH "\techo \"End Contra\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/contra.log\n\n";

    print CONTRA_SH "else\n";
    print CONTRA_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >> $log_dir/contra.log\n";
    print CONTRA_SH "fi\n";

    close CONTRA_SH;

    ## Run job
    if ( @running_jobs ){
	system "qsub -q $opt{CONTRA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CONTRA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{CONTRA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CONTRA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
}

sub runContraVisualization {
    ###
    # Run contra visualization using the contra plotscript create by Annelies Smouter
    ###
    my ($sample_tumor, $out_dir, $job_dir, $log_dir, $contra_job, $opt) = (@_);
    my %opt = %{$opt};
    my $contra_out_dir = "$out_dir/contra";

    ## Skip Contra if .done file exist
    if (-e "$log_dir/contra_visualization.done"){
	print "WARNING: $log_dir/contra_visualization.done exists, skipping \n";
	return;
    }
    my $job_id = "CNTRVIS_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    open CONTRAVIS_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print CONTRAVIS_SH "#!/bin/bash\n\n";
    print CONTRAVIS_SH "if [ -f $log_dir/contra.done ]\n";
    print CONTRAVIS_SH "\tthen\n";
    print CONTRAVIS_SH "\tperl $opt{CONTRA_PLOTSCRIPT} -input $contra_out_dir/table/$sample_tumor.CNATable.10rd.10bases.20bins.txt -d $opt{CONTRA_PLOTDESIGN}\n";
    print CONTRAVIS_SH "\tchmod 644 $contra_out_dir/table/*\n";
    print CONTRAVIS_SH "touch $log_dir/contra_visualization.done\n";
    print CONTRAVIS_SH "fi\n";
    close CONTRAVIS_SH;

    ## Run job
    if ( $contra_job ){
	system "qsub -q $opt{CONTRA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CONTRA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid $contra_job $bash_file";
    } else {
	system "qsub -q $opt{CONTRA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CONTRA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
}


############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
