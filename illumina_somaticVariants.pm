#!/usr/bin/perl -w

#########################################################
### illumina_somaticVariants.pm
### - Run somatic variant callers
###   - Varscan, Strelka, FreeBayes
### - Merge and annotate somatic high confidence calls.
###
### Author: R.F.Ernst
#########################################################

package illumina_somaticVariants;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);

sub parseSamples {
    ###
    # Parse sample names
    # Expects CPCT samples (CPCT........T/R)
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my %somatic_samples;
    my @pileupJobs;

    foreach my $sample (@{$opt{SAMPLES}}){
	# Parse cpct samples based on regular expression defining two groups, sample name and sample origin.
	my ($cpct_name,$origin) = ($sample =~ /$opt{SOMVAR_REGEX}/);

	if ( (! $cpct_name) || (! $origin) ){
	    print "WARNING: $sample is not passing somatic samplename parsing, skipping \n\n";
	    next;
	# Run pileup for varscan
	} else {
	    if($opt{SOMVAR_VARSCAN} eq "yes"){
		print "Creating pileup for: $sample\n";
		my $pileup_job = runPileup($sample, \%opt);
		push(@pileupJobs, $pileup_job);
	    }
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
    $opt{RUNNING_JOBS}->{'pileup'} = \@pileupJobs;
    return \%opt;
}

### Run and merge
sub runSomaticVariantCallers {
    ###
    # Run somatic variant callers
    # Merge and annotate high confidence calls
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @merge_somvar_jobs;
    ### Loop over tumor samples
    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){

	# Check correct sample ref
	if (! $opt{SOMATIC_SAMPLES}{$sample}{'ref'}){
	    print "WARNING: No ref sample for $sample, skipping \n";
	    next;
	}

	foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){
	    foreach my $sample_ref (@{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}){
		my @somvar_jobs;
		## Create output, log and job directories
		my $sample_tumor_name = "$sample_ref\_$sample_tumor";
		my $sample_tumor_out_dir = "$opt{OUTPUT_DIR}/somaticVariants/$sample_tumor_name";
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

		## Skip Somatic Callers if .done file exist
		if (-e "$sample_tumor_log_dir/$sample_tumor_name.done"){
		    print "WARNING: $sample_tumor_log_dir/$sample_tumor_name.done exists, skipping \n";
		    next;
		}

		## Run somatic callers
		if($opt{SOMVAR_STRELKA} eq "yes"){
		    print "\n###SCHEDULING STRELKA####\n";
		    my $strelka_job = runStrelka($sample_tumor, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		    if($strelka_job){push(@somvar_jobs, $strelka_job)};
		}
		if($opt{SOMVAR_VARSCAN} eq "yes"){
		    print "\n###SCHEDULING VARSCAN####\n";
		    my $varscan_job = runVarscan($sample_tumor, $sample_tumor_name, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		    if($varscan_job){push(@somvar_jobs, $varscan_job)};
		}
		if($opt{SOMVAR_FREEBAYES} eq "yes"){
		    print "\n###SCHEDULING FREEBAYES####\n";
		    my $freebayes_job = runFreeBayes($sample_tumor, $sample_tumor_name, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		    if($freebayes_job){push(@somvar_jobs, $freebayes_job)};
		}
		if($opt{SOMVAR_MUTECT} eq "yes"){
		    print "\n###SCHEDULING MUTECT####\n";
		    my $mutect_job = runMutect($sample_tumor, $sample_tumor_name, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		    if($mutect_job){push(@somvar_jobs, $mutect_job)};
		}
		## Merge somatic vcfs
		print "\n###SCHEDULING MERGE SOMATIC VCFS####\n";

		my $job_id = "MERGE_".$sample_tumor."_".get_job_id();
		my $bash_file = $sample_tumor_job_dir."/".$job_id.".sh";

		open MERGE_SH, ">$bash_file" or die "cannot open file $bash_file \n";
		print MERGE_SH "#!/bin/bash\n\n";
		print MERGE_SH "echo \"Start Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";

		# Merge vcfs
		my $invcf;
		my $outvcf = "$sample_tumor_out_dir/$sample_tumor_name\_merged_somatics.vcf";
		print MERGE_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} -o $outvcf --genotypemergeoption uniquify ";
		if($opt{SOMVAR_STRELKA} eq "yes"){ print MERGE_SH "-V:strelka $sample_tumor_out_dir/strelka/passed.somatic.merged.vcf "; }
		if($opt{SOMVAR_VARSCAN} eq "yes"){ print MERGE_SH "-V:varscan $sample_tumor_out_dir/varscan/$sample_tumor_name.merged.Somatic.hc.vcf "; }
		if($opt{SOMVAR_FREEBAYES} eq "yes"){ print MERGE_SH "-V:freebayes $sample_tumor_out_dir/freebayes/$sample_tumor_name\_somatic_filtered.vcf "; }
		if($opt{SOMVAR_MUTECT} eq "yes"){ print MERGE_SH "-V:mutect $sample_tumor_out_dir/mutect/$sample_tumor_name\_mutect_passed.vcf ";}

		# Filter vcf on target
		if($opt{SOMVAR_TARGETS}){
		    $invcf = $outvcf;
		    $outvcf = "$sample_tumor_out_dir/$sample_tumor_name\_filtered_merged_somatics.vcf";
		    print MERGE_SH "\n\njava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T SelectVariants -R $opt{GENOME} -L $opt{SOMVAR_TARGETS} -V $invcf -o $outvcf\n";
		    print MERGE_SH "rm $invcf*";
		}

		# Annotate somatic vcf
		if($opt{SOMVAR_ANNOTATE} eq "yes"){
		    $invcf = $outvcf;
		    my $preAnnotateVCF = $invcf;
		    $outvcf =~ s/.vcf/_snpEff.vcf/;
		    print MERGE_SH "\n\njava -Xmx6G -jar $opt{SNPEFF_PATH}/snpEff.jar -c $opt{SNPEFF_PATH}/snpEff.config $opt{ANNOTATE_DB} -v $invcf $opt{ANNOTATE_FLAGS} > $outvcf\n";
		
		    ## dbsnp
		    $invcf = $outvcf;
		    my $suffix = "_dbSNP.vcf";
		    $outvcf =~ s/.vcf/$suffix/;
		    print MERGE_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T VariantAnnotator -nt $opt{SOMVARMERGE_THREADS} -R $opt{GENOME} -o $outvcf --variant $invcf --dbsnp $opt{CALLING_DBSNP} --alwaysAppendDbsnpId\n";
		    print MERGE_SH "if [ -s $outvcf ]\nthen\n\trm $invcf $invcf.idx \nfi\n";
		
		    ## cosmic
		    $invcf = $outvcf;
		    $suffix = "_$opt{ANNOTATE_IDNAME}.vcf";
		    $outvcf =~ s/.vcf/$suffix/;
		    print MERGE_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T VariantAnnotator -nt $opt{SOMVARMERGE_THREADS} -R $opt{GENOME} -o $outvcf --variant $invcf --dbsnp $opt{ANNOTATE_IDDB} --alwaysAppendDbsnpId\n";
		    print MERGE_SH "if [ -s $outvcf ]\nthen\n\trm $invcf $invcf.idx \nfi\n";
		
		    ## Check annotated vcf using the last position
		    print MERGE_SH "\nif [ \"\$(tail -n 1 $preAnnotateVCF | cut -f 1,2)\" = \"\$(tail -n 1 $outvcf | cut -f 1,2)\" -a -s $preAnnotateVCF -a -s $outvcf ]\n";
		    print MERGE_SH "then\n";
		    print MERGE_SH "\ttouch $sample_tumor_log_dir/$sample_tumor_name.done\n";
		    print MERGE_SH "fi\n";
		    print MERGE_SH "echo \"END Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";
		    close MERGE_SH;
		} else {
		    print MERGE_SH "\nif [ -s $outvcf ]\n";
		    print MERGE_SH "then\n";
		    print MERGE_SH "\ttouch $sample_tumor_log_dir/$sample_tumor_name.done\n";
		    print MERGE_SH "fi\n";
		    print MERGE_SH "echo \"END Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";
		    close MERGE_SH;
		}

		# Run job
		if ( @somvar_jobs ){
		    system "qsub -q $opt{SOMVARMERGE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{SOMVARMERGE_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id -hold_jid ".join(",",@somvar_jobs)." $bash_file";
		} else {
		    system "qsub -q $opt{SOMVARMERGE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{SOMVARMERGE_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id $bash_file";
		}
		push(@merge_somvar_jobs, $job_id);
	    }
	}
    }
    return \@merge_somvar_jobs;
}

### Somatic Variant Callers
sub runStrelka {
    my ($sample_tumor, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $strelka_out_dir = "$out_dir/strelka";

    ## Skip Strelka if .done file exist
    if (-e "$log_dir/strelka.done"){
	print "WARNING: $log_dir/strelka.done, skipping \n";
	return;
    }

    ## Create strelka bash script
    my $job_id = "STR_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    open STRELKA_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print STRELKA_SH "#!/bin/bash\n\n";
    print STRELKA_SH "if [ -s $sample_tumor_bam -a -s $sample_ref_bam ]\n";
    print STRELKA_SH "then\n";
    print STRELKA_SH "\techo \"Start Strelka\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/strelka.log\n\n";

    # Run Strelka
    print STRELKA_SH "\t$opt{STRELKA_PATH}/bin/configureStrelkaWorkflow.pl --tumor $sample_tumor_bam --normal $sample_ref_bam --ref $opt{GENOME} --config $opt{STRELKA_INI} --output-dir $strelka_out_dir\n\n";

    print STRELKA_SH "\tcd $strelka_out_dir\n";
    print STRELKA_SH "\tmake -j 8\n\n";

    # Check strelka completed
    print STRELKA_SH "\tif [ -f $strelka_out_dir/task.complete ]\n";
    print STRELKA_SH "\tthen\n";
    print STRELKA_SH "\t\tjava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} --genotypemergeoption unsorted -o passed.somatic.merged.vcf -V results/passed.somatic.snvs.vcf -V results/passed.somatic.indels.vcf \n";
    print STRELKA_SH "\t\tperl -p -e 's/\\t([A-Z][A-Z]:)/\\tGT:\$1/g' passed.somatic.merged.vcf | perl -p -e 's/(:T[UO]R?)\\t/\$1\\t0\\/0:/g' | perl -p -e 's/(:\\d+,\\d+)\\t/\$1\\t0\\/1:/g' | perl -p -e 's/(#CHROM.*)/##StrelkaGATKCompatibility=Added GT fields to strelka calls for gatk compatibility.\\n\$1/g' > temp.vcf\n";
    print STRELKA_SH "\t\tmv temp.vcf passed.somatic.merged.vcf\n";
    print STRELKA_SH "\t\trm -r chromosomes/ \n";
    print STRELKA_SH "\t\ttouch $log_dir/strelka.done\n";
    print STRELKA_SH "\tfi\n\n";
    print STRELKA_SH "\techo \"End Strelka\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/strelka.log\n\n";
    
    print STRELKA_SH "else\n";
    print STRELKA_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print STRELKA_SH "fi\n";
    
    close STRELKA_SH;

    ## Run job
    if ( @running_jobs ){
	system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }

    return $job_id;
}

sub runPileup {
    my ($sample, $configuration) = (@_);
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];

    my $bam = $opt{BAM_FILES}->{$sample};
    (my $pileup = $bam) =~ s/\.bam/\.pileup/;
    my $jobID = "PILEUP_$sample\_".get_job_id();
    
    ## Check for pileup.done file
    if (-e "$opt{OUTPUT_DIR}/$sample/logs/Pileup_$sample.done"){
	print "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/Pileup_$sample.done exists, skipping\n";
	return $jobID;
    }

    ## Pileup command
    my $pileup_command = "$opt{SAMBAMBA_PATH}/sambamba mpileup -t $opt{PILEUP_THREADS} --tmpdir=$opt{OUTPUT_DIR}/$sample/tmp/ ";
    if ( $opt{SOMVAR_TARGETS} ) {
	$pileup_command .= "-L $opt{SOMVAR_TARGETS} ";
    }
    $pileup_command .= "$opt{OUTPUT_DIR}/$sample/mapping/$bam --samtools \"-q 1 -f $opt{GENOME}\" | $opt{TABIX_PATH}/bgzip -c > $pileup.gz";
    
    # TABIX
    my $tabix_command = "$opt{TABIX_PATH}/tabix -s 1 -b 2 -e 2 $pileup.gz";

    ## Create pileup bash script
    my $logDir = $opt{OUTPUT_DIR}."/".$sample."/logs";
    my $bashFile = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobID.".sh";
    
    open PILEUP_SH,">$bashFile" or die "Couldn't create $bashFile\n";
    print PILEUP_SH "\#!/bin/sh\n\n";
    print PILEUP_SH "cd $opt{OUTPUT_DIR}/$sample/tmp\n";
    print PILEUP_SH "echo \"Start pileup\t\" `date` \"\t$bam\t\" `uname -n` >> $logDir/$sample.log\n\n";
    print PILEUP_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$bam ]\n";
    print PILEUP_SH "then\n";
    print PILEUP_SH "\tPATH=$opt{SAMTOOLS_PATH}:\$PATH\n";
    print PILEUP_SH "\texport PATH\n";
    print PILEUP_SH "\t$pileup_command\n";
    print PILEUP_SH "\t$tabix_command\n";
    print PILEUP_SH "\tif [ \"\$($opt{TABIX_PATH}/tabix $pileup.gz MT | tail -n 1 | cut -f 1)\" = \"MT\" ]\n";
    print PILEUP_SH "\tthen\n";
    print PILEUP_SH "\t\tmv $pileup.gz* $opt{OUTPUT_DIR}/$sample/mapping/\n";
    print PILEUP_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/Pileup_$sample.done\n";
    print PILEUP_SH "\telse\n";
    print PILEUP_SH "\t\techo \"ERROR: $pileup seems incomplete, it does not end with MT\" >&2\n";
    print PILEUP_SH "\tfi\n";
    print PILEUP_SH "else\n";
    print PILEUP_SH "\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$bam does not exist.\" >&2\n";
    print PILEUP_SH "fi\n\n";
    print PILEUP_SH "echo \"END pileup\t\" `date` \"\t$bam\t\" `uname -n` >> $logDir/$sample.log\n";
    close PILEUP_SH;
    
    ### Submit realign bash script
    if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
	system "qsub -q $opt{PILEUP_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{PILEUP_THREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/Pileup_$sample.out -e $logDir/Pileup_$sample.err -N $jobID -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." $bashFile";
    } else {
	system "qsub -q $opt{PILEUP_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{PILEUP_THREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/Pileup_$sample.out -e $logDir/Pileup_$sample.err -N $jobID $bashFile";
    }
    return $jobID;
}

sub runVarscan {
    my ($sample_tumor, $sample_tumor_name, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my %opt = %{$opt};
    my @running_jobs = @{$running_jobs};
    push(@running_jobs, @{$opt{RUNNING_JOBS}->{'pileup'}});
    my $varscan_out_dir = "$out_dir/varscan";
    (my $sample_tumor_pileup = $sample_tumor_bam) =~ s/\.bam/\.pileup\.gz/;
    (my $sample_ref_pileup = $sample_ref_bam) =~ s/\.bam/\.pileup\.gz/;

    ## Create output dir
    if( ! -e $varscan_out_dir ){
	make_path($varscan_out_dir) or die "Couldn't create directory: $varscan_out_dir\n";
    }

    ## Skip varscan if .done file exist
    if ( -e "$log_dir/varscan.done" ){
	print "WARNING: $log_dir/varscan.done exists, skipping \n";
	return;
    }

    ## Run varscan per chromosome
    my $dictFile = $opt{GENOME};
    $dictFile =~ s/.fasta$/.dict/;
    my @chrs = @{get_chrs_from_dict($dictFile)};
    my @varscan_jobs;

    foreach my $chr (@chrs){
	## ADD: Chunk done check and skip if done.
	my $job_id = "VS_".$sample_tumor."_".$chr."_".get_job_id();
	my $bash_file = $job_dir."/".$job_id.".sh";
	my $output_name = $sample_tumor_name."_".$chr;
	open VARSCAN_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	print VARSCAN_SH "#!/bin/bash\n\n";
	print VARSCAN_SH "cd $varscan_out_dir\n";
	print VARSCAN_SH "if [ -s $sample_ref_pileup -a -s $sample_tumor_pileup ]\n";
	print VARSCAN_SH "then\n";

	print VARSCAN_SH "\techo \"Start Varscan\t\" `date` \"\t $chr \t $sample_ref_pileup \t $sample_tumor_pileup\t\" `uname -n` >> $log_dir/varscan.log\n";
	print VARSCAN_SH "\tjava -Xmx12g -jar $opt{VARSCAN_PATH} somatic <($opt{TABIX_PATH}/tabix $sample_ref_pileup $chr) <($opt{TABIX_PATH}/tabix $sample_tumor_pileup $chr) $output_name $opt{VARSCAN_SETTINGS} --output-vcf 1\n\n";
	print VARSCAN_SH "\techo \"End Varscan\t\" `date` \"\t $chr $sample_ref_pileup \t $sample_tumor_pileup\t\" `uname -n` >> $log_dir/varscan.log\n";
	print VARSCAN_SH "else\n";
	print VARSCAN_SH "\techo \"ERROR: $sample_tumor_pileup or $sample_ref_pileup does not exist.\" >&2\n";
	print VARSCAN_SH "fi\n";
	close VARSCAN_SH;

	## Run job
	if ( @running_jobs ){
	    system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
	} else {
	    system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
	}
	
	push(@varscan_jobs,$job_id);
    }

    ## Concat chromosome vcfs 
    my $job_id = "VS_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    # Setup test, concat and rm of chr chunks
    my $file_test = "if [ -s $sample_ref_bam -a -s $sample_tumor_bam ";
    my $snp_concat_command = "$opt{VCFTOOLS_PATH}/vcf-concat ";
    my $indel_concat_command = "$opt{VCFTOOLS_PATH}/vcf-concat ";
    my $rm_command = "rm ";

    foreach my $chr (@chrs){
	my $snp_output = $sample_tumor_name."_".$chr.".snp.vcf";
	my $indel_output = $sample_tumor_name."_".$chr.".indel.vcf";
	$file_test .= "-a -s $snp_output -a -s $indel_output ";
	$snp_concat_command .= "$snp_output ";
	$indel_concat_command .= "$indel_output ";
	$rm_command .= "$snp_output $indel_output ";
    }
    $file_test .= "]";
    $snp_concat_command .= "> $sample_tumor_name.snp.vcf";
    $indel_concat_command .= "> $sample_tumor_name.indel.vcf";

    # Create bash script
    open VARSCAN_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print VARSCAN_SH "#!/bin/bash\n\n";

    print VARSCAN_SH "cd $varscan_out_dir\n";
    print VARSCAN_SH "$file_test\n";
    print VARSCAN_SH "then\n";
    print VARSCAN_SH "\techo \"Start concat and postprocess Varscan\t\" `date` \"\t $sample_ref_pileup \t $sample_tumor_pileup\t\" `uname -n` >> $log_dir/varscan.log\n";
    print VARSCAN_SH "\t$snp_concat_command\n";
    print VARSCAN_SH "\t$indel_concat_command\n\n";

    # postprocessing
    print VARSCAN_SH "\tjava -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $sample_tumor_name.indel.vcf $opt{VARSCAN_POSTSETTINGS}\n";
    print VARSCAN_SH "\tjava -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $sample_tumor_name.snp.vcf $opt{VARSCAN_POSTSETTINGS}\n\n";

    # merge varscan hc snps and indels
    print VARSCAN_SH "\tjava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} --genotypemergeoption unsorted -o $sample_tumor_name.merged.Somatic.hc.vcf -V $sample_tumor_name.snp.Somatic.hc.vcf -V $sample_tumor_name.indel.Somatic.hc.vcf\n";
    print VARSCAN_SH "\tsed -i 's/SSC/VS_SSC/' $sample_tumor_name.merged.Somatic.hc.vcf\n\n"; # to resolve merge conflict with FB vcfs

    # Check varscan completed
    print VARSCAN_SH "\tif [ -s $sample_tumor_name.merged.Somatic.hc.vcf ]\n";
    print VARSCAN_SH "\tthen\n";
    print VARSCAN_SH "\t\t$rm_command\n"; #remove tmp chr vcf files
    print VARSCAN_SH "\t\ttouch $log_dir/varscan.done\n";
    print VARSCAN_SH "\tfi\n\n";
    print VARSCAN_SH "\techo \"END concat and postprocess Varscan\t\" `date` \"\t $sample_ref_pileup \t $sample_tumor_pileup\t\" `uname -n` >> $log_dir/varscan.log\n";

    print VARSCAN_SH "else\n";
    print VARSCAN_SH "\techo \"ERROR: $sample_tumor_pileup or $sample_ref_pileup does not exist.\" >&2\n";
    print VARSCAN_SH "fi\n";
    close VARSCAN_SH;

    ## Run job
    if ( @varscan_jobs ){
        system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@varscan_jobs)." $bash_file";
    } else {
        system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
}

sub runFreeBayes {
    my ($sample_tumor, $sample_tumor_name, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $freebayes_out_dir = "$out_dir/freebayes";

    ## Create output dir
    if(! -e $freebayes_out_dir){
	make_path($freebayes_out_dir) or die "Couldn't create directory: $freebayes_out_dir\n";
    }

    ## Skip freebayes if .done file exist
    if (-e "$log_dir/freebayes.done"){
	print "WARNING: $log_dir/freebayes.done exists, skipping \n";
	return;
    }

    ## Run freebayes per chromosome
    my $dictFile = $opt{GENOME};
    $dictFile =~ s/.fasta$/.dict/;
    my @chrs = @{get_chrs_from_dict($dictFile)};
    my @freebayes_jobs;

    foreach my $chr (@chrs){
	## ADD: Chunk done check and skip if done.
	my $job_id = "FB_".$sample_tumor."_".$chr."_".get_job_id();
	my $bash_file = $job_dir."/".$job_id.".sh";
	my $output_name = $sample_tumor_name."_".$chr;

	## Create freebayes command
	my $freebayes_command = "$opt{FREEBAYES_PATH}/freebayes -f $opt{GENOME} -r $chr ";
	$freebayes_command .= "$opt{FREEBAYES_SETTINGS} $sample_ref_bam $sample_tumor_bam > $freebayes_out_dir/$output_name.vcf";

	## Sort vcf, remove duplicate lines and filter on target
	my $sort_uniq_filter_command = "$opt{VCFTOOLS_PATH}/vcf-sort -c $freebayes_out_dir/$output_name.vcf | $opt{VCFLIB_PATH}/vcfuniq > $freebayes_out_dir/$output_name.sorted_uniq.vcf";
	my $mv_command;
	# Filter vcf on target
	if($opt{SOMVAR_TARGETS}){
	    $sort_uniq_filter_command .= "\n\tjava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T SelectVariants -R $opt{GENOME} -L $opt{SOMVAR_TARGETS} -V $freebayes_out_dir/$output_name.sorted_uniq.vcf -o $freebayes_out_dir/$output_name.sorted_uniq_targetfilter.vcf\n";
	    $mv_command = "mv $freebayes_out_dir/$output_name.sorted_uniq_targetfilter.vcf $freebayes_out_dir/$output_name.vcf";
	} else {
	    $mv_command = "mv $freebayes_out_dir/$output_name.sorted_uniq.vcf $freebayes_out_dir/$output_name.vcf";
	}
	## Create bashscript
	open FREEBAYES_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	print FREEBAYES_SH "#!/bin/bash\n\n";
	print FREEBAYES_SH "cd $freebayes_out_dir\n";
	print FREEBAYES_SH "if [ -s $sample_tumor_bam -a -s $sample_ref_bam ]\n";
	print FREEBAYES_SH "then\n";
	print FREEBAYES_SH "\techo \"Start Freebayes\t\" `date` \"\t $chr \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n\n";
	print FREEBAYES_SH "\t$freebayes_command\n";
	print FREEBAYES_SH "\t$sort_uniq_filter_command\n";
	print FREEBAYES_SH "\t$mv_command\n\n";
	print FREEBAYES_SH "\techo \"End Freebayes\t\" `date` \"\t $chr $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n";
	print FREEBAYES_SH "else\n";
	print FREEBAYES_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
	print FREEBAYES_SH "fi\n";
	close FREEBAYES_SH;

	## Run job
	if ( @running_jobs ){
	    system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
	} else {
	    system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
	}
	
	push(@freebayes_jobs,$job_id);
    }

    ## Concat chromosome vcfs and postprocess vcf
    my $job_id = "FB_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    # Setup test, concat and rm of chr chunks
    my $file_test = "if [ -s $sample_ref_bam -a -s $sample_tumor_bam ";
    my $concat_command = "$opt{VCFTOOLS_PATH}/vcf-concat ";
    my $rm_command = "rm ";
    foreach my $chr (@chrs){
	my $snp_output = $sample_tumor_name."_".$chr;
	$file_test .= "-a -s $snp_output\.vcf ";
	$concat_command .= "$snp_output\.vcf ";
	$rm_command .= "$snp_output\* ";
    }
    $file_test .= "]";
    $concat_command .= "> $sample_tumor_name.vcf";

    # Create bash script
    open FREEBAYES_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print FREEBAYES_SH "#!/bin/bash\n\n";

    print FREEBAYES_SH "cd $freebayes_out_dir\n";
    print FREEBAYES_SH "$file_test\n";
    print FREEBAYES_SH "then\n";
    print FREEBAYES_SH "\techo \"Start concat and postprocess Freebayes\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n";
    print FREEBAYES_SH "\t$concat_command\n\n";

    # Uniqify freebayes output 
    print FREEBAYES_SH "\tuniq $freebayes_out_dir/$sample_tumor_name.vcf > $freebayes_out_dir/$sample_tumor_name.uniq.vcf\n";
    print FREEBAYES_SH "\tmv $freebayes_out_dir/$sample_tumor_name.uniq.vcf $freebayes_out_dir/$sample_tumor_name.vcf\n\n";

    # get sample ids
    print FREEBAYES_SH "\tsample_R=`grep -P \"^#CHROM\" $freebayes_out_dir/$sample_tumor_name.vcf | cut -f 10`\n";
    print FREEBAYES_SH "\tsample_T=`grep -P \"^#CHROM\" $freebayes_out_dir/$sample_tumor_name.vcf | cut -f 11`\n\n";

    # annotate somatic and germline scores
    print FREEBAYES_SH "\t$opt{VCFSAMPLEDIFF_PATH}/vcfsamplediff VT \$sample_R \$sample_T $freebayes_out_dir/$sample_tumor_name.vcf > $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n";
    print FREEBAYES_SH "\tsed -i 's/SSC/FB_SSC/' $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n"; # to resolve merge conflicts with varscan vcfs
    print FREEBAYES_SH "\tgrep -P \"^#\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf > $freebayes_out_dir/$sample_tumor_name\_germline.vcf\n";
    print FREEBAYES_SH "\tgrep -P \"^#\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf > $freebayes_out_dir/$sample_tumor_name\_somatic.vcf\n";
    print FREEBAYES_SH "\tgrep -i \"VT=germline\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf >> $freebayes_out_dir/$sample_tumor_name\_germline.vcf\n";
    print FREEBAYES_SH "\tgrep -i \"VT=somatic\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf >> $freebayes_out_dir/$sample_tumor_name\_somatic.vcf\n";
    print FREEBAYES_SH "\trm $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n\n";

    # Filter
    print FREEBAYES_SH "\tcat $freebayes_out_dir/$sample_tumor_name\_somatic.vcf | $opt{BIOVCF_PATH}/bio-vcf $opt{FREEBAYES_SOMATICFILTER} > $freebayes_out_dir/$sample_tumor_name\_somatic_filtered.vcf\n";
    print FREEBAYES_SH "\tcat $freebayes_out_dir/$sample_tumor_name\_germline.vcf | $opt{BIOVCF_PATH}/bio-vcf $opt{FREEBAYES_GERMLINEFILTER} > $freebayes_out_dir/$sample_tumor_name\_germline_filtered.vcf\n\n";
    
    #Check freebayes completed
    print FREEBAYES_SH "\tif [ -s $freebayes_out_dir/$sample_tumor_name\_somatic_filtered.vcf -a -s $freebayes_out_dir/$sample_tumor_name\_germline_filtered.vcf ]\n";
    print FREEBAYES_SH "\tthen\n";
    print FREEBAYES_SH "\t\t$rm_command\n";
    print FREEBAYES_SH "\t\ttouch $log_dir/freebayes.done\n\n"; ## Check on complete output!!!
    print FREEBAYES_SH "\tfi\n";
    
    print FREEBAYES_SH "\techo \"End concat and postprocess Freebayes\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n";
    print FREEBAYES_SH "else\n";
    print FREEBAYES_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print FREEBAYES_SH "fi\n";

    close FREEBAYES_SH;

    # Run job
    if ( @freebayes_jobs ){
	system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@freebayes_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
}

sub runMutect {
    my ($sample_tumor, $sample_tumor_name, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $mutect_out_dir = "$out_dir/mutect";
    my $mutect_tmp_dir = "$mutect_out_dir/tmp";

    ## Create output and tmp dir
    if(! -e $mutect_out_dir){
	make_path($mutect_out_dir) or die "Couldn't create directory: $mutect_out_dir\n";
    }
    if(! -e $mutect_tmp_dir){
	make_path($mutect_tmp_dir) or die "Couldn't create directory: $mutect_tmp_dir\n";
    }
    
    ## Skip Mutect if .done file exist
    if (-e "$log_dir/mutect.done"){
	print "WARNING: $log_dir/mutect.done, skipping \n";
	return;
    }

    ## Build Queue command 
    ## wait for GATK-MuTect integration, see http://gatkforums.broadinstitute.org/discussion/comment/24614#Comment_24614
    #my $command = "java -Xmx".$javaMem."G -Xms".$opt{MUTECT_MASTERMEM}."G -jar $opt{QUEUE_PATH}/Queue.jar ";
    #$command .= "-jobQueue $opt{MUTECT_QUEUE} -jobNative \"-pe threaded $opt{MUTECT_THREADS} -P $opt{CLUSTER_PROJECT}\" -jobRunner GridEngine -jobReport $log_dir/mutect.jobReport.txt "; #Queue options
    #$command .= "-S $opt{MUTECT_SCALA} ";
    #$command .= "-R $opt{GENOME} -O $sample_tumor_name -mem $opt{MUTECT_MEM} -nsc $opt{MUTECT_SCATTER} ";
    #$command .= "-tb $sample_tumor_bam -nb $sample_ref_bam ";
    #$command .= "-D $opt{CALLING_DBSNP} -C $opt{MUTECT_COSMIC} ";
    
    ### Optional settings
    #if ( $opt{SOMVAR_TARGETS} ) {
	#$command .= "-L $opt{SOMVAR_TARGETS} ";
	#if ( $opt{CALLING_INTERVALPADDING} ) {
	    #$command .= "-ip $opt{CALLING_INTERVALPADDING} ";
	#}
    #}
    #if($opt{QUEUE_RETRY} eq 'yes'){
	#$command  .= "-retry 1 ";
    #}
    
    ## Set run option
    #$command .= "-run";
    
    ### Mutect .jar command
    my $command = "java -Xmx".$opt{MUTECT_MEM}."G -jar $opt{MUTECT_PATH}/mutect.jar -T MuTect ";
    $command .= "-R $opt{GENOME} --cosmic $opt{MUTECT_COSMIC} --dbsnp $opt{CALLING_DBSNP} ";
    if ( $opt{SOMVAR_TARGETS} ) {
	$command .= "--intervals $opt{SOMVAR_TARGETS} ";
    }
    $command .= "--input_file:normal $sample_ref_bam --input_file:tumor $sample_tumor_bam ";
    $command .= "--out call_stats.out --vcf $sample_tumor_name\_mutect.vcf";

    ## Create mutect bash script
    my $job_id = "MUT_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";
    open MUTECT_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print MUTECT_SH "#!/bin/bash\n\n";
    print MUTECT_SH "if [ -s $sample_tumor_bam -a -s $sample_ref_bam ]\n";
    print MUTECT_SH "then\n";
    print MUTECT_SH "\techo \"Start Mutect\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/mutect.log\n\n";

    # Run Mutect
    print MUTECT_SH "\tcd $mutect_tmp_dir\n";
    print MUTECT_SH "\t$command\n";
    
    # Filter Mutect result
    $command = "cat $sample_tumor_name\_mutect.vcf | java -Xmx".$opt{MUTECT_MEM}."G -jar $opt{SNPEFF_PATH}/SnpSift.jar filter \"( na FILTER ) | (FILTER = 'PASS')\" > $sample_tumor_name\_mutect_passed.vcf \n";
    print MUTECT_SH "\t$command\n\n";
    # Check Mutect completed
    print MUTECT_SH "\tif [ -s $sample_tumor_name\_mutect.vcf -a -s $sample_tumor_name\_mutect_passed.vcf ]\n";
    print MUTECT_SH "\tthen\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect.vcf $mutect_out_dir/\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect.vcf.idx $mutect_out_dir/\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect_passed.vcf $mutect_out_dir/\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect_passed.vcf.idx $mutect_out_dir/\n";
    print MUTECT_SH "\t\tcd $mutect_out_dir/\n";
    print MUTECT_SH "\t\trm -r tmp/\n";
    print MUTECT_SH "\t\ttouch $log_dir/mutect.done\n";
    print MUTECT_SH "\tfi\n\n";
    print MUTECT_SH "\techo \"End Mutect\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/mutect.log\n\n";

    print MUTECT_SH "else\n";
    print MUTECT_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print MUTECT_SH "fi\n";

    close MUTECT_SH;

    ## Run job
    if ( @running_jobs ){
	#system "qsub -q $opt{MUTECT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_MASTERTHREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
	system "qsub -q $opt{MUTECT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	#system "qsub -q $opt{MUTECT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_MASTERTHREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
	system "qsub -q $opt{MUTECT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }

    return $job_id;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}

sub get_chrs_from_dict {
    my $dictFile = shift;
    #my %chrs;
    my @chrs;
    open DICT, $dictFile;
    while(<DICT>) {
	chomp;
	my ($chr, $length) = ($1, $2) if $_ =~ /SN:(\w+)\s*LN:(\d+)/;
	#$chrs{$chr} = $length if $chr;
	push(@chrs, $chr) if $chr;
    }
    close DICT;
    
    return \@chrs;
}
############

1;
