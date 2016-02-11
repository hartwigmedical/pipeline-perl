#!/usr/bin/perl -w

##################################################################
### illumina_calling.pm
### - Run gatk variant callers, depending on qscript
###  - Haplotype caller: normal and gvcf mode
###  - Unified genotyper 
### - VCF Prep function if pipeline is started with a vcf file
###
### Author: R.F.Ernst
##################################################################

package illumina_calling;

use strict;
use POSIX qw(tmpnam);

sub runVariantCalling {
    ###
    # Run variant callers
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @sampleBams;
    my @runningJobs;
    my $jobID = "VC_".get_job_id();

    ### Skip variant calling if .raw_variants.vcf already exists
    if (-e "$opt{OUTPUT_DIR}/logs/VariantCaller.done"){
	print "WARNING: $opt{OUTPUT_DIR}/logs/VariantCaller.done exists, skipping \n";
	return \%opt;
    }
    
    ### Create gvcf folder if CALLING_GVCF eq yes
    if((! -e "$opt{OUTPUT_DIR}/gvcf" && $opt{CALLING_GVCF} eq 'yes')){
	mkdir("$opt{OUTPUT_DIR}/gvcf") or die "Couldn't create directory: $opt{OUTPUT_DIR}/gvcf\n";
    }
    
    ### Build Queue command
    my $command = "java -Xmx".$opt{CALLING_MASTERMEM}."G -jar $opt{QUEUE_PATH}/Queue.jar ";
    $command .= "-jobQueue $opt{CALLING_QUEUE} -jobNative \"-pe threaded $opt{CALLING_THREADS} -P $opt{CLUSTER_PROJECT}\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/logs/VariantCaller.jobReport.txt "; #Queue options

    ### Add caller and UG specific settings
    $command .= "-S $opt{CALLING_SCALA} ";
    if ($opt{CALLING_UGMODE}) {
	$command .= " -glm $opt{CALLING_UGMODE} ";
    }

    ### Common settings
    $command .= "-R $opt{GENOME} -O $runName -mem $opt{CALLING_MEM} -nct $opt{CALLING_THREADS} -nsc $opt{CALLING_SCATTER} -stand_call_conf $opt{CALLING_STANDCALLCONF} -stand_emit_conf $opt{CALLING_STANDEMITCONF} ";

    ### Add all bams
    foreach my $sample (@{$opt{SAMPLES}}){
	my $sampleBam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";

	$command .= "-I $sampleBam ";
	push( @sampleBams, $sampleBam);
	## Running jobs
	if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
	}
    }

    ### Optional settings
    if ( $opt{CALLING_DBSNP} ) {
	$command .= "-D $opt{CALLING_DBSNP} ";
    }
    if ( $opt{CALLING_TARGETS} ) {
	$command .= "-L $opt{CALLING_TARGETS} ";
	if ( $opt{CALLING_INTERVALPADDING} ) {
	    $command .= "-ip $opt{CALLING_INTERVALPADDING} ";
	}
    }
    if ( $opt{CALLING_PLOIDY} ) {
	$command .= "-ploidy $opt{CALLING_PLOIDY} ";
    }
    
    ### retry option
    if($opt{QUEUE_RETRY} eq 'yes'){
	$command  .= "-retry 1 ";
    }
    $command .= "-run";

    #Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/VariantCalling_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";

    open CALLING_SH, ">$bashFile" or die "cannot open file $bashFile \n";
    print CALLING_SH "#!/bin/bash\n\n";
    print CALLING_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
    print CALLING_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
    print CALLING_SH "echo \"Start variant caller\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n\n";
    
    print CALLING_SH "if [ -s ".shift(@sampleBams)." ";
    foreach my $sampleBam (@sampleBams){
	print CALLING_SH "-a -s $sampleBam ";
    }
    print CALLING_SH "]\n";
    print CALLING_SH "then\n";
    print CALLING_SH "\t$command\n";
    print CALLING_SH "else\n";
    print CALLING_SH "\techo \"ERROR: One or more input bam files do not exist.\" >&2\n";
    print CALLING_SH "fi\n\n";
    
    print CALLING_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.raw_variants.vcf.done ]\n";
    print CALLING_SH "then\n";
    print CALLING_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.raw_variants.vcf $opt{OUTPUT_DIR}/\n";
    print CALLING_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.raw_variants.vcf.idx $opt{OUTPUT_DIR}/\n";
    if($opt{CALLING_GVCF} eq 'yes'){
	print CALLING_SH "\tmv $opt{OUTPUT_DIR}/tmp/*.g.vcf.gz $opt{OUTPUT_DIR}/gvcf/\n";
	print CALLING_SH "\tmv $opt{OUTPUT_DIR}/tmp/*.g.vcf.gz.tbi $opt{OUTPUT_DIR}/gvcf/\n";
	
    }
    print CALLING_SH "\ttouch $opt{OUTPUT_DIR}/logs/VariantCaller.done\n";
    print CALLING_SH "fi\n\n";
    print CALLING_SH "echo \"Finished variant caller\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
    
    #Start main bash script
    if (@runningJobs){
	system "qsub -q $opt{CALLING_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CALLING_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/VariantCaller_$runName.out -e $logDir/VariantCaller_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{CALLING_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CALLING_MASTERTHREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/VariantCaller_$runName.out -e $logDir/VariantCaller_$runName.err -N $jobID $bashFile";
    }
    
    ### Store jobID
    foreach my $sample (@{$opt{SAMPLES}}){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $jobID);
    }
    return \%opt;
}

sub runVcfPrep {
    ###
    # Run vcf prep when starting pipeline with a vcf file.
    ##
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];

    symlink($opt{VCF},"$opt{OUTPUT_DIR}/$runName.raw_variants.vcf");
    @{$opt{SAMPLES}} = ($runName);
    
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