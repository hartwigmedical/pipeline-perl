#!/usr/bin/perl -w

#############################################################
### illumina_filterVariants.pm
###
### - Filter variants according to the gatk best practices
### 
### Authors: R.F.Ernst & H.H.D.Kerstens
#############################################################

package illumina_filterVariants;

use strict;
use POSIX qw(tmpnam);
use illumina_sge;
use illumina_template;

sub runFilterVariants {
    ###
    # Run vcf filters
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "FV_".get_job_id();

    ### Skip variant calling if .raw_variants.vcf already exists
    ### add .done file checking?
    if (-e "$opt{OUTPUT_DIR}/logs/VariantFilter.done"){
	print "WARNING: $opt{OUTPUT_DIR}/logs/VariantFilter.done exists, skipping \n";
	return $jobID;
    }

    ### Build Queue command
    my $command = "java -Xmx".$opt{FILTER_MASTER_MEM}."G -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp -jar $opt{QUEUE_PATH}/Queue.jar ";
    my $jobNative = &jobNative(\%opt,"FILTER");
    $command .= "-jobQueue $opt{FILTER_QUEUE} -jobNative \"$jobNative\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/logs/VariantFilter.jobReport.txt ";

    ### Common settings
    $command .= "-S $opt{FILTER_SCALA} -R $opt{GENOME} -V $opt{OUTPUT_DIR}/$runName\.raw_variants.vcf -O $runName -mem $opt{FILTER_MEM} -nsc $opt{FILTER_SCATTER} -mode both ";

	my @SNPFilterNames = split("\t",$opt{FILTER_SNPNAME});
	my @SNPFilterExprs = split("\t",$opt{FILTER_SNPEXPR});
	my @snpTypes = split(",",$opt{FILTER_SNPTYPES});

	## Add SNP Types
	foreach my $snpType (@snpTypes){
		$command.= "-snpType $snpType ";
	}

	## SNP Filter Settings
	if (scalar(@SNPFilterNames) ne scalar(@SNPFilterExprs)) {
		die "FILTER_SNPNAME and FILTER_SNPEXPR do not have the same length.";
	}
	foreach my $i (0 .. scalar(@SNPFilterNames)-1 ){
		$command .= "-snpFilterName $SNPFilterNames[$i] -snpFilterExpression \"$SNPFilterExprs[$i]\" ";
	}
	if ($opt{FILTER_CLUSTERSIZE} and $opt{FILTER_CLUSTERWINDOWSIZE}){
		$command .= "-cluster $opt{FILTER_CLUSTERSIZE} -window $opt{FILTER_CLUSTERWINDOWSIZE} ";
	}
	}

	my @INDELFilterNames = split("\t",$opt{FILTER_INDELNAME});
	my @INDELFilterExprs = split("\t",$opt{FILTER_INDELEXPR});
	my @indelTypes = split(",",$opt{FILTER_INDELTYPES});

	## Add SNP Types
	foreach my $indelType (@indelTypes){
		$command.= "-indelType $indelType ";
	}

	if (scalar(@INDELFilterNames) ne scalar(@INDELFilterExprs)) {
		die "FILTER_INDELNAME and FILTER_INDELEXPR do not have the same length.";
	}

	foreach my $i (0 .. scalar(@INDELFilterNames)-1 ){
		$command .= "-indelFilterName $INDELFilterNames[$i] -indelFilterExpression \"$INDELFilterExprs[$i]\" ";
	}

    ### retry option
    if($opt{QUEUE_RETRY} eq 'yes'){
        $command  .= "-retry 1 ";
    }

    $command .= "-run";

    ### Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/FilterVariants_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";
    from_template("FilterVariants.sh.tt", $bashFile, runName => $runName, command => $command, opt => \%opt);

    ### Process runningjobs
    foreach my $sample (@{$opt{SAMPLES}}){
	if( exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}} ) {
	    push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	}

    ### Start main bash script
    my $qsub = &qsubJava(\%opt,"FILTER_MASTER");
    if (@runningJobs){
	system "$qsub -o $logDir/VariantFilter_$runName.out -e $logDir/VariantFilter_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "$qsub -o $logDir/VariantFilter_$runName.out -e $logDir/VariantFilter_$runName.err -N $jobID $bashFile";
    }

    return $jobID;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
