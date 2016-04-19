#!/usr/bin/perl -w

##################################################################################################################################################
### illumina_annotateVariants.pm
### - Annotate vcf file using multiple tools:
###   - SnpEff
###   - SnpSift -> DBNSFP
###   - ID from vcf file, for example Cosmic
###   - AC and AF from a vcf file, for example GoNL
### Authors: R.F.Ernst & H.H.D.Kerstens
##################################################################################################################################################

package illumina_annotateVariants;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_sge;
use illumina_template;

sub runAnnotateVariants {
    ###
    # Run annotation tools
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $command;
    my $jobID = "AV_".get_job_id();

    ### Skip variant annotation if .done file exists.
    if (-e "$opt{OUTPUT_DIR}/logs/VariantAnnotation.done"){
	print "WARNING: $opt{OUTPUT_DIR}/logs/VariantAnnotation.done exists, skipping \n";
	return $jobID;
    }

    ### vcf file
    my $invcf;
    my $outvcf;
    if ( $opt{FILTER_VARIANTS} eq "yes" ) {
	if ( $opt{FILTER_MODE} eq "BOTH" ) { $invcf = $runName.".filtered_variants.vcf"; }
	if ( $opt{FILTER_MODE} eq "SNP" ) { $invcf = $runName.".filtered_snps.vcf"; }
	if ( $opt{FILTER_MODE} eq "INDEL" ) { $invcf = $runName.".filtered_indels.vcf"; }
    } elsif ($opt{FILTER_VARIANTS} eq "no") { $invcf = $runName.".raw_variants.vcf"; }
    my $preAnnotateVCF = $invcf;

    ### Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/AnnotateVariants_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";
    from_template("AnnotateVariants.sh.tt", $bashFile, runName => $runName, invcf => $invcf, preAnnotateVCF => $preAnnotateVCF, opt => \%opt);

    ### Process runningjobs
    foreach my $sample (@{$opt{SAMPLES}}){
	if( exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}} ) {
	    push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	}
    }

    ### Start main bash script
    my $qsub = &qsubJava(\%opt,"ANNOTATE");
    if (@runningJobs){
	system "$qsub -o $logDir/VariantAnnotation_$runName.out -e $logDir/VariantAnnotation_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "$qsub -o $logDir/VariantAnnotation_$runName.out -e $logDir/VariantAnnotation_$runName.err -N $jobID $bashFile";
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
