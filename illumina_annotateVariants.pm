#!/usr/bin/perl -w

package illumina_annotateVariants;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_sge;
use illumina_template;

sub runAnnotateVariants {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "AnnotateVariants_".get_job_id();

    if (-e "$opt{OUTPUT_DIR}/logs/VariantAnnotation.done"){
        print "WARNING: $opt{OUTPUT_DIR}/logs/VariantAnnotation.done exists, skipping \n";
        return $jobID;
    }

    my $invcf = $runName.".filtered_snps.vcf";
    my $preAnnotateVCF = $invcf;
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/AnnotateVariants_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";

    from_template("AnnotateVariants.sh.tt", $bashFile, runName => $runName, invcf => $invcf, preAnnotateVCF => $preAnnotateVCF, opt => \%opt);

    foreach my $sample (@{$opt{SAMPLES}}){
        if( exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}} ) {
            push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
        }
    }

    my $qsub = &qsubJava(\%opt, "ANNOTATE");
    if (@runningJobs) {
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
