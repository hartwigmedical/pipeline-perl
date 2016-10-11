package illumina_germlineAnnotation;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge qw(getJobId qsubJava);
use illumina_template qw(from_template);


sub runAnnotateVariants {
    my ($opt) = @_;
    my $runName = basename($opt->{OUTPUT_DIR});

    my @runningJobs;
    my $jobID = "GermlineAnnotation_" . getJobId();

    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-f "$opt->{OUTPUT_DIR}/logs/GermlineAnnotation.done" || -f "$opt->{OUTPUT_DIR}/logs/VariantAnnotation.done") {
        say "WARNING: $opt->{OUTPUT_DIR}/logs/GermlineAnnotation.done exists, skipping";
        return;
    }

    my $invcf = "${runName}.filtered_variants.vcf";
    my $preAnnotateVCF = $invcf;
    my $bashFile = catfile($opt->{OUTPUT_DIR}, "jobs", "${jobID}.sh");
    my $logDir = catfile($opt->{OUTPUT_DIR}, "logs");

    from_template("GermlineAnnotation.sh.tt", $bashFile, runName => $runName, invcf => $invcf, preAnnotateVCF => $preAnnotateVCF, opt => $opt);

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if (exists $opt->{RUNNING_JOBS}->{$sample} && @{$opt->{RUNNING_JOBS}->{$sample}}) {
            push @runningJobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    my $qsub = qsubJava($opt, "ANNOTATE");
    if (@runningJobs) {
	    system "$qsub -o $logDir/GermlineAnnotation_$runName.out -e $logDir/GermlineAnnotation_$runName.err -N $jobID -hold_jid " . join(",", @runningJobs) . " $bashFile";
    } else {
	    system "$qsub -o $logDir/GermlineAnnotation_$runName.out -e $logDir/GermlineAnnotation_$runName.err -N $jobID $bashFile";
    }

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $jobID;
    }
    return;
}

1;
