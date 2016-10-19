package illumina_germlineAnnotation;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_sge qw(qsubJava);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);


sub runAnnotateVariants {
    my ($opt) = @_;

    my @runningJobs;
    my $jobID = "GermlineAnnotation_" . getJobId();

    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-f "$opt->{OUTPUT_DIR}/logs/GermlineAnnotation.done" || -f "$opt->{OUTPUT_DIR}/logs/VariantAnnotation.done") {
        say "WARNING: $opt->{OUTPUT_DIR}/logs/GermlineAnnotation.done exists, skipping";
        return;
    }

    my $invcf = "$opt->{RUN_NAME}.filtered_variants.vcf";
    my $preAnnotateVCF = $invcf;
    my $bashFile = catfile($opt->{OUTPUT_DIR}, "jobs", "${jobID}.sh");
    my $logDir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $stdout = catfile($logDir, "GermlineAnnotation_$opt->{RUN_NAME}.out");
    my $stderr = catfile($logDir, "GermlineAnnotation_$opt->{RUN_NAME}.err");

    from_template("GermlineAnnotation.sh.tt", $bashFile,
                  invcf => $invcf,
                  preAnnotateVCF => $preAnnotateVCF,
                  opt => $opt);

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if (exists $opt->{RUNNING_JOBS}->{$sample} && @{$opt->{RUNNING_JOBS}->{$sample}}) {
            push @runningJobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    my $qsub = qsubJava($opt, "ANNOTATE");
    if (@runningJobs) {
	    system "$qsub -o $stdout -e $stderr -N $jobID -hold_jid " . join(",", @runningJobs) . " $bashFile";
    } else {
	    system "$qsub -o $stdout -e $stderr -N $jobID $bashFile";
    }

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $jobID;
    }
    return;
}

1;
