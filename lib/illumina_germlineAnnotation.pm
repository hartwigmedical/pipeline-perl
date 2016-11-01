package illumina_germlineAnnotation;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_sge qw(qsubJava);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);
use illumina_metadata;


sub runAnnotateVariants {
    my ($opt) = @_;

    say "\n### SCHEDULING VARIANT ANNOTATION ###";

    my @running_jobs;
    my $job_id = "GermlineAnnotation_" . getJobId();

    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-f "$opt->{OUTPUT_DIR}/logs/GermlineAnnotation.done" || -f "$opt->{OUTPUT_DIR}/logs/VariantAnnotation.done") {
        say "WARNING: $opt->{OUTPUT_DIR}/logs/GermlineAnnotation.done exists, skipping";
        return;
    }

    my $in_vcf = "$opt->{RUN_NAME}.filtered_variants.vcf";
    my $annotated_vcf = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.filtered_variants.annotated.vcf");
    my $pre_annotated_vcf = $in_vcf;
    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "${job_id}.sh");
    my $log_dir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $stdout = catfile($log_dir, "GermlineAnnotation_$opt->{RUN_NAME}.out");
    my $stderr = catfile($log_dir, "GermlineAnnotation_$opt->{RUN_NAME}.err");

    from_template("GermlineAnnotation.sh.tt", $bash_file,
                  in_vcf => $in_vcf,
                  pre_annotated_vcf => $pre_annotated_vcf,
                  annotated_vcf => $annotated_vcf,
                  opt => $opt);

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if (exists $opt->{RUNNING_JOBS}->{$sample} && @{$opt->{RUNNING_JOBS}->{$sample}}) {
            push @running_jobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    my $qsub = qsubJava($opt, "ANNOTATE");
    if (@running_jobs) {
        system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
    }

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }

    illumina_metadata::linkArtefact($annotated_vcf, "germline_vcf", $opt);
    illumina_metadata::linkArtefact("${annotated_vcf}.idx", "germline_vcf_index", $opt);

    return;
}

1;
