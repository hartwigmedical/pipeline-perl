package illumina_germlineFiltering;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_sge qw(jobNative qsubJava);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);
use illumina_metadata;


sub runFilterVariants {
    my ($opt) = @_;

    say "\n### SCHEDULING VARIANT FILTRATION ###";

    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-f "$opt->{OUTPUT_DIR}/logs/GermlineFilter.done" || -f "$opt->{OUTPUT_DIR}/logs/VariantFilter.done") {
        say "WARNING: $opt->{OUTPUT_DIR}/logs/GermlineFilter.done exists, skipping";
        return;
    }

    my @snp_types = split ",", $opt->{FILTER_SNPTYPES};
    my @snp_filter_names = split "\t", $opt->{FILTER_SNPNAME};
    my @snp_filter_exprs = split "\t", $opt->{FILTER_SNPEXPR};
    if (scalar @snp_filter_names ne scalar @snp_filter_exprs) {
        die "FILTER_SNPNAME and FILTER_SNPEXPR do not have the same length";
    }

    my @indel_types = split ",", $opt->{FILTER_INDELTYPES};
    my @indel_filter_names = split "\t", $opt->{FILTER_INDELNAME};
    my @indel_filter_exprs = split "\t", $opt->{FILTER_INDELEXPR};
    if (scalar @indel_filter_names ne scalar @indel_filter_exprs) {
        die "FILTER_INDELNAME and FILTER_INDELEXPR do not have the same length";
    }

    my $job_id = "GermlineFilter_" . getJobId();
    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "${job_id}.sh");
    my $log_dir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $stdout = catfile($log_dir, "GermlineFiltering_$opt->{RUN_NAME}.out");
    my $stderr = catfile($log_dir, "GermlineFiltering_$opt->{RUN_NAME}.err");

    from_template("GermlineFiltering.sh.tt", $bash_file,
                  snp_types => \@snp_types,
                  snp_filter_names => \@snp_filter_names,
                  snp_filter_exprs => \@snp_filter_exprs,
                  indel_types => \@indel_types,
                  indel_filter_names => \@indel_filter_names,
                  indel_filter_exprs => \@indel_filter_exprs,
                  job_native => jobNative($opt, "FILTER"),
                  opt => $opt);

    my @running_jobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if (exists $opt->{RUNNING_JOBS}->{$sample} && @{$opt->{RUNNING_JOBS}->{$sample}}) {
            push @running_jobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    my $qsub = qsubJava($opt, "FILTER_MASTER");
    if (@running_jobs) {
        system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
    }

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }

    # dependent GermlineFilter.scala, could fix to be explicit
    my $germline_vcf_path = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.filtered_variants.vcf");
    illumina_metadata::linkArtefact($germline_vcf_path, "germline_vcf", $opt);
    illumina_metadata::linkArtefact("${germline_vcf_path}.idx", "germline_vcf_index", $opt);

    return;
}

1;
