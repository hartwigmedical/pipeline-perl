package HMF::Pipeline::GermlineCalling;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs refSampleBamAndJobs);
use HMF::Pipeline::Functions::Sge qw(jobNative qsubJava);
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Metadata qw(linkVcfArtefacts);

use parent qw(Exporter);
our @EXPORT_OK = qw(run runGermlineRerun);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE CALLING ###";

    my ($caller_job_id, $caller_vcf) = runCaller($opt);
    my ($filter_job_id, $filtered_vcf) = runFiltering($opt, $caller_job_id, $caller_vcf);
    my $annotated_vcf = runAnnotation($opt, $filter_job_id, $filtered_vcf);

    HMF::Pipeline::Functions::Metadata::linkVcfArtefacts($annotated_vcf, "germline_variant", $opt);

    return;
}

sub runGermlineRerun {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE RERUN PROCESSING ###";

    my (undef, undef, $running_jobs) = refSampleBamAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $final_vcf = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.annotated.vcf.gz");

    my $job_id = fromTemplate(
        "GermlineRerunProcess",
        undef,
        1,
        qsubJava($opt, "GERMLINE_CALLING"),
        $running_jobs,
        $dirs,
        $opt,
        final_vcf => $final_vcf,
        # SABR: comment to avoid perltidy putting on one line
    );

    push @{$opt->{RUNNING_JOBS}->{germline}}, [$job_id] if $job_id;

    HMF::Pipeline::Functions::Metadata::linkVcfArtefacts($final_vcf, "germline_variant", $opt);
    return;
}


sub runCaller {
    my ($opt) = @_;

    my ($ref_sample, $ref_sample_bam, $running_jobs) = refSampleBamAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR}, gvcf => "gvcf");

    my $calling_vcf = catfile($dirs->{out}, "$opt->{RUN_NAME}.raw_variants.vcf");
    my $final_gvcf = catfile($dirs->{gvcf}, $ref_sample . ".g.vcf.gz");
    my $ref_sample_bam_name = fileparse($ref_sample_bam);
    (my $tmp_scala_gvcf = $ref_sample_bam_name) =~ s/\.bam/\.g\.vcf\.gz/;

    my $job_id = fromTemplate(
        "GermlineCalling",
        undef,
        1,
        qsubJava($opt, "GERMLINE_CALLING_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        ref_sample_bam => $ref_sample_bam,
        final_vcf => $calling_vcf,
        final_gvcf => $final_gvcf,
        tmp_scala_gvcf => $tmp_scala_gvcf,
        job_native => jobNative($opt, "GERMLINE_CALLING"),
    );

    push @{$opt->{RUNNING_JOBS}->{germline}}, [$job_id] if $job_id;

    return ($job_id, $calling_vcf);
}

sub runFiltering {
    my ($opt, $caller_job_id, $input_vcf) = @_;

    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $filtered_vcf = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.filtered_variants.vcf");
    my $job_id = fromTemplate(
        "GermlineFiltering",
        undef,
        1,
        qsubJava($opt, "GERMLINE_FILTER_MASTER"),
        [$caller_job_id],
        $dirs,
        $opt,
        input_vcf => $input_vcf,
        final_vcf => $filtered_vcf,
        snp_config => snpConfig($opt),
        indel_config => indelConfig($opt),
        job_native => jobNative($opt, "GERMLINE_FILTER"),
    );

    push @{$opt->{RUNNING_JOBS}->{germline}}, [$job_id] if $job_id;

    return ($job_id, $filtered_vcf);
}

sub runAnnotation {
    my ($opt, $filter_job_id, $input_vcf) = @_;

    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $annotated_vcf = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.annotated.vcf");
    my $job_id = fromTemplate(
        "GermlineAnnotation",
        undef,
        1,
        qsubJava($opt, "GERMLINE_ANNOTATE"),
        [$filter_job_id],
        $dirs,
        $opt,
        input_vcf => $input_vcf,
        final_vcf => $annotated_vcf,
    );

    push @{$opt->{RUNNING_JOBS}->{germline}}, $job_id if $job_id;

    # KODU: We implicitly gzip the final vcf in the annotation bash script.
    return join "", $annotated_vcf, ".gz";
}

# NB: SABR: these functions have not been made more general in order to maintain easy search of config key names
sub snpConfig {
    my ($opt) = @_;

    my @snp_filter_names = split "\t", $opt->{GERMLINE_FILTER_SNPNAME};
    my @snp_filter_exprs = split "\t", $opt->{GERMLINE_FILTER_SNPEXPR};
    if (scalar @snp_filter_names ne scalar @snp_filter_exprs) {
        die "GERMLINE_FILTER_SNPNAME and GERMLINE_FILTER_SNPEXPR do not have the same length";
    }

    my %config;
    $config{types} = [ split ",", $opt->{GERMLINE_FILTER_SNPTYPES} ];
    @{$config{filters}}{@snp_filter_names} = @snp_filter_exprs;
    return \%config;
}

sub indelConfig {
    my ($opt) = @_;

    my @indel_filter_names = split "\t", $opt->{GERMLINE_FILTER_INDELNAME};
    my @indel_filter_exprs = split "\t", $opt->{GERMLINE_FILTER_INDELEXPR};
    if (scalar @indel_filter_names ne scalar @indel_filter_exprs) {
        die "GERMLINE_FILTER_INDELNAME and GERMLINE_FILTER_INDELEXPR do not have the same length";
    }

    my %config;
    $config{types} = [ split ",", $opt->{GERMLINE_FILTER_INDELTYPES} ];
    @{$config{filters}}{@indel_filter_names} = @indel_filter_exprs;
    return \%config;
}

1;
