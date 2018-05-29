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
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE CALLING ###";

    runCaller($opt);
    runFiltering($opt);
    runAnnotation($opt);

    HMF::Pipeline::Functions::Metadata::linkVcfArtefacts($opt->{GERMLINE_VCF_FILE}, "germline", $opt);

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

    $opt->{GERMLINE_VCF_FILE} = $calling_vcf;

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

    return;
}

sub runFiltering {
    my ($opt) = @_;

    my (undef, undef, $running_jobs) = refSampleBamAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $filtered_vcf = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.filtered_variants.vcf");
    my $job_id = fromTemplate(
        "GermlineFiltering",
        undef,
        1,
        qsubJava($opt, "GERMLINE_FILTER_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        input_vcf => $opt->{GERMLINE_VCF_FILE},
        final_vcf => $filtered_vcf,
        snp_config => snpConfig($opt),
        indel_config => indelConfig($opt),
        job_native => jobNative($opt, "GERMLINE_FILTER"),
    );

    $opt->{GERMLINE_VCF_FILE} = $filtered_vcf;

    push @{$opt->{RUNNING_JOBS}->{germline}}, [$job_id] if $job_id;

    return;
}

sub runAnnotation {
    my ($opt) = @_;

    my (undef, undef, $running_jobs) = refSampleBamAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $annotated_vcf = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.annotated.vcf");
    my $job_id = fromTemplate(
        "GermlineAnnotation",
        undef,
        1,
        qsubJava($opt, "ANNOTATE"),
        $running_jobs,
        $dirs,
        $opt,
        input_vcf => $opt->{GERMLINE_VCF_FILE},
        final_vcf => $annotated_vcf,
    );

    $opt->{GERMLINE_VCF_FILE} = catfile($annotated_vcf, ".gz");

    push @{$opt->{RUNNING_JOBS}->{germline}}, $job_id if $job_id;

    return;
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
