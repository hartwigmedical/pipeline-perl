package HMF::Pipeline::GermlineFiltering;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs sampleBamsAndJobs recordAllSampleJob);
use HMF::Pipeline::Functions::Sge qw(jobNative qsubJava);
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE FILTERING ###";

    my (undef, $running_jobs) = sampleBamsAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $germline_vcf_path = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.filtered_variants.vcf");
    my $job_id = fromTemplate(
        "GermlineFiltering",
        undef,
        1,
        qsubJava($opt, "GERMLINE_FILTER_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        input_vcf => $opt->{GERMLINE_VCF_FILE},
        final_vcf => $germline_vcf_path,
        snp_config => snpConfig($opt),
        indel_config => indelConfig($opt),
        job_native => jobNative($opt, "GERMLINE_FILTER"),
    );

    $opt->{GERMLINE_VCF_FILE} = $germline_vcf_path;
    return unless $job_id;

    HMF::Pipeline::Functions::Metadata::linkVcfArtefacts($opt->{GERMLINE_VCF_FILE}, "germline", $opt);
    recordAllSampleJob($opt, $job_id);
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