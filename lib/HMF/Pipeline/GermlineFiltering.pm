package HMF::Pipeline::GermlineFiltering;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamsAndJobs recordAllSampleJob);
use HMF::Pipeline::Sge qw(jobNative qsubJava);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Metadata;

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
        qsubJava($opt, "FILTER_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        input_vcf => $opt->{GERMLINE_VCF_FILE},
        final_vcf => $germline_vcf_path,
        snp_config => snpConfig($opt),
        indel_config => indelConfig($opt),
        job_native => jobNative($opt, "FILTER"),
    );

    $opt->{GERMLINE_VCF_FILE} = $germline_vcf_path;
    return unless $job_id;

    HMF::Pipeline::Metadata::linkVcfArtefacts($opt->{GERMLINE_VCF_FILE}, "germline", $opt);
    recordAllSampleJob($opt, $job_id);
    return;
}

# NB: these functions have not been made more general in order to maintain easy search of config key names
sub snpConfig {
    my ($opt) = @_;

    my @snp_filter_names = split "\t", $opt->{FILTER_SNPNAME};
    my @snp_filter_exprs = split "\t", $opt->{FILTER_SNPEXPR};
    if (scalar @snp_filter_names ne scalar @snp_filter_exprs) {
        die "FILTER_SNPNAME and FILTER_SNPEXPR do not have the same length";
    }

    my %config;
    $config{types} = [ split ",", $opt->{FILTER_SNPTYPES} ];
    @{$config{filters}}{@snp_filter_names} = @snp_filter_exprs;
    return \%config;
}

sub indelConfig {
    my ($opt) = @_;

    my @indel_filter_names = split "\t", $opt->{FILTER_INDELNAME};
    my @indel_filter_exprs = split "\t", $opt->{FILTER_INDELEXPR};
    if (scalar @indel_filter_names ne scalar @indel_filter_exprs) {
        die "FILTER_INDELNAME and FILTER_INDELEXPR do not have the same length";
    }

    my %config;
    $config{types} = [ split ",", $opt->{FILTER_INDELTYPES} ];
    @{$config{filters}}{@indel_filter_names} = @indel_filter_exprs;
    return \%config;
}

1;