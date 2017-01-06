package HMF::Pipeline::GermlineFiltering;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamsAndJobs recordPerSampleJob);
use HMF::Pipeline::Sge qw(jobNative qsubJava);
use HMF::Pipeline::Job qw(getId fromTemplate);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE FILTERING ###";

    my (undef, $running_jobs) = sampleBamsAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $job_id = fromTemplate(
        "GermlineFiltering",
        undef,
        1,
        qsubJava($opt, "FILTER_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        snp_config => snpConfig($opt),
        indel_config => indelConfig($opt),
        job_native => jobNative($opt, "FILTER"),
    );

    # dependent GermlineFilter.scala, could fix to be explicit
    my $germline_vcf_path = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.filtered_variants.vcf");
    HMF::Pipeline::Metadata::linkArtefact($germline_vcf_path, "germline_vcf", $opt);
    HMF::Pipeline::Metadata::linkArtefact("${germline_vcf_path}.idx", "germline_vcf_index", $opt);

    recordPerSampleJob($opt, $job_id) if $job_id;
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
