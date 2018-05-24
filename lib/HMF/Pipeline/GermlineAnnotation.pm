package HMF::Pipeline::GermlineAnnotation;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs sampleBamsAndJobs recordAllSampleJob);
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Sge qw(qsubJava);
use HMF::Pipeline::Functions::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE ANNOTATION ###";

    my (undef, $running_jobs) = sampleBamsAndJobs($opt);
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

    $opt->{GERMLINE_VCF_FILE} = $annotated_vcf;
    push @{$opt->{RUNNING_JOBS}->{germline}}, $job_id;
    return unless $job_id;

    HMF::Pipeline::Metadata::linkVcfArtefacts($opt->{GERMLINE_VCF_FILE}, "germline", $opt);
    recordAllSampleJob($opt, $job_id);
    return;
}

1;
