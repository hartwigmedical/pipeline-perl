package HMF::Pipeline::GermlineCalling;

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

    say "\n### SCHEDULING GERMLINE CALLING ###";

    my ($sample_bams, $running_jobs) = sampleBamsAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR}, gvcf => "gvcf");
    $opt->{GERMLINE_VCF_FILE} = catfile($dirs->{out}, "$opt->{RUN_NAME}.raw_variants.vcf");
    $opt->{GVCF_FILES} = {map { $_ => catfile($dirs->{gvcf}, "${_}.g.vcf.gz") } keys %{$opt->{SAMPLES}}} if $opt->{CALLING_GVCF} eq "yes";

    my $job_id = fromTemplate(
        "GermlineCalling",
        undef,
        1,
        qsubJava($opt, "CALLING_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        sample_bams => $sample_bams,
        final_vcf => $opt->{GERMLINE_VCF_FILE},
        final_gvcfs => $opt->{GVCF_FILES},
        job_native => jobNative($opt, "CALLING"),
    );
    return unless $job_id;

    linkArtefacts($opt);
    recordAllSampleJob($opt, $job_id);
    return;
}

sub linkArtefacts {
    my ($opt) = @_;

    if ($opt->{CALLING_GVCF} eq "yes") {
        while (my ($sample, $gvcf_path) = each %{$opt->{GVCF_FILES}}) {
            my $meta_sample_name = HMF::Pipeline::Metadata::metaSampleName($sample, $opt);
            HMF::Pipeline::Metadata::linkArtefact($gvcf_path, "${meta_sample_name}_gvcf", $opt);
            HMF::Pipeline::Metadata::linkArtefact("${gvcf_path}.tbi", "${meta_sample_name}_gvcf_index", $opt);
        }
    }
    HMF::Pipeline::Metadata::linkVcfArtefacts($opt->{GERMLINE_VCF_FILE}, "germline", $opt);
    return;
}

1;