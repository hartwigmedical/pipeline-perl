package HMF::Pipeline::GermlineCalling;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs refSampleBamAndJobs recordAllSampleJob);
use HMF::Pipeline::Functions::Sge qw(jobNative qsubJava);
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Metadata qw(linkArtefact linkVcfArtefacts);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE CALLING ###";

    my ($ref_sample, $ref_sample_bam, $running_jobs) = refSampleBamAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR}, gvcf => "gvcf");

    my $final_vcf = catfile($dirs->{out}, "$opt->{RUN_NAME}.raw_variants.vcf");
    my $final_gvcf = catfile($dirs->{gvcf}, $ref_sample . ".g.vcf.gz");
    my $ref_sample_bam_name = fileparse($ref_sample_bam);
    (my $tmp_scala_gvcf = $ref_sample_bam_name) =~ s/\.bam/\.g\.vcf\.gz/;

    $opt->{GERMLINE_VCF_FILE} = $final_vcf;
    $opt->{GERMLINE_GVCF_FILE} = $final_gvcf;

    my $job_id = fromTemplate(
        "GermlineCalling",
        undef,
        1,
        qsubJava($opt, "GERMLINE_CALLING_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        ref_sample_bam => $ref_sample_bam,
        final_vcf => $final_vcf,
        final_gvcf => $final_gvcf,
        tmp_scala_gvcf => $tmp_scala_gvcf,
        job_native => jobNative($opt, "GERMLINE_CALLING"),
    );
    return unless $job_id;

    linkVcfArtefacts($final_vcf, "germline", $opt);
    linkArtefact($final_gvcf, 'ref_sample_gvcf', $opt);
    linkArtefact($final_gvcf . '.tbi', 'ref_sample_gvcf_index', $opt);

    recordAllSampleJob($opt, $job_id);
    return;
}

1;
