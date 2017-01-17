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

    my $job_id = fromTemplate(
        "GermlineCalling",
        undef,
        1,
        qsubJava($opt, "CALLING_MASTER"),
        $running_jobs,
        $dirs,
        $opt,
        sample_bams => [ values %{$sample_bams} ],
        final_vcf => $opt->{GERMLINE_VCF_FILE},
        job_native => jobNative($opt, "CALLING"),
    );

    linkArtefacts($dirs->{gvcf}, $opt);
    recordAllSampleJob($opt, $job_id) if $job_id;
    return;
}

# naming dependent on GermlineCaller.scala, could fix to be explicit
sub linkArtefacts {
    my ($gvcf_dir, $opt) = @_;

    if ($opt->{CALLING_GVCF} eq "yes") {
        foreach my $sample (keys %{$opt->{SAMPLES}}) {
            my $bam_file = $opt->{BAM_FILES}->{$sample};
            (my $gvcf_file = $bam_file) =~ s/\.bam$/.g.vcf.gz/;
            my $gvcf_path = catfile($gvcf_dir, $gvcf_file);
            my $sample_name = HMF::Pipeline::Metadata::metaSampleName($sample, $opt);
            HMF::Pipeline::Metadata::linkArtefact($gvcf_path, "${sample_name}_gvcf", $opt);
            HMF::Pipeline::Metadata::linkArtefact("${gvcf_path}.tbi", "${sample_name}_gvcf_index", $opt);
        }
    }
    HMF::Pipeline::Metadata::linkVcfArtefacts($opt->{GERMLINE_VCF_FILE}, "germline", $opt);
    return;
}

1;
