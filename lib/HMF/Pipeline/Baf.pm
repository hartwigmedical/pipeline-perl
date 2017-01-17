package HMF::Pipeline::Baf;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamAndJobs);
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Metadata qw(linkExtraArtefact);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING BAF ANALYSIS ###";

    my @baf_jobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $sample));
        my ($sample_bam, $running_jobs) = sampleBamAndJobs($sample, $opt);

        my $output_vcf = "${sample}_BAF_SNPS.vcf";
        my $output_baf = "${sample}_BAF.txt";
        my $output_bafplot = "${sample}_BAF.pdf";

        my $job_id = fromTemplate(
            "BAF",
            $sample,
            1,
            qsubJava($opt, "BAF"),
            $running_jobs,
            $dirs,
            $opt,
            sample => $sample,
            sample_bam => $sample_bam,
            output_vcf => $output_vcf,
            output_baf => $output_baf,
            output_bafplot => $output_bafplot,
        );

        foreach my $artefact ($output_vcf, "${output_vcf}.idx", $output_baf, $output_bafplot) {
            linkExtraArtefact(catfile($dirs->{out}, $artefact), $opt);
        }

        push @baf_jobs, $job_id if $job_id;
    }
    $opt->{RUNNING_JOBS}->{'baf'} = \@baf_jobs;
    return;
}

1;
