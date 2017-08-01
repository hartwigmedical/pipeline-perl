package HMF::Pipeline::CallableLoci;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamAndJobs);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Metadata qw(linkExtraArtefact);
use HMF::Pipeline::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING CALLABLE LOCI ANALYSIS ###";

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $sample), mapping => "mapping");
        my ($sample_bam, $running_jobs) = sampleBamAndJobs($sample, $opt);

        my $output_bed = "${sample}_CallableLoci.bed";
        my $output_summary = "${sample}_CallableLoci.txt";

        my $job_id = fromTemplate(
            "CallableLoci",
            $sample,
            1,
            qsubJava($opt, "CALLABLE_LOCI"),
            $running_jobs,
            $dirs,
            $opt,
            sample => $sample,
            sample_bam => $sample_bam,
            output_bed => $output_bed,
            output_summary => $output_summary,
        );
        next unless $job_id;

        foreach my $artefact ($output_bed, $output_summary) {
            linkExtraArtefact(catfile($dirs->{out}, $artefact), $opt);
        }

        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }
    return;
}

1;