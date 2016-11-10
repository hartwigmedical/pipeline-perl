package HMF::Pipeline::CallableLoci;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Job qw(getId fromTemplate);
use HMF::Pipeline::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING CALLABLE LOCI ANALYSIS ###";

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $out_dir = catfile($opt->{OUTPUT_DIR}, $sample);
        my $dirs = createDirs($out_dir, mapping => "mapping");
        my $sample_bam = catfile($dirs->{mapping}, $opt->{BAM_FILES}->{$sample});

        my $job_id = fromTemplate(
            "CallableLoci",
            $sample,
            qsubJava($opt, "CALLABLE_LOCI"),
            $opt->{RUNNING_JOBS}->{$sample},
            $dirs,
            $opt,
            sample => $sample,
            sample_bam => $sample_bam,
            output_bed => "${sample}_CallableLoci.bed",
            output_summary => "${sample}_CallableLoci.txt",
        );

        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id if $job_id;
    }
    return;
}

1;
