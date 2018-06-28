package HMF::Pipeline::Realignment;

use FindBin::libs;
use discipline;

use File::Spec::Functions;
use List::MoreUtils qw(uniq);

use HMF::Pipeline::Functions::Bam;
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Config qw(createDirs);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING INDEL REALIGNMENT ###";
    say "Running single sample indel realignment for the following BAM-files:";

    my $known_files = "";
    $known_files = join " ", map { "-known $_" } split '\t', $opt->{REALIGNMENT_KNOWN} if $opt->{REALIGNMENT_KNOWN};

    my $dirs = createDirs($opt->{OUTPUT_DIR});

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $original_bam_path = $opt->{BAM_FILES}->{$sample};

        HMF::Pipeline::Functions::Bam::operationWithSliceChecks("Realignment", $sample, $known_files, "realigned", "realign", $opt);

        # KODU (DEV-440): Cleanup after realignment needs to wait on poststats.
        my $dependent_jobs = [ uniq $opt->{RUNNING_JOBS}->{$sample}, @{$opt->{RUNNING_JOBS}->{poststats}} ];
        my $realign_cleanup_job = fromTemplate(
            "RealignmentCleanup",
            $sample,
            1,
            qsubTemplate($opt, "REALIGNMENT_MASTER"),
            $dependent_jobs,
            $dirs,
            $opt,
            sample => $sample,
            cleanup_name => "RealignmentCleanup_" . $sample,
            original_bam_path => $original_bam_path,
            # KODU: comment to avoid perltidy putting on one line
        );

        push @{$opt->{RUNNING_JOBS}->{realigncleanup}}, $realign_cleanup_job;

    }
    return;
}

1;