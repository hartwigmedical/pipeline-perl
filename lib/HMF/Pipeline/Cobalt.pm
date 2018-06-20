package HMF::Pipeline::Cobalt;

use FindBin::libs;
use discipline;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Functions::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Functions::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING COBALT ###";
    $opt->{RUNNING_JOBS}->{cobalt} = [];

    my $sub_dir = "cobalt";
    my $dirs = createDirs($opt->{OUTPUT_DIR}, cobalt => $sub_dir);
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $done_file = checkReportedDoneFile("Cobalt_$joint_name", undef, $dirs, $opt) or return;

    my @cobalt_jobs;
    push @cobalt_jobs, runCobalt($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $running_jobs, $dirs, $opt);
    push @cobalt_jobs, markDone($done_file, \@cobalt_jobs, $dirs, $opt);
    push @{$opt->{RUNNING_JOBS}->{cobalt}}, @cobalt_jobs;

    return;
}

sub runCobalt {
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $running_jobs, $dirs, $opt) = @_;

    my $job_id = fromTemplate(
        "Cobalt",
        undef,
        1,
        qsubJava($opt, "COBALT"),
        $running_jobs,
        $dirs,
        $opt,
        ref_sample => $ref_sample,
        ref_bam_path => $ref_bam_path,
        tumor_sample => $tumor_sample,
        tumor_bam_path => $tumor_bam_path,
    );

    return $job_id;
}

1;
