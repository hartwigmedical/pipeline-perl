package HMF::Pipeline::Cobalt;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING COBALT ANALYSIS ###";

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR}, cobalt => "cobalt");

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

    push @{$opt->{RUNNING_JOBS}->{'cobalt'}}, $job_id;
    return $job_id;
}

1;