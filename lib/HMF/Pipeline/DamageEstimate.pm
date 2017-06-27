package HMF::Pipeline::DamageEstimate;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "damageEstimate", $joint_name));
    my $done_file = checkReportedDoneFile($joint_name, undef, $dirs, $opt) or return;

    say "\n### SCHEDULING DAMAGE ESTIMATE ###";

    my @job_ids;
    my $ref_job_id = fromTemplate(
        "DamageEstimate",
        undef,
        1,
        qsubJava($opt, "DAMAGE_ESTIMATE"),
        $running_jobs,
        $dirs,
        $opt,
        damage_estimate_bam_path => $ref_bam_path,
        joint_name => $joint_name
    );
    push @job_ids, $ref_job_id;

    my $tumor_job_id = fromTemplate(
        "DamageEstimate",
        undef,
        1,
        qsubJava($opt, "DAMAGE_ESTIMATE"),
        $running_jobs,
        $dirs,
        $opt,
        damage_estimate_bam_path => $tumor_bam_path,
        joint_name => $joint_name
    );
    push @job_ids, $tumor_job_id;

    my $job_id = markDone($done_file, [@job_ids], $dirs, $opt);

    return;
}

1;
