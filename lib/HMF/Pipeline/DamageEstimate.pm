package HMF::Pipeline::DamageEstimate;

use FindBin::libs;
use discipline;

use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    my $dirs = createDirs($opt->{OUTPUT_DIR}, damage_estimate => "damage_estimate");
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    #$dirs->{strelka}->{out} = addSubDir($dirs, "damage_estimate");
    #my $done_file = checkReportedDoneFile("damage_estimate", undef, $dirs, $opt) or return;

    say "\n### SCHEDULING DAMAGE ESTIMATE ###";

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

    #my $job_id = markDone($done_file, [ $ref_job_id, $tumor_job_id ], $dirs, $opt);

    return;
}

1;
