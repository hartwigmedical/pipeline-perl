package HMF::Pipeline::Amber;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Functions::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Template qw(writeFromTemplate);
use HMF::Pipeline::Functions::Metadata qw(linkArtefact);

use List::Util qw[min max];

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING AMBER ###";
    $opt->{RUNNING_JOBS}->{'amber'} = [];

    my $sub_dir = "amber";
    my $dirs = createDirs($opt->{OUTPUT_DIR}, amber => $sub_dir);
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $done_file = checkReportedDoneFile("Amber_$joint_name", undef, $dirs, $opt) or return;

    my @amber_jobs;
    my $ref_threads = max(1, int($opt->{AMBER_THREADS} * 1 / 4));
    my $tumor_threads = max(1, $opt->{AMBER_THREADS} - $ref_threads);
    say "Using $ref_threads threads for reference pileup and $tumor_threads threads for tumor";

    push @amber_jobs, runAmberPileup($ref_sample, $ref_bam_path, $ref_threads, $running_jobs, $dirs, $opt);
    push @amber_jobs, runAmberPileup($tumor_sample, $tumor_bam_path, $tumor_threads, $running_jobs, $dirs, $opt);
    push @amber_jobs, runAmber($tumor_sample, $ref_bam_path, $tumor_bam_path, \@amber_jobs, $dirs, $opt);
    push @amber_jobs, markDone($done_file, \@amber_jobs, $dirs, $opt);
    push @{$opt->{RUNNING_JOBS}->{'amber'}}, @amber_jobs;

    my $baf_txt_path = "${sub_dir}/${tumor_sample}.amber.baf";
    linkArtefact($baf_txt_path, 'somatic_amber_baf', $opt);

    return;
}

sub runAmber {
    my ($tumor_sample, $ref_bam_path, $tumor_bam_path, $running_jobs, $dirs, $opt) = @_;

    say "\n### SCHEDULING AMBER ###";
    my $job_id = fromTemplate(
        "Amber",
        undef,
        1,
        qsubTemplate($opt, "AMBER"),
        $running_jobs,
        $dirs,
        $opt,
        tumor_sample => $tumor_sample,
        ref_bam_path => $ref_bam_path,
        tumor_bam_path => $tumor_bam_path,
    );

    return $job_id;
}

sub runAmberPileup {
    my ($sample, $sample_bam, $threads, $running_jobs, $dirs, $opt) = @_;

    say "\n### SCHEDULING AMBER PILEUP ON $sample ###";
    my $job_id = fromTemplate(
        "AmberPileup",
        $sample,
        1,
        qsubTemplate($opt, "AMBER"),
        $running_jobs,
        $dirs,
        $opt,
        sample => $sample,
        sample_bam => $sample_bam,
        threads => $threads,
    );

    return $job_id;
}

1;
