package HMF::Pipeline::Purple;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamsAndJobs);
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Job qw(fromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING PURPLE ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR}, purple => "purple");
    my $dependent_jobs = copyNumberAndGermlineJobs($opt);
    my $job_id = fromTemplate(
        "Purple",
        undef,
        1,
        qsubJava($opt, "PURPLE"),
        $dependent_jobs,
        $dirs,
        $opt,
        # comment to avoid perltidy putting on one line
    );

    $opt->{RUNNING_JOBS}->{purple} = [$job_id] if $job_id;
    return;
}

sub copyNumberAndGermlineJobs {
    my ($opt) = @_;

    my (undef, $running_jobs) = sampleBamsAndJobs($opt);
    my @jobs;
    push @jobs, @{$opt->{RUNNING_JOBS}->{cnv}};
    push @jobs, @{$running_jobs};
    return \@jobs;
}


1;
