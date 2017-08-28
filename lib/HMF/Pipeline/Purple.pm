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
    my $dependent_jobs = dependencies($opt);
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

sub dependencies {
    my ($opt) = @_;

    my @jobs;
    push @jobs, @{$opt->{RUNNING_JOBS}->{amber}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{cobalt}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{sv}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{somvar}};
    return \@jobs;
}

1;