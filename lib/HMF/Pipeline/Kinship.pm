package HMF::Pipeline::Kinship;

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

    say "\n### SCHEDULING KINSHIP ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my (undef, $running_jobs) = sampleBamsAndJobs($opt);
    my $job_id = fromTemplate(
        "Kinship",
        undef,
        1,
        qsubJava($opt, "KINSHIP"),
        $running_jobs,
        $dirs,
        $opt,
        vcf_path => $opt->{GERMLINE_VCF_FILE},
        # comment to avoid perltidy putting on one line
    );

    $opt->{RUNNING_JOBS}->{kinship} = [$job_id] if $job_id;
    return;
}

1;
