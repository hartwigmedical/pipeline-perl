package HMF::Pipeline::Gender;

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

    say "\n### SCHEDULING GENDER ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my (undef, $running_jobs) = sampleBamsAndJobs($opt);
    my $job_id = fromTemplate(
        "Gender",
        undef,
        1,
        qsubJava($opt, "GENDER"),
        $running_jobs,
        $dirs,
        $opt,
        vcf_path => $opt->{GERMLINE_VCF_FILE},
        output_file => catfile($dirs->{out}, "$opt->{RUN_NAME}.gender"),
    );

    push @{$opt->{RUNNING_JOBS}->{gender}}, $job_id;
    return;
}

1;
