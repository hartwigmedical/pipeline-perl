package HMF::Pipeline::PreStats;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Config::Validate qw(parseFastqName);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(fromTemplate getId);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "### SCHEDULING PRESTATS ###";

    say "Creating FASTQC report for the following fastq.gz files:";
    foreach my $input_path (keys %{$opt->{FASTQ}}) {
        say "\t$input_path";

        my $fastq = parseFastqName($input_path);
        my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $fastq->{sampleName}), qc => "QCStats");

        my $job_id = fromTemplate(
            "PreStats",
            $fastq->{coreName},
            1,
            qsubTemplate($opt, "PRESTATS"),
            [],
            $dirs,
            $opt,
            sample_name => $fastq->{sampleName},
            core_name => $fastq->{coreName},
            input_path => $input_path,
        );

        push @{$opt->{RUNNING_JOBS}->{prestats}}, $job_id if $job_id;
    }
    return;
}

1;
