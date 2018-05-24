package HMF::Pipeline::PreStats;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs);
use HMF::Pipeline::Functions::Validate qw(parseFastqName);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Job qw(fromTemplate);

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

        my $fastq_name = fileparse($input_path);
        my $job_id = fromTemplate(
            "PreStats",
            $fastq_name,
            1,
            qsubTemplate($opt, "PRESTATS"),
            [],
            $dirs,
            $opt,
            step => $fastq_name,
            sample_name => $fastq->{sampleName},
            input_path => $input_path,
        );

        push @{$opt->{RUNNING_JOBS}->{prestats}}, $job_id if $job_id;
    }
    return;
}

1;