package HMF::Pipeline::Prestats;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "### SCHEDULING PRESTATS ###";

    say "Creating FASTQC report for the following fastq.gz files:";
    foreach my $input_file (keys %{$opt->{FASTQ}}) {
        my $core_name = fileparse($input_file);
        $core_name =~ s/\.fastq.gz$//;
        my ($sample_name) = split "_", $core_name;
        say "\t$input_file";

        my $out_dir = catfile($opt->{OUTPUT_DIR}, $sample_name);
        my $dirs = createDirs($out_dir, qc => "QCStats");

        my $done_file = catfile($dirs->{log}, "PreStats_${sample_name}.done");
        if (!-f $done_file) {
            my $job_id = "PreStats_${core_name}_" . getId();
            my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

            writeFromTemplate(
                "PreStats.sh.tt", $bash_file,
                sample_name => $sample_name,
                core_name => $core_name,
                input => $input_file,
                dirs => $dirs,
                opt => $opt,
            );

            my $qsub = qsubTemplate($opt, "PRESTATS");
            system("$qsub -o $dirs->{log}/PreStat_${core_name}.out -e $dirs->{log}/PreStats_${core_name}.err -N $job_id $bash_file");
            push @{$opt->{RUNNING_JOBS}->{prestats}}, $job_id;
        } else {
            say "\tWARNING: $done_file exists, skipping";
        }
    }
    return;
}

1;
