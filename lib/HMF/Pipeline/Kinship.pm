package HMF::Pipeline::Kinship;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING KINSHIP ###";

    my $job_id = "Kinship_" . getId();
    my $log_dir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $done_file = catfile($log_dir, "Kinship.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return;
    }

    my $vcf = "$opt->{RUN_NAME}.filtered_variants.vcf";
    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "${job_id}.sh");
    my $stdout = catfile($log_dir, "Kinship_$opt->{RUN_NAME}.out");
    my $stderr = catfile($log_dir, "Kinship_$opt->{RUN_NAME}.err");

    writeFromTemplate(
        "Kinship.sh.tt", $bash_file,
        vcf => $vcf,
        opt => $opt,
    );

    my @runningJobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if (exists $opt->{RUNNING_JOBS}->{$sample} && @{$opt->{RUNNING_JOBS}->{$sample}}) {
            push @runningJobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    my $qsub = qsubJava($opt, "KINSHIP");
    if (@runningJobs) {
        system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @runningJobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
    }

    $opt->{RUNNING_JOBS}->{kinship} = [$job_id];
    return;
}

1;
