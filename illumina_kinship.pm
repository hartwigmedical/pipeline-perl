package illumina_kinship;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge qw(qsubJava);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);


sub runKinship {
    my ($opt) = @_;

    my $log_dir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $done_file = catfile($log_dir, "Kinship.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my $vcf = "$opt->{RUN_NAME}.filtered_variants.vcf";
    my $job_id = "Kinship_" . getJobId();
    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "${job_id}.sh");
    my $stdout = catfile($log_dir, "Kinship_$opt->{RUN_NAME}.out");
    my $stderr = catfile($log_dir, "Kinship_$opt->{RUN_NAME}.err");

    from_template("Kinship.sh.tt", $bash_file,
                  vcf => $vcf,
                  opt => $opt);

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
