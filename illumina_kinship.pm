package illumina_kinship;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge;
use illumina_template;


sub runKinship {
    my ($opt) = @_;
    my $run_name = basename($opt->{OUTPUT_DIR});

    my $log_dir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $done_file = catfile($log_dir, "Kinship.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my $vcf = "${run_name}.filtered_variants.vcf";
    my $job_id = "Kinship_" . getJobId();
    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "${job_id}.sh");

    from_template("Kinship.sh.tt", $bash_file,
                  vcf => $vcf,
                  opt => $opt,
                  run_name => $run_name);

    my @runningJobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if (exists $opt->{RUNNING_JOBS}->{$sample} && @{$opt->{RUNNING_JOBS}->{$sample}}) {
            push @runningJobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    my $qsub = qsubJava($opt, "KINSHIP");
    if (@runningJobs) {
        system "$qsub -o $log_dir/Kinship_${run_name}.out -e $log_dir/Kinship_${run_name}.err -N $job_id -hold_jid " . join(",", @runningJobs) . " $bash_file";
    } else {
        system "$qsub -o $log_dir/Kinship_${run_name}.out -e $log_dir/Kinship_${run_name}.err -N $job_id $bash_file";
    }

    $opt->{RUNNING_JOBS}->{kinship} = [$job_id];
}

1;
