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
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = basename($opt{OUTPUT_DIR});
    my $jobID = "Kinship_" . getJobId();

    my $done_file = catfile($opt{OUTPUT_DIR}, "logs", "Kinship.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return $jobID;
    }

    my $vcf = "${runName}.filtered_variants.vcf";
    my $bashFile = catfile($opt{OUTPUT_DIR}, "jobs", "${jobID}.sh");
    my $logDir = catfile($opt{OUTPUT_DIR}, "logs");

    from_template("Kinship.sh.tt", $bashFile,
                  vcf => $vcf,
                  opt => \%opt,
                  runName => $runName);

    my @runningJobs;
    foreach my $sample (keys $opt{SAMPLES}) {
        if (exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}}) {
            push @runningJobs, join(",", @{$opt{RUNNING_JOBS}->{$sample}});
        }
    }

    my $qsub = qsubJava(\%opt, "KINSHIP");
    if (@runningJobs) {
        system "$qsub -o $logDir/Kinship_$runName.out -e $logDir/Kinship_$runName.err -N $jobID -hold_jid " . join(",", @runningJobs) . " $bashFile";
    } else {
        system "$qsub -o $logDir/Kinship_$runName.out -e $logDir/Kinship_$runName.err -N $jobID $bashFile";
    }

    return [$jobID];
}

1;
