package illumina_finalize;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge;
use illumina_template;
use illumina_metadataParser;


sub runFinalize {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $run_name = basename($opt{OUTPUT_DIR});

    my $job_id = "${run_name}_" . getJobId();
    my $bash_file = catfile($opt{OUTPUT_DIR}, "jobs", "Finalize_${job_id}.sh");
    my $log_file = catfile($opt{OUTPUT_DIR}, "logs", "PipelineCheck.log");

    my $joint_name = "";
    if ($opt{SOMATIC_VARIANTS} eq "yes" || ($opt{COPY_NUMBER} eq "yes" && $opt{CNV_MODE} eq "sample_control")) {
        my $metadata = metadataParse($opt{OUTPUT_DIR});
        my $ref_sample = $metadata->{ref_sample};
        my $tumor_sample = $metadata->{tumor_sample};
        $joint_name = "${ref_sample}_${tumor_sample}";
    }

    my @runningJobs = map { @$_ }
        grep { defined } @{$opt{RUNNING_JOBS}}{"baf",
                                               keys %{$opt{SAMPLES}},
                                               "postStats",
                                               "somVar",
                                               "CNV",
                                               "Kinship"
                                           };

    from_template("Finalize.sh.tt", $bash_file,
                  joint_name => $joint_name,
                  log_file => $log_file,
                  run_name => $run_name,
                  opt => $configuration);

    my $qsub = qsubTemplate(\%opt, "FINALIZE");
    if (@runningJobs) {
        system "$qsub -o /dev/null -e /dev/null -N Finalize_${job_id} -hold_jid " . join(",", @runningJobs) . " $bash_file";
    } else {
        system "$qsub -o /dev/null -e /dev/null -N Finalize_${job_id} $bash_file";
    }
}

1;
