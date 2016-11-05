package illumina_finalize;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_sge qw(qsubTemplate);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);
use illumina_metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(runFinalize);


sub runFinalize {
    my ($opt) = @_;

    say "\n### SCHEDULING PIPELINE FINALIZE ####";

    my $job_id = "$opt->{RUN_NAME}_" . getJobId();
    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "Finalize_${job_id}.sh");
    my $log_file = catfile($opt->{OUTPUT_DIR}, "logs", "PipelineCheck.log");

    my $joint_name = "";
    if ($opt->{SOMATIC_VARIANTS} eq "yes" || ($opt->{COPY_NUMBER} eq "yes" && $opt->{CNV_MODE} eq "sample_control")) {
        my $metadata = illumina_metadata::parse($opt);
        my $ref_sample = $metadata->{ref_sample};
        my $tumor_sample = $metadata->{tumor_sample};
        $joint_name = "${ref_sample}_${tumor_sample}";
    }

    my @runningJobs = map { @$_ }
        grep { defined } @{$opt->{RUNNING_JOBS}}{"baf", "prestats", keys %{$opt->{SAMPLES}}, "slicing", "poststats", "somvar", "cnv", "kinship"};

    from_template(
        "Finalize.sh.tt", $bash_file,
        joint_name => $joint_name,
        log_file => $log_file,
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "FINALIZE");
    if (@runningJobs) {
        system "$qsub -o /dev/null -e /dev/null -N Finalize_${job_id} -hold_jid " . join(",", @runningJobs) . " $bash_file";
    } else {
        system "$qsub -o /dev/null -e /dev/null -N Finalize_${job_id} $bash_file";
    }
    return;
}

1;
