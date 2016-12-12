package HMF::Pipeline::Finalize;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING PIPELINE FINALIZE ###";

    my $job_id = "$opt->{RUN_NAME}_" . getId();
    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my $bash_file = catfile($dirs->{job}, "Finalize_${job_id}.sh");
    my $log_file = catfile($dirs->{log}, "PipelineCheck.log");
    my $extras_tar = catfile($dirs->{out}, "$opt->{RUN_NAME}_extras.tar.gz");
    my $extras_zip = catfile($dirs->{out}, "$opt->{RUN_NAME}_extras.zip");

    my $joint_name = "";
    if ($opt->{SOMATIC_VARIANTS} eq "yes" || ($opt->{COPY_NUMBER} eq "yes" && $opt->{CNV_MODE} eq "sample_control")) {
        my $metadata = HMF::Pipeline::Metadata::parse($opt);
        my $ref_sample = $metadata->{ref_sample};
        my $tumor_sample = $metadata->{tumor_sample};
        $joint_name = "${ref_sample}_${tumor_sample}";
    }

    my @running_jobs = map { @$_ } grep { defined } @{$opt->{RUNNING_JOBS}}{
        # flatten all registered jobs into one list
        "baf",
        "prestats",
        keys %{$opt->{SAMPLES}},
        "slicing",
        "poststats",
        "somvar",
        "cnv",
        "sv",
        "kinship",
    };

    writeFromTemplate(
        "Finalize.sh.tt", $bash_file,
        joint_name => $joint_name,
        extras_tar => $extras_tar,
        extras_zip => $extras_zip,
        log_file => $log_file,
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "FINALIZE");
    my $stdout = catfile($dirs->{log}, "Finalize_$opt->{RUN_NAME}.out");
    my $stderr = catfile($dirs->{log}, "Finalize_$opt->{RUN_NAME}.err");
    if (@running_jobs) {
        system "$qsub -o $stdout -e $stderr -N Finalize_${job_id} -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N Finalize_${job_id} $bash_file";
    }

    HMF::Pipeline::Metadata::linkArtefact($extras_tar, "extras_tar", $opt) if $opt->{EXTRAS};
    HMF::Pipeline::Metadata::linkArtefact($extras_zip, "extras_zip", $opt) if $opt->{EXTRAS};

    return;
}

1;
