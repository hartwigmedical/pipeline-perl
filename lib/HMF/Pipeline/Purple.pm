package HMF::Pipeline::Purple;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamsAndJobs);
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Metadata qw(parse linkArtefact);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;
    my $metadata = parse($opt);
    my $tumor_sample = $metadata->{tumor_sample} or die "metadata missing tumor_sample";

    say "\n### SCHEDULING PURPLE ###";

    my $sub_dir = "purple";
    my $dirs = createDirs($opt->{OUTPUT_DIR}, purple => $sub_dir);
    my $dependent_jobs = dependencies($opt);
    my $job_id = fromTemplate(
        "Purple",
        undef,
        1,
        qsubJava($opt, "PURPLE"),
        $dependent_jobs,
        $dirs,
        $opt,
        # SABR: comment to avoid perltidy putting on one line
    );

    $opt->{RUNNING_JOBS}->{purple} = [$job_id] if $job_id;

    my $cir_png_path = "${sub_dir}/plot/${tumor_sample}.circos.png";
    my $cnv_csv_path = "${sub_dir}/${tumor_sample}.purple.cnv";
    my $pur_txt_path = "${sub_dir}/${tumor_sample}.purple.purity";

    linkArtefact($cir_png_path, 'somatic_circos_plot', $opt);
    linkArtefact($cnv_csv_path, 'somatic_copynumber_calls', $opt);
    linkArtefact($pur_txt_path, 'somatic_tumor_purity', $opt);

    return;
}

sub dependencies {
    my ($opt) = @_;

    my @jobs;
    push @jobs, @{$opt->{RUNNING_JOBS}->{amber}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{cobalt}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{sv}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{strelka}};
    return \@jobs;
}

1;
