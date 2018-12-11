package HMF::Pipeline::Purple;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs);
use HMF::Pipeline::Functions::Sge qw(qsubJava);
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Metadata qw(parse linkArtefact);

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

    my $circos_path = "${sub_dir}/plot/${tumor_sample}.circos.png";
    my $purple_cnv = "${sub_dir}/${tumor_sample}.purple.cnv";
    my $purple_gene_cnv = "${sub_dir}/${tumor_sample}.purple.gene.cnv";
    my $purple_germline_cnv = "${sub_dir}/${tumor_sample}.purple.germline.cnv";
    my $purple_purity = "${sub_dir}/${tumor_sample}.purple.purity";

    my $job_id = fromTemplate(
        "Purple",
        undef,
        1,
        qsubJava($opt, "PURPLE"),
        $dependent_jobs,
        $dirs,
        $opt,
        purple_purity_path => $purple_purity,
        # SABR: comment to avoid perltidy putting on one line
    );

    $opt->{RUNNING_JOBS}->{purple} = [$job_id] if $job_id;

    linkArtefact($circos_path, 'circos_plot', $opt);
    linkArtefact($purple_cnv, 'purple_cnv', $opt);
    linkArtefact($purple_gene_cnv, 'purple_gene_cnv', $opt);
    linkArtefact($purple_germline_cnv, 'purple_germline_cnv', $opt);
    linkArtefact($purple_purity, 'purple_purity', $opt);

    return;
}

sub dependencies {
    my ($opt) = @_;

    my @jobs;
    push @jobs, @{$opt->{RUNNING_JOBS}->{amber}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{cobalt}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{gridss}};
    push @jobs, @{$opt->{RUNNING_JOBS}->{strelka}};
    return \@jobs;
}

1;
