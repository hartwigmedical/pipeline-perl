package HMF::Pipeline::GermlineCalling;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Sge qw(jobNative qsubJava);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING GERMLINE CALLING ###";

    my $job_id = "GermlineCalling_" . getId();
    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-f "$opt->{OUTPUT_DIR}/logs/GermlineCaller.done" || -f "$opt->{OUTPUT_DIR}/logs/VariantCaller.done") {
        say "WARNING: $opt->{OUTPUT_DIR}/logs/GermlineCaller.done exists, skipping $job_id";
        return;
    }

    my $gvcf_dir = catfile($opt->{OUTPUT_DIR}, "gvcf");
    if ((!-d $gvcf_dir && $opt->{CALLING_GVCF} eq "yes")) {
        mkdir($gvcf_dir) or die "Couldn't create directory $gvcf_dir: $!";
    }

    my @sample_bams;
    my @running_jobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $sample_bam = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $opt->{BAM_FILES}->{$sample});
        push @sample_bams, $sample_bam;
        push @running_jobs, @{$opt->{RUNNING_JOBS}->{$sample}} if @{$opt->{RUNNING_JOBS}->{$sample}};
    }

    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "${job_id}.sh");
    my $log_dir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $stdout = catfile($log_dir, "GermlineCaller_$opt->{RUN_NAME}.out");
    my $stderr = catfile($log_dir, "GermlineCaller_$opt->{RUN_NAME}.err");

    writeFromTemplate(
        "GermlineCalling.sh.tt", $bash_file,
        gvcf_dir => $gvcf_dir,
        sample_bams => \@sample_bams,
        job_native => jobNative($opt, "CALLING"),
        opt => $opt,
    );

    my $qsub = qsubJava($opt, "CALLING_MASTER");
    if (@running_jobs) {
        system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
    }

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }

    linkArtefacts($gvcf_dir, $opt);
    return;
}

# naming dependent on GermlineCaller.scala, could fix to be explicit
sub linkArtefacts {
    my ($gvcf_dir, $opt) = @_;

    if ($opt->{CALLING_GVCF} eq "yes") {
        foreach my $sample (keys %{$opt->{SAMPLES}}) {
            my $bam_file = $opt->{BAM_FILES}->{$sample};
            (my $gvcf_file = $bam_file) =~ s/\.bam$/.g.vcf.gz/;
            my $gvcf_path = catfile($gvcf_dir, $gvcf_file);
            my $sample_name = HMF::Pipeline::Metadata::metaSampleName($sample, $opt);
            HMF::Pipeline::Metadata::linkArtefact($gvcf_path, "${sample_name}_gvcf", $opt);
            HMF::Pipeline::Metadata::linkArtefact("${gvcf_path}.tbi", "${sample_name}_gvcf_index", $opt);
        }
    }
    my $germline_vcf_path = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.raw_variants.vcf");
    HMF::Pipeline::Metadata::linkArtefact($germline_vcf_path, "germline_vcf", $opt);
    HMF::Pipeline::Metadata::linkArtefact("${germline_vcf_path}.idx", "germline_vcf_index", $opt);
    return;
}

1;
