package HMF::Pipeline;

use FindBin::libs;
use discipline;

use Fcntl qw/O_WRONLY O_CREAT O_EXCL/;
use File::Basename;
use File::Path qw(make_path);
use File::Spec::Functions;

use HMF::Pipeline::Metadata;
use HMF::Pipeline::Prestats;
use HMF::Pipeline::Mapping;
use HMF::Pipeline::Poststats;
use HMF::Pipeline::Realignment;
use HMF::Pipeline::BaseRecalibration;
use HMF::Pipeline::GermlineCalling;
use HMF::Pipeline::GermlineFiltering;
use HMF::Pipeline::GermlineAnnotation;
use HMF::Pipeline::SomaticVariants;
use HMF::Pipeline::CopyNumber;
use HMF::Pipeline::Baf;
use HMF::Pipeline::Kinship;
use HMF::Pipeline::Finalize;

use parent qw(Exporter);
our @EXPORT_OK = qw(createOutputDirs lockRun run);
our $VERSION = 'v1.11';


sub run {
    my ($opt) = @_;

    if ($opt->{FASTQ}) {
        HMF::Pipeline::Prestats::run($opt) if $opt->{PRESTATS} eq "yes";
        HMF::Pipeline::Mapping::run($opt) if $opt->{MAPPING} eq "yes";
    } elsif ($opt->{BAM}) {
        HMF::Pipeline::Mapping::runBamPrep($opt);
    }

    if ($opt->{FASTQ} or $opt->{BAM}) {
        HMF::Pipeline::Poststats::run($opt) if $opt->{POSTSTATS} eq "yes";
        HMF::Pipeline::Realignment::run($opt) if $opt->{INDELREALIGNMENT} eq "yes";
        HMF::Pipeline::BaseRecalibration::run($opt) if $opt->{BASEQUALITYRECAL} eq "yes";
        linkBamArtefacts($opt);
        HMF::Pipeline::SomaticVariants::run($opt) if $opt->{SOMATIC_VARIANTS} eq "yes";
        HMF::Pipeline::CopyNumber::run($opt) if $opt->{COPY_NUMBER} eq "yes";
        HMF::Pipeline::Baf::run($opt) if $opt->{BAF} eq "yes";
        HMF::Pipeline::GermlineCalling::run($opt) if $opt->{VARIANT_CALLING} eq "yes";
        HMF::Pipeline::GermlineFiltering::run($opt) if $opt->{FILTER_VARIANTS} eq "yes";
        HMF::Pipeline::GermlineAnnotation::run($opt) if $opt->{ANNOTATE_VARIANTS} eq "yes";
        HMF::Pipeline::Kinship::run($opt) if $opt->{KINSHIP} eq "yes";
        HMF::Pipeline::Finalize::run($opt) if $opt->{FINALIZE} eq "yes";
        HMF::Pipeline::Metadata::writeLinks($opt);
    }
    return;
}

sub createOutputDirs {
    my ($output_dir, $samples) = @_;

    if (!-d $output_dir) {
        make_path($output_dir) or die "Couldn't create directory ${output_dir}: $!";
    }

    if (!-d "${output_dir}/QCStats") {
        mkdir("${output_dir}/QCStats") or die "Couldn't create directory ${output_dir}/QCStats: $!";
    }

    if (!-d "${output_dir}/jobs") {
        mkdir("${output_dir}/jobs") or die "Couldn't create directory ${output_dir}/jobs: $!";
    }

    if (!-d "${output_dir}/logs") {
        mkdir("${output_dir}/logs") or die "Couldn't create directory ${output_dir}/logs: $!";
    }

    if (!-d "${output_dir}/tmp") {
        mkdir("${output_dir}/tmp") or die "Couldn't create directory ${output_dir}/tmp: $!";
    }

    foreach my $sample (keys %{$samples}) {
        if (!-d "${output_dir}/$sample") {
            mkdir("${output_dir}/$sample") or die "Couldn't create directory ${output_dir}/$sample: $!";
        }

        if (!-d "${output_dir}/$sample/mapping") {
            mkdir("${output_dir}/$sample/mapping") or die "Couldn't create directory ${output_dir}/$sample/mapping: $!";
        }

        if (!-d "${output_dir}/$sample/QCStats") {
            mkdir("${output_dir}/$sample/QCStats") or die "Couldn't create directory ${output_dir}/$sample/QCStats: $!";
        }

        if (!-d "${output_dir}/$sample/jobs") {
            mkdir("${output_dir}/$sample/jobs") or die "Couldn't create directory ${output_dir}/$sample/jobs: $!";
        }

        if (!-d "${output_dir}/$sample/logs") {
            mkdir("${output_dir}/$sample/logs") or die "Couldn't create directory ${output_dir}/$sample/logs: $!";
        }

        if (!-d "${output_dir}/$sample/tmp") {
            mkdir("${output_dir}/$sample/tmp") or die "Couldn't create directory ${output_dir}/$sample/tmp: $!";
        }
    }
    return;
}

sub lockRun {
    my ($dir) = @_;
    my $lock_file = catfile($dir, "run.lock");
    ## no critic (Bangs::ProhibitBitwiseOperators)
    sysopen my $fh, $lock_file, O_WRONLY | O_CREAT | O_EXCL or die "Couldn't obtain lock file, are you *sure* there are no more jobs running? (error: $!)";
    ## use critic
    close $fh;
    return;
}

sub linkBamArtefacts {
    my ($opt) = @_;

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $bam_path = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $opt->{BAM_FILES}->{$sample});
        my $sample_name = HMF::Pipeline::Metadata::metaSampleName($sample, $opt);
        HMF::Pipeline::Metadata::linkArtefact($bam_path, "${sample_name}_bam", $opt);
        HMF::Pipeline::Metadata::linkArtefact("${bam_path}.bai", "${sample_name}_bai", $opt);
    }
    return;
}

1;

__END__

=pod

=encoding utf8

=head1 NAME

HMF::Pipeline - Analysis pipeline for whole genome sequencing of DNA

=head1 AUTHORS

=over

=item *

Robert Ernst <r.f.ernst-3@umcutrecht.nl>

=item *

Sam Brightman <s.brightman@hartwigmedicalfoundation.nl>

=item *

Korneel Duyvesteyn <korneel.duyvesteyn@gmail.com>

=item *

Thijs Houtenbos <thoutenbos@schubergphilis.com>

=item *

Arjen Wolfs <awolfs@schubergphilis.com>

=item *

Andrew Repton <arepton@schubergphilis.com>

=item *

Sander Boymans <s.boymans@hubrecht.eu>

=item *

Stef van Lieshout <stefvanlieshout@fastmail.fm>

=item *

Joep de Ligt <joepio@gmail.com>

=item *

Mark van Roosmalen <m.vanroosmalen-2@umcutrecht.nl>

=item *

Hindrik Kerstens <hindrik.kerstens@gmail.com>

=back

=cut
