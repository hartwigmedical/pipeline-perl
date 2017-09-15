package HMF::Pipeline;

use FindBin::libs;
use discipline;

use Fcntl qw/O_WRONLY O_CREAT O_EXCL/;
use File::Spec::Functions;

use HMF::Pipeline::Amber;
use HMF::Pipeline::Baf;
use HMF::Pipeline::BaseRecalibration;
use HMF::Pipeline::CallableLoci;
use HMF::Pipeline::Cobalt;
use HMF::Pipeline::CopyNumber;
use HMF::Pipeline::Finalize;
use HMF::Pipeline::GermlineAnnotation;
use HMF::Pipeline::GermlineCalling;
use HMF::Pipeline::GermlineFiltering;
use HMF::Pipeline::Gender;
use HMF::Pipeline::Kinship;
use HMF::Pipeline::Mapping;
use HMF::Pipeline::Metadata;
use HMF::Pipeline::PostStats;
use HMF::Pipeline::PreCalling;
use HMF::Pipeline::PreStats;
use HMF::Pipeline::Purple;
use HMF::Pipeline::Realignment;
use HMF::Pipeline::SomaticVariants;
use HMF::Pipeline::StructuralVariants;
use HMF::Pipeline::DamageEstimate;
use HMF::Pipeline::PipelineCheck;
use HMF::Pipeline::HealthCheck;

use parent qw(Exporter);
our @EXPORT_OK = qw(lockRun run);
our $VERSION = 'v3.1';


sub run {
    my ($opt) = @_;

    if ($opt->{FASTQ}) {
        HMF::Pipeline::PreStats::run($opt) if $opt->{PRESTATS} eq "yes";
        HMF::Pipeline::Mapping::run($opt) if $opt->{MAPPING} eq "yes";
    } elsif ($opt->{BAM}) {
        HMF::Pipeline::Mapping::runBamPrep($opt);
    }

    if ($opt->{FASTQ} or $opt->{BAM}) {
        HMF::Pipeline::PostStats::run($opt) if $opt->{POSTSTATS} eq "yes";
        HMF::Pipeline::Realignment::run($opt) if $opt->{INDELREALIGNMENT} eq "yes";
        HMF::Pipeline::Metadata::linkBamArtefacts($opt);
        HMF::Pipeline::BaseRecalibration::run($opt) if $opt->{BASEQUALITYRECAL} eq "yes";
        HMF::Pipeline::Amber::run($opt) if $opt->{AMBER} eq "yes";
        HMF::Pipeline::Cobalt::run($opt) if $opt->{COBALT} eq "yes";

        HMF::Pipeline::PreCalling::run($opt);

        HMF::Pipeline::DamageEstimate::run($opt) if $opt->{DAMAGE_ESTIMATE} eq "yes";
        HMF::Pipeline::SomaticVariants::run($opt) if $opt->{SOMATIC_VARIANTS} eq "yes";
        HMF::Pipeline::CopyNumber::run($opt) if $opt->{COPY_NUMBER} eq "yes";
        HMF::Pipeline::StructuralVariants::run($opt) if $opt->{SV_CALLING} eq "yes";
        HMF::Pipeline::Baf::run($opt) if $opt->{BAF} eq "yes";
        HMF::Pipeline::CallableLoci::run($opt) if $opt->{CALLABLE_LOCI} eq "yes";
        HMF::Pipeline::GermlineCalling::run($opt) if $opt->{VARIANT_CALLING} eq "yes";
        HMF::Pipeline::GermlineFiltering::run($opt) if $opt->{FILTER_VARIANTS} eq "yes";
        HMF::Pipeline::GermlineAnnotation::run($opt) if $opt->{ANNOTATE_VARIANTS} eq "yes";
        HMF::Pipeline::Purple::run($opt) if $opt->{PURPLE} eq "yes";
        HMF::Pipeline::Gender::run($opt) if $opt->{GENDER} eq "yes";
        HMF::Pipeline::Kinship::run($opt) if $opt->{KINSHIP} eq "yes";
        HMF::Pipeline::PipelineCheck::run($opt);
        HMF::Pipeline::HealthCheck::run($opt) if $opt->{HEALTHCHECK} eq "yes";
        HMF::Pipeline::Finalize::run($opt) if $opt->{FINALIZE} eq "yes";

        HMF::Pipeline::Metadata::writeLinks($opt);
    }
    return;
}

sub lockRun {
    my ($dir) = @_;
    my $lock_file = catfile($dir, "run.lock");
    ## no critic (Bangs::ProhibitBitwiseOperators)
    sysopen my $fh, $lock_file, O_WRONLY | O_CREAT | O_EXCL
        or die "Couldn't obtain lock file (error: $!), are you *sure* there are no more jobs running? You MUST check qstat!
Jobs could have been scheduled before the error, or by external tools (GATK etc.)";
    ## use critic
    close $fh;
    return $lock_file;
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
