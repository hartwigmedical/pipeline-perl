package HMF::Pipeline;

use FindBin::libs;
use discipline;

use Fcntl qw/O_WRONLY O_CREAT O_EXCL/;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Metadata;

use HMF::Pipeline::PreStats;
use HMF::Pipeline::Mapping;
use HMF::Pipeline::Realignment;
use HMF::Pipeline::DamageEstimate;
use HMF::Pipeline::PostStats;
use HMF::Pipeline::Amber;
use HMF::Pipeline::Cobalt;
use HMF::Pipeline::GermlineCalling;
use HMF::Pipeline::Strelka;
use HMF::Pipeline::StructuralVariants;
use HMF::Pipeline::Purple;
use HMF::Pipeline::HealthCheck;
use HMF::Pipeline::PipelineCheck;

use parent qw(Exporter);
our @EXPORT_OK = qw(lockRun run);
our $VERSION = 'v4.1';

sub run {
    my ($opt) = @_;

    if ($opt->{FASTQ}) {
        HMF::Pipeline::PreStats::run($opt) if $opt->{PRESTATS} eq "yes";
        HMF::Pipeline::Mapping::run($opt) if $opt->{MAPPING} eq "yes";
    } elsif ($opt->{BAM}) {
        HMF::Pipeline::Mapping::runBamPrep($opt);
    }

    if (($opt->{FASTQ} and $opt->{MAPPING} eq "yes") or $opt->{BAM}) {
        HMF::Pipeline::PostStats::run($opt) if $opt->{POSTSTATS} eq "yes";
        HMF::Pipeline::Realignment::run($opt) if $opt->{INDEL_REALIGNMENT} eq "yes";
        # KODU: If we run from BAM and don't realign, we don't actually have a (new) BAM file to link up.
        HMF::Pipeline::Functions::Metadata::linkBamArtefacts($opt) if $opt->{INDEL_REALIGNMENT} eq "yes" or $opt->{FASTQ};

        HMF::Pipeline::Amber::run($opt) if $opt->{AMBER} eq "yes";
        # KODU: Amber BAF Segmentation is for rerunning only and expects (previously run) amber output.
        HMF::Pipeline::Amber::runAmberRerun($opt) if $opt->{AMBER_RERUN} eq "yes";
        HMF::Pipeline::Cobalt::run($opt) if $opt->{COBALT} eq "yes" or $opt->{COBALT_RERUN} eq "yes";
        HMF::Pipeline::DamageEstimate::run($opt) if $opt->{DAMAGE_ESTIMATE} eq "yes";

        HMF::Pipeline::GermlineCalling::run($opt) if $opt->{GERMLINE_CALLING} eq "yes";
        HMF::Pipeline::GermlineCalling::runGermlineRerun($opt) if $opt->{GERMLINE_RERUN} eq "yes";

        HMF::Pipeline::Strelka::run($opt) if $opt->{STRELKA} eq "yes";
        HMF::Pipeline::StructuralVariants::run($opt) if $opt->{STRUCTURAL_VARIANT_CALLING} eq "yes";
        HMF::Pipeline::Purple::run($opt) if $opt->{PURPLE} eq "yes";

        HMF::Pipeline::HealthCheck::run($opt) if $opt->{HEALTHCHECK} eq "yes";
        HMF::Pipeline::PipelineCheck::run($opt);
    }

    HMF::Pipeline::Functions::Metadata::writeLinks($opt);

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