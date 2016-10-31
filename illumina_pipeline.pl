#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use Getopt::Long;
use File::Path qw(make_path);
use File::Basename;
use File::Spec::Functions;
use Fcntl qw/O_WRONLY O_CREAT O_EXCL/;

use illumina_metadata;
use illumina_config qw(readConfig checkConfig setupLogging recordGitVersion copyConfigAndScripts);
use illumina_prestats;
use illumina_mapping;
use illumina_poststats;
use illumina_realign;
use illumina_baseRecal;
use illumina_germlineCalling;
use illumina_germlineFiltering;
use illumina_germlineAnnotation;
use illumina_somaticVariants;
use illumina_copyNumber;
use illumina_baf;
use illumina_kinship;
use illumina_finalize;


my $opt = {};
die usage() if @ARGV == 0;
readConfig($ARGV[0], $opt);
checkConfig($opt);
getSamples($opt);
createOutputDirs($opt->{OUTPUT_DIR}, $opt->{SAMPLES});
setupLogging($opt->{OUTPUT_DIR});
lockRun($opt->{OUTPUT_DIR});
recordGitVersion($opt);
copyConfigAndScripts($opt);
runPipeline($opt);


sub runPipeline {
    if ($opt->{FASTQ}) {
        illumina_prestats::runPreStats($opt) if $opt->{PRESTATS} eq "yes";
        illumina_mapping::runMapping($opt) if $opt->{MAPPING} eq "yes";
    } elsif ($opt->{BAM}) {
        illumina_mapping::runBamPrep($opt);
    }

    if ($opt->{FASTQ} or $opt->{BAM}) {
        illumina_poststats::runPostStats($opt) if $opt->{POSTSTATS} eq "yes";
        illumina_realign::runRealignment($opt) if $opt->{INDELREALIGNMENT} eq "yes";
        illumina_baseRecal::runBaseRecalibration($opt) if $opt->{BASEQUALITYRECAL} eq "yes";
        linkBamArtefacts($opt);
        illumina_somaticVariants::runSomaticVariantCallers($opt) if $opt->{SOMATIC_VARIANTS} eq "yes";
        illumina_copyNumber::runCopyNumberTools($opt) if $opt->{COPY_NUMBER} eq "yes";
        illumina_baf::runBAF($opt) if $opt->{BAF} eq "yes";
        illumina_germlineCalling::runVariantCalling($opt) if $opt->{VARIANT_CALLING} eq "yes";
        illumina_germlineFiltering::runFilterVariants($opt) if $opt->{FILTER_VARIANTS} eq "yes";
        illumina_germlineAnnotation::runAnnotateVariants($opt) if $opt->{ANNOTATE_VARIANTS} eq "yes";
        illumina_kinship::runKinship($opt) if $opt->{KINSHIP} eq "yes";
        illumina_finalize::runFinalize($opt) if $opt->{FINALIZE} eq "yes";
        illumina_metadata::writeLinks($opt);
    }
    return;
}

sub usage {
    warn <<"END";
    Usage: perl illumina_pipeline.pl configurationFile.conf
END
    exit;
}

sub getSamples {
    my ($opt) = @_;

    $opt->{SAMPLES} = {};

    if ($opt->{FASTQ}) {
        foreach my $input_file (keys %{$opt->{FASTQ}}) {
            my $fastqFile = fileparse($input_file);
            my ($sampleName) = split "_", $fastqFile;
            $opt->{SAMPLES}->{$sampleName} = $input_file;
            @{$opt->{RUNNING_JOBS}->{$sampleName}} = ();
        }
    }

    if ($opt->{BAM}) {
        foreach my $input_file (keys %{$opt->{BAM}}) {
            my $sampleName = illumina_mapping::verifyBam($input_file, $opt);
            not exists $opt->{SAMPLES}->{$sampleName} or die "sample $sampleName from $input_file already used by $opt->{SAMPLES}->{$sampleName}";
            $opt->{SAMPLES}->{$sampleName} = $input_file;
            @{$opt->{RUNNING_JOBS}->{$sampleName}} = ();
        }
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
    sysopen my $fh, $lock_file, O_WRONLY | O_CREAT | O_EXCL or die "Couldn't obtain lock file, are you *sure* there are no more jobs running? (error: $!)";
    close $fh;
    return;
}

sub linkBamArtefacts {
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $bam_path = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $opt->{BAM_FILES}->{$sample});
        my $sample_name = illumina_metadata::metaSampleName($sample, $opt);
        illumina_metadata::linkArtefact($bam_path, "${sample_name}_bam", $opt);
        illumina_metadata::linkArtefact("${bam_path}.bai", "${sample_name}_bai", $opt);
    }
    return;
}
