#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use Getopt::Long;
use Cwd qw(abs_path);
use File::Path qw(make_path);
use File::Copy::Recursive qw(rcopy);
use File::Basename;
use File::Spec::Functions;
use Fcntl qw/O_WRONLY O_CREAT O_EXCL/;
use IO::Pipe;
use POSIX qw(strftime);
use Time::HiRes qw(gettimeofday);

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
use illumina_metadata;


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

sub setupLogging {
    my ($output_dir) = @_;

    my ($seconds, $microseconds) = gettimeofday;
    my $datetime = strftime('%Y%m%d_%H%M%S_', localtime $seconds) . sprintf('%.6d', $microseconds);
    my $out_file = catfile($output_dir, "logs", "submitlog_${datetime}.out");
    my $err_file = catfile($output_dir, "logs", "submitlog_${datetime}.err");
    my $out_fh = IO::Pipe->new()->writer("tee $out_file") or die "Couldn't tee to $out_file: $!";
    my $err_fh = IO::Pipe->new()->writer("tee $err_file >&2") or die "Couldn't tee to $err_file: $!";
    open STDOUT, ">&", $out_fh or die "STDOUT redirection failed: $!";
    open STDERR, ">&", $err_fh or die "STDERR redirection failed: $!";
    STDOUT->autoflush(1);
    STDERR->autoflush(1);
    $out_fh->autoflush(1);
    return;
}

sub usage {
    warn <<"END";
    Usage: perl illumina_pipeline.pl configurationFile.conf
END
    exit;
}

sub readConfig {
    my ($configurationFile, $opt) = @_;

    open my $fh, "<", $configurationFile or die "Couldn't open $configurationFile: $!";
    while (<$fh>) {
        chomp;
        next if m/^#/ or not $_;
        my ($key, $val) = split "\t", $_, 2;
        warn "Key '$key' is missing a value - is it badly formatted?" if not defined $val;

        if ($key eq 'INIFILE') {
            $val = catfile(dirname(abs_path($0)), $val) unless file_name_is_absolute($val);
            push @{$opt->{$key}}, $val;
            readConfig($val, $opt);
        } elsif ($key eq 'FASTQ' or $key eq 'BAM') {
            $opt->{$key}->{$val} = 1;
        } else {
            $opt->{$key} = $val;
        }
    }
    close $fh;
    return;
}

sub checkConfig {
    my ($opt) = @_;
    my %opt = %{$opt};

    my $checkFailed = 0;

    ### Input and Output
    if (!$opt{INIFILE}) { say "ERROR: No INIFILE option found in config files."; $checkFailed = 1; }
    if (!$opt{OUTPUT_DIR}) { say "ERROR: No OUTPUT_DIR found in config files."; $checkFailed = 1; }
    if (!($opt{FASTQ} || $opt{BAM})) { say "ERROR: No FASTQ or BAM files found in config files."; $checkFailed = 1; }
    if (!$opt{MAIL}) { say "ERROR: No MAIL address specified in config files."; $checkFailed = 1; }

    ### Cluster settings
    if (!$opt{CLUSTER_PATH}) { say "ERROR: No CLUSTER_PATH option found in config files."; $checkFailed = 1; }
    if (!$opt{CLUSTER_TMP}) { say "ERROR: No CLUSTER_TMP option found in config files."; $checkFailed = 1; }
    if (!$opt{CLUSTER_RESERVATION}) { say "ERROR: No CLUSTER_RESERVATION option found in config files."; $checkFailed = 1; }
    if (!$opt{CLUSTER_PROJECT}) { say "ERROR: No CLUSTER_PROJECT option found in config files."; $checkFailed = 1; }

    ### Module Yes or No
    if (!$opt{PRESTATS}) { say "ERROR: No PRESTATS option found in config files."; $checkFailed = 1; }
    if (!$opt{MAPPING}) { say "ERROR: No MAPPING option found in config files."; $checkFailed = 1; }
    if (!$opt{POSTSTATS}) { say "ERROR: No POSTSTATS option found in config files."; $checkFailed = 1; }
    if (!$opt{INDELREALIGNMENT}) { say "ERROR: No INDELREALIGNMENT option found in config files."; $checkFailed = 1; }
    if (!$opt{BASEQUALITYRECAL}) { say "ERROR: No BASEQUALITYRECAL option found in config files."; $checkFailed = 1; }
    if (!$opt{VARIANT_CALLING}) { say "ERROR: No VARIANT_CALLING option found in config files."; $checkFailed = 1; }
    if (!$opt{FILTER_VARIANTS}) { say "ERROR: No FILTER_VARIANTS option found in config files."; $checkFailed = 1; }
    if (!$opt{SOMATIC_VARIANTS}) { say "ERROR: No SOMATIC_VARIANTS option found in config files."; $checkFailed = 1; }
    if (!$opt{COPY_NUMBER}) { say "ERROR: No COPY_NUMBER option found in config files."; $checkFailed = 1; }
    if (!$opt{BAF}) { say "ERROR: No BAF option found in config files."; $checkFailed = 1; }
    if (!$opt{ANNOTATE_VARIANTS}) { say "ERROR: No ANNOTATE_VARIANTS option found in config files."; $checkFailed = 1; }
    if (!$opt{KINSHIP}) { say "ERROR: No KINSHIP option found in config files."; $checkFailed = 1; }
    if (!$opt{FINALIZE}) { say "ERROR: No FINALIZE option found in config files."; $checkFailed = 1; }

    ### Module Settings / tools
    if (!$opt{GENOME}) { say "ERROR: No GENOME option found in config files."; $checkFailed = 1; }
    elsif (!-f $opt{GENOME}) { say "ERROR: $opt{GENOME} does not exist"}
    if (!$opt{SAMBAMBA_PATH}) { say "ERROR: No SAMBAMBA_PATH option found in config files."; $checkFailed = 1; }
    if (!$opt{QUEUE_PATH}) { say "ERROR: No QUEUE_PATH option found in config files."; $checkFailed = 1; }

    ## PRESTATS
    if ($opt{PRESTATS} && $opt{PRESTATS} eq "yes") {
        if (!$opt{FASTQC_PATH}) { say "ERROR: No FASTQC_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{PRESTATS_THREADS}) { say "ERROR: No PRESTATS_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{PRESTATS_MEM}) { say "ERROR: No PRESTATS_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{PRESTATS_QUEUE}) { say "ERROR: No PRESTATS_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{PRESTATS_TIME}) { say "ERROR: No PRESTATS_TIME option found in config files."; $checkFailed = 1; }
    }

    ## MAPPING
    if ($opt{MAPPING} && $opt{MAPPING} eq "yes") {
        if (!$opt{BWA_PATH}) { say "ERROR: No BWA_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{MAPPING_THREADS}) { say "ERROR: No MAPPING_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{MAPPING_MEM}) { say "ERROR: No MAPPING_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{MAPPING_QUEUE}) { say "ERROR: No MAPPING_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{MAPPING_TIME}) { say "ERROR: No MAPPING_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{MAPPING_SETTINGS}) { say "ERROR: No MAPPING_SETTINGS option found in config files."; $checkFailed = 1; }

        if (!$opt{MARKDUP_QUEUE}) { say "ERROR: No MARKDUP_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{MARKDUP_TIME}) { say "ERROR: No MARKDUP_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{MARKDUP_THREADS}) { say "ERROR: No MARKDUP_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{MARKDUP_MEM}) { say "ERROR: No MARKDUP_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{MARKDUP_OVERFLOW_LIST_SIZE}) { say "ERROR: No MARKDUP_OVERFLOW_LIST_SIZE option found in config files."; $checkFailed = 1; }
    }

    if (!$opt{FLAGSTAT_QUEUE}) { say "ERROR: No FLAGSTAT_QUEUE option found in config files."; $checkFailed = 1; }
    if (!$opt{FLAGSTAT_THREADS}) { say "ERROR: No FLAGSTAT_THREADS option found in config files."; $checkFailed = 1; }
    if (!$opt{FLAGSTAT_MEM}) { say "ERROR: No FLAGSTAT_MEM option found in config files."; $checkFailed = 1; }
    if (!$opt{FLAGSTAT_TIME}) { say "ERROR: No FLAGSTAT_TIME option found in config files."; $checkFailed = 1; }

    ## POSTSTATS
    if ($opt{POSTSTATS} && $opt{POSTSTATS} eq "yes") {
        if (!$opt{BAMMETRICS_PATH}) { say "ERROR: No BAMMETRICS_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{PICARD_PATH}) { say "ERROR: No PICARD_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{POSTSTATS_THREADS}) { say "ERROR: No POSTSTATS_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{POSTSTATS_MEM}) { say "ERROR: No POSTSTATS_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{POSTSTATS_QUEUE}) { say "ERROR: No POSTSTATS_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{POSTSTATS_TIME}) { say "ERROR: No POSTSTATS_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{EXONCALLCOV}) { say "ERROR: No EXONCALLCOV option found in config files."; $checkFailed = 1; }
        if ($opt{EXONCALLCOV} eq "yes") {
            if (!$opt{EXONCALLCOV_QUEUE}) { say "ERROR: No EXONCALLCOV_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{EXONCALLCOV_TIME}) { say "ERROR: No EXONCALLCOV_TIME option found in config files."; $checkFailed = 1; }
            if (!$opt{EXONCALLCOV_MEM}) { say "ERROR: No EXONCALLCOV_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{EXONCALLCOV_PATH}) { say "ERROR: No EXONCALLCOV_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{EXONCALLCOV_BED}) { say "ERROR: No EXONCALLCOV_BED option found in config files."; $checkFailed = 1; }
            if ($opt{EXONCALLCOV_BED} && !-f $opt{EXONCALLCOV_BED}) { say "ERROR: $opt{EXONCALLCOV_BED} does not exist."; $checkFailed = 1; }
            if (!$opt{EXONCALLCOV_PREF}) { say "ERROR: No EXONCALLCOV_PREF option found in config files."; $checkFailed = 1; }
            if ($opt{EXONCALLCOV_PREF} && !-f $opt{EXONCALLCOV_PREF}) { say "ERROR: $opt{EXONCALLCOV_PREF} does not exist."; $checkFailed = 1; }
            if (!$opt{EXONCALLCOV_PANEL}) { say "ERROR: No EXONCALLCOV_PANEL option found in config files."; $checkFailed = 1; }
            if ($opt{EXONCALLCOV_PANEL} && !-f $opt{EXONCALLCOV_PANEL}) { say "ERROR: $opt{EXONCALLCOV_PANEL} does not exist."; $checkFailed = 1; }
            if (!$opt{EXONCALLCOV_ENS}) { say "ERROR: No EXONCALLCOV_ENS option found in config files."; $checkFailed = 1; }
            if ($opt{EXONCALLCOV_ENS} && !-f $opt{EXONCALLCOV_ENS}) { say "ERROR: $opt{EXONCALLCOV_ENS} does not exist."; $checkFailed = 1; }
        }
    }

    ## INDELREALIGNMENT
    if ($opt{INDELREALIGNMENT} && $opt{INDELREALIGNMENT} eq "yes") {
        if (!$opt{BAMUTIL_PATH}) { say "ERROR: No BAMUTIL_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_MASTER_QUEUE}) { say "ERROR: No REALIGNMENT_MASTER_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_MASTER_THREADS}) { say "ERROR: No REALIGNMENT_MASTER_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_MASTER_TIME}) { say "ERROR: No REALIGNMENT_MASTER_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_MASTER_MEM}) { say "ERROR: No REALIGNMENT_MASTER_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_QUEUE}) { say "ERROR: No REALIGNMENT_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_THREADS}) { say "ERROR: No REALIGNMENT_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_MEM}) { say "ERROR: No REALIGNMENT_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_TIME}) { say "ERROR: No REALIGNMENT_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_MERGETHREADS}) { say "ERROR: No REALIGNMENT_MERGETHREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_SCALA}) { say "ERROR: No REALIGNMENT_SCALA option found in config files."; $checkFailed = 1; }
        if (!$opt{REALIGNMENT_SCATTER}) { say "ERROR: No REALIGNMENT_SCATTER option found in config files."; $checkFailed = 1; }
        if ($opt{REALIGNMENT_KNOWN} && grep { !-f } split "\t", $opt{REALIGNMENT_KNOWN}) { say "ERROR: Some of $opt{REALIGNMENT_KNOWN} do not exist."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_QUEUE}) { say "ERROR: No FLAGSTAT_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_THREADS}) { say "ERROR: No FLAGSTAT_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_MEM}) { say "ERROR: No FLAGSTAT_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_TIME}) { say "ERROR: No FLAGSTAT_TIME option found in config files."; $checkFailed = 1; }
    }

    ## BASEQUALITYRECAL
    if ($opt{BASEQUALITYRECAL} && $opt{BASEQUALITYRECAL} eq "yes") {
        if (!$opt{BASERECALIBRATION_MASTER_QUEUE}) { say "ERROR: No BASERECALIBRATION_MASTER_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_MASTER_TIME}) { say "ERROR: No BASERECALIBRATION_MASTER_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_MASTER_THREADS}) { say "ERROR: No BASERECALIBRATION_MASTER_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_MASTER_MEM}) { say "ERROR: No BASERECALIBRATION_MASTER_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_QUEUE}) { say "ERROR: No BASERECALIBRATION_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_THREADS}) { say "ERROR: No BASERECALIBRATION_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_MEM}) { say "ERROR: No BASERECALIBRATION_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_TIME}) { say "ERROR: No BASERECALIBRATION_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_SCALA}) { say "ERROR: No BASERECALIBRATION_SCALA option found in config files."; $checkFailed = 1; }
        if (!$opt{BASERECALIBRATION_SCATTER}) { say "ERROR: No BASERECALIBRATION_SCATTER option found in config files."; $checkFailed = 1; }
        if ($opt{BASERECALIBRATION_KNOWN} && grep { !-f } split "\t", $opt{BASERECALIBRATION_KNOWN}) { say "ERROR: Some of $opt{BASERECALIBRATION_KNOWN} do not exist."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_QUEUE}) { say "ERROR: No FLAGSTAT_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_THREADS}) { say "ERROR: No FLAGSTAT_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_MEM}) { say "ERROR: No FLAGSTAT_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{FLAGSTAT_TIME}) { say "ERROR: No FLAGSTAT_TIME option found in config files."; $checkFailed = 1; }
    }

    ## VARIANT_CALLING
    if ($opt{VARIANT_CALLING} && $opt{VARIANT_CALLING} eq "yes") {
        if (!$opt{CALLING_MASTER_QUEUE}) { say "ERROR: No CALLING_MASTER_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_MASTER_TIME}) { say "ERROR: No CALLING_MASTER_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_MASTER_THREADS}) { say "ERROR: No CALLING_MASTER_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_MASTER_MEM}) { say "ERROR: No CALLING_MASTER_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_QUEUE}) { say "ERROR: No CALLING_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_THREADS}) { say "ERROR: No CALLING_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_MEM}) { say "ERROR: No CALLING_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_TIME}) { say "ERROR: No CALLING_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_SCATTER}) { say "ERROR: No CALLING_SCATTER option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_GVCF}) { say "ERROR: No CALLING_GVCF option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_SCALA}) { say "ERROR: No CALLING_SCALA option found in config files."; $checkFailed = 1; }
        if ($opt{CALLING_UGMODE}) {
            if ($opt{CALLING_UGMODE} ne "SNP" and $opt{CALLING_UGMODE} ne "INDEL" and $opt{CALLING_UGMODE} ne "BOTH") { say "ERROR: UGMODE: $opt{CALLING_UGMODE} does not exist use SNP, INDEL or BOTH"; $checkFailed = 1; }
        }
        if (!$opt{CALLING_STANDCALLCONF}) { say "ERROR: No CALLING_STANDCALLCONF option found in config files."; $checkFailed = 1; }
        if (!$opt{CALLING_STANDEMITCONF}) { say "ERROR: No CALLING_STANDEMITCONF option found in config files."; $checkFailed = 1; }
        if ($opt{CALLING_TARGETS} && !-f $opt{CALLING_TARGETS}) { say "ERROR: $opt{CALLING_TARGETS} does not exist"; $checkFailed = 1; }
        if ($opt{CALLING_DBSNP} && !-f $opt{CALLING_DBSNP}) { say "ERROR: $opt{CALLING_DBSNP} does not exist"; $checkFailed = 1; }
    }

    ## FILTER_VARIANTS
    if ($opt{FILTER_VARIANTS} && $opt{FILTER_VARIANTS} eq "yes") {
        if (!$opt{FILTER_MASTER_QUEUE}) { say "ERROR: No FILTER_MASTER_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_MASTER_TIME}) { say "ERROR: No FILTER_MASTER_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_MASTER_THREADS}) { say "ERROR: No FILTER_MASTER_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_MASTER_MEM}) { say "ERROR: No FILTER_MASTER_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_QUEUE}) { say "ERROR: No FILTER_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_THREADS}) { say "ERROR: No FILTER_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_MEM}) { say "ERROR: No FILTER_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_TIME}) { say "ERROR: No FILTER_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_SCATTER}) { say "ERROR: No FILTER_SCATTER option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_SCALA}) { say "ERROR: No FILTER_SCALA option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_SNPTYPES}) { say "ERROR: No FILTER_SNPTYPES option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_SNPNAME}) { say "ERROR: No FILTER_SNPNAME option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_SNPEXPR}) { say "ERROR: No FILTER_SNPEXPR  option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_INDELTYPES}) { say "ERROR: No FILTER_INDELTYPES option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_INDELNAME}) { say "ERROR: No FILTER_INDELNAME option found in config files."; $checkFailed = 1; }
        if (!$opt{FILTER_INDELEXPR}) { say "ERROR: No FILTER_INDELEXPR option found in config files."; $checkFailed = 1; }
    }

    ## SOMATIC_VARIANTS
    if ($opt{SOMATIC_VARIANTS} && $opt{SOMATIC_VARIANTS} eq "yes") {
        if (!$opt{VCFTOOLS_PATH}) { say "ERROR: No VCFTOOLS_PATH found in .ini file"; $checkFailed = 1; }
        if (!$opt{SAMTOOLS_PATH}) { say "ERROR: No SAMTOOLS_PATH option found in config files."; $checkFailed = 1; }
        if ($opt{SOMVAR_TARGETS} && !-f $opt{SOMVAR_TARGETS}) { say "ERROR: $opt{SOMVAR_TARGETS} does not exist"; $checkFailed = 1; }
        if (!$opt{SOMVAR_STRELKA}) { say "ERROR: No SOMVAR_STRELKA option found in config files."; $checkFailed = 1; }
        if ($opt{SOMVAR_STRELKA} && $opt{SOMVAR_STRELKA} eq "yes") {
            if (!$opt{STRELKA_PATH}) { say "ERROR: No STRELKA_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{STRELKA_INI}) { say "ERROR: No STRELKA_INI option found in config files."; $checkFailed = 1; }
            if (!$opt{STRELKA_QUEUE}) { say "ERROR: No STRELKA_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{STRELKA_THREADS}) { say "ERROR: No STRELKA_THREADS option found in config files."; $checkFailed = 1; }
            if (!$opt{STRELKA_MEM}) { say "ERROR: No STRELKA_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{STRELKA_TIME}) { say "ERROR: No STRELKA_TIME option found in config files."; $checkFailed = 1; }
        }
        if (!$opt{SOMVAR_VARSCAN}) { say "ERROR: No SOMVAR_VARSCAN option found in config files."; $checkFailed = 1; }
        if ($opt{SOMVAR_VARSCAN} && $opt{SOMVAR_VARSCAN} eq "yes") {
            if (!$opt{VARSCAN_PATH}) { say "ERROR: No VARSCAN_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{PBGZIP_PATH}) { say "ERROR: No PBGZIP_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{TABIX_PATH}) { say "ERROR: No TABIX_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{VARSCAN_QUEUE}) { say "ERROR: No VARSCAN_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{VARSCAN_THREADS}) { say "ERROR: No VARSCAN_THREADS option found in config files."; $checkFailed = 1; }
            if (!$opt{VARSCAN_TIME}) { say "ERROR: No VARSCAN_TIME option found in config files."; $checkFailed = 1; }
            if (!$opt{VARSCAN_MEM}) { say "ERROR: No VARSCAN_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{VARSCAN_SETTINGS}) { say "ERROR: No VARSCAN_SETTINGS option found in config files."; $checkFailed = 1; }
            if (!$opt{VARSCAN_POSTSETTINGS}) { say "ERROR: No VARSCAN_POSTSETTINGS option found in config files."; $checkFailed = 1; }
            if (!$opt{PILEUP_QUEUE}) { say "ERROR: No PILEUP_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{PILEUP_DIVISOR}) { say "ERROR: No PILEUP_DIVISOR option found in config files."; $checkFailed = 1; }
            if (!$opt{PILEUP_THREADS}) { say "ERROR: No PILEUP_THREADS option found in config files."; $checkFailed = 1; }
            elsif ($opt{PILEUP_THREADS} < $opt{PILEUP_DIVISOR}) { say "ERROR: PILEUP_THREADS ($opt{PILEUP_THREADS}) must be at least PILEUP_DIVISOR ($opt{PILEUP_DIVISOR})."; $checkFailed = 1; }
            if (!$opt{PILEUP_MEM}) { say "ERROR: No PILEUP_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{PILEUP_TIME}) { say "ERROR: No PILEUP_TIME option found in config files."; $checkFailed = 1; }
            if (!$opt{FINALIZE_KEEP_PILEUP}) { say "ERROR: No FINALIZE_KEEP_PILEUP found in .ini file"; $checkFailed = 1; }
        }
        if (!$opt{SOMVAR_FREEBAYES}) { say "ERROR: No SOMVAR_FREEBAYES option found in config files."; $checkFailed = 1; }
        if ($opt{SOMVAR_FREEBAYES} && $opt{SOMVAR_FREEBAYES} eq "yes") {
            if (!$opt{FREEBAYES_PATH}) { say "ERROR: No FREEBAYES_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{VCFLIB_PATH}) { say "ERROR: No VCFLIB_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEBAYES_QUEUE}) { say "ERROR: No FREEBAYES_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEBAYES_THREADS}) { say "ERROR: No FREEBAYES_THREADS option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEBAYES_MEM}) { say "ERROR: No FREEBAYES_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEBAYES_TIME}) { say "ERROR: No FREEBAYES_TIME option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEBAYES_SETTINGS}) { say "ERROR: No FREEBAYES_SETTINGS option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEBAYES_SOMATICFILTER}) { say "ERROR: No FREEBAYES_SOMATICFILTER option found in config files."; $checkFailed = 1; }
        }
        if (!$opt{SOMVAR_MUTECT}) { say "ERROR: No SOMVAR_MUTECT option found in config files."; $checkFailed = 1; }
        if ($opt{SOMVAR_MUTECT} && $opt{SOMVAR_MUTECT} eq "yes") {
            if (!$opt{MUTECT_PATH}) { say "ERROR: No MUTECT_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{MUTECT_QUEUE}) { say "ERROR: No MUTECT_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{MUTECT_THREADS}) { say "ERROR: No MUTECT_THREADS option found in config files."; $checkFailed = 1; }
            if (!$opt{MUTECT_MEM}) { say "ERROR: No MUTECT_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{MUTECT_TIME}) { say "ERROR: No MUTECT_TIME option found in config files."; $checkFailed = 1; }
            if (!$opt{MUTECT_COSMIC}) { say "ERROR: No MUTECT_COSMIC option found in config files."; $checkFailed = 1; }
        }
        if (!$opt{SOMVARMERGE_QUEUE}) { say "ERROR: No SOMVARMERGE_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{SOMVARMERGE_THREADS}) { say "ERROR: No SOMVARMERGE_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{SOMVARMERGE_MEM}) { say "ERROR: No SOMVARMERGE_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{SOMVARMERGE_TIME}) { say "ERROR: No SOMVARMERGE_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{SOMVAR_ANNOTATE}) { say "ERROR: No SOMVAR_ANNOTATE option found in config files."; $checkFailed = 1; }
        if ($opt{SOMVAR_ANNOTATE} && $opt{SOMVAR_ANNOTATE} eq "yes") {
            if (!$opt{ANNOTATE_DB}) { say "ERROR: No ANNOTATE_DB option found in config files."; $checkFailed = 1; }
            if (!$opt{ANNOTATE_FLAGS}) { say "ERROR: No ANNOTATE_FLAGS option found in config files."; $checkFailed = 1; }
            if (!$opt{ANNOTATE_IDNAME}) { say "ERROR: No ANNOTATE_IDNAME option found in config files."; $checkFailed = 1; }
            if (!$opt{ANNOTATE_IDDB}) { say "ERROR: No ANNOTATE_IDDB option found in config files."; $checkFailed = 1; }
            if (!$opt{CALLING_DBSNP}) { say "ERROR: No CALLING_DBSNP option found in config files."; $checkFailed = 1; }
        }
    }

    ## COPY_NUMBER
    if ($opt{COPY_NUMBER} && $opt{COPY_NUMBER} eq "yes") {
        if (!$opt{CNVCHECK_QUEUE}) { say "ERROR: No CNVCHECK_QUEUE in config files."; $checkFailed = 1; }
        if (!$opt{CNVCHECK_THREADS}) { say "ERROR: No CNVCHECK_THREADS  in config files."; $checkFailed = 1; }
        if (!$opt{CNVCHECK_MEM}) { say "ERROR: No CNVCHECK_MEM in config files."; $checkFailed = 1; }
        if (!$opt{CNVCHECK_TIME}) { say "ERROR: No CNVCHECK_TIME in config files."; $checkFailed = 1; }
        if (!$opt{CNV_MODE}) { say "ERROR: No CNV_MODE in config files."; $checkFailed = 1; }
        if (!$opt{CNV_FREEC}) { say "ERROR: No CNV_FREEC  in config files."; $checkFailed = 1; }
        if ($opt{CNV_FREEC} eq "yes") {
            if (!$opt{FREEC_PATH}) { say "ERROR: No FREEC_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_QUEUE}) { say "ERROR: No FREEC_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_THREADS}) { say "ERROR: No FREEC_THREADS option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_MEM}) { say "ERROR: No FREEC_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_TIME}) { say "ERROR: No FREEC_TIME option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_CHRLENFILE}) { say "ERROR: No FREEC_CHRLENFILE option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_CHRFILES}) { say "ERROR: No FREEC_CHRFILES option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_PLOIDY}) { say "ERROR: No FREEC_PLOIDY option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_WINDOW}) { say "ERROR: No FREEC_WINDOW option found in config files."; $checkFailed = 1; }
            if (!$opt{FREEC_TELOCENTROMERIC}) { say "ERROR: No FREEC_TELOCENTROMERIC option found in config files."; $checkFailed = 1; }
        }
        if (!$opt{CNV_QDNASEQ}) { say "ERROR: No CNV_QDNASEQ in config files."; $checkFailed = 1; }
        if ($opt{CNV_QDNASEQ} eq "yes") {
            if (!$opt{QDNASEQ_PATH}) { say "ERROR: No QDNASEQ_PATH option found in config files."; $checkFailed = 1; }
            if (!$opt{QDNASEQ_QUEUE}) { say "ERROR: No QDNASEQ_QUEUE option found in config files."; $checkFailed = 1; }
            if (!$opt{QDNASEQ_THREADS}) { say "ERROR: No QDNASEQ_THREADS option found in config files."; $checkFailed = 1; }
            if (!$opt{QDNASEQ_MEM}) { say "ERROR: No QDNASEQ_MEM option found in config files."; $checkFailed = 1; }
            if (!$opt{QDNASEQ_TIME}) { say "ERROR: No QDNASEQ_TIME option found in config files."; $checkFailed = 1; }
        }
    }

    ##BAF Analysis
    if ($opt{BAF} && $opt{BAF} eq "yes") {
        if (!$opt{BAF_QUEUE}) { say "ERROR: No BAF_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{BAF_THREADS}) { say "ERROR: No BAF_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{BAF_MEM}) { say "ERROR: No BAF_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{BAF_TIME}) { say "ERROR: No BAF_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{BIOVCF_PATH}) { say "ERROR: No BIOVCF_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{BAF_SNPS}) { say "ERROR: No BAF_SNPS option found in config files."; $checkFailed = 1; }
    }

    ## ANNOTATE_VARIANTS
    if ($opt{ANNOTATE_VARIANTS} && $opt{ANNOTATE_VARIANTS} eq "yes") {
        if (!$opt{SNPEFF_PATH}) { say "ERROR: No SNPEFF_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{IGVTOOLS_PATH}) { say "ERROR: No IGVTOOLS_PATH option found in config files."; $checkFailed = 1; }
        if (!$opt{ANNOTATE_QUEUE}) { say "ERROR: No ANNOTATE_QUEUE option found in config files."; $checkFailed = 1; }
        if (!$opt{ANNOTATE_THREADS}) { say "ERROR: No ANNOTATE_THREADS option found in config files."; $checkFailed = 1; }
        if (!$opt{ANNOTATE_MEM}) { say "ERROR: No ANNOTATE_MEM option found in config files."; $checkFailed = 1; }
        if (!$opt{ANNOTATE_TIME}) { say "ERROR: No ANNOTATE_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{ANNOTATE_SNPEFF}) { say "ERROR: No ANNOTATE_SNPEFF option found in config files."; $checkFailed = 1; }
        if ($opt{ANNOTATE_SNPEFF} eq "yes") {
            if (!$opt{ANNOTATE_DB}) { say "ERROR: No ANNOTATE_DB option found in config files."; $checkFailed = 1; }
            if (!$opt{ANNOTATE_FLAGS}) { say "ERROR: No ANNOTATE_FLAGS option found in config files."; $checkFailed = 1; }
        }
        if (!$opt{ANNOTATE_SNPSIFT}) { say "ERROR: No ANNOTATE_SNPSIFT option found in config files."; $checkFailed = 1; }
        if ($opt{ANNOTATE_SNPSIFT} eq "yes") {
            if (!$opt{ANNOTATE_DBNSFP}) { say "ERROR: No ANNOTATE_DBNSFP option found in config files."; $checkFailed = 1; }
            elsif ($opt{ANNOTATE_DBNSFP} && !-f $opt{ANNOTATE_DBNSFP}) { say "ERROR: $opt{ANNOTATE_DBNSFP} does not exist"; $checkFailed = 1; }
            if (!$opt{ANNOTATE_FIELDS}) { say "ERROR: No ANNOTATE_FIELDS option found in config files."; $checkFailed = 1; }
        }
        if (!$opt{ANNOTATE_FREQUENCIES}) { say "ERROR: No ANNOTATE_FREQUENCIES option found in config files."; $checkFailed = 1; }
        if ($opt{ANNOTATE_FREQUENCIES} eq "yes") {
            if (!$opt{ANNOTATE_FREQNAME}) { say "ERROR: No ANNOTATE_FREQNAME option found in config files."; $checkFailed = 1; }
            if (!$opt{ANNOTATE_FREQDB}) { say "ERROR: No ANNOTATE_FREQDB option found in config files."; $checkFailed = 1; }
            elsif ($opt{ANNOTATE_FREQDB} && !-f $opt{ANNOTATE_FREQDB}) { say "ERROR: $opt{ANNOTATE_FREQDB} does not exist"; $checkFailed = 1; }
            if (!$opt{ANNOTATE_FREQINFO}) { say "ERROR: No ANNOTATE_FREQINFO option found in config files."; $checkFailed = 1; }
        }
        if (!$opt{ANNOTATE_IDFIELD}) { say "ERROR: No ANNOTATE_IDFIELD option found in config files."; $checkFailed = 1; }
            if ($opt{ANNOTATE_IDFIELD} eq "yes") {
            if (!$opt{ANNOTATE_IDNAME}) { say "ERROR: No ANNOTATE_IDNAME option found in config files."; $checkFailed = 1; }
            if (!$opt{ANNOTATE_IDDB}) { say "ERROR: No ANNOTATE_IDDB option found in config files."; $checkFailed = 1; }
            elsif ($opt{ANNOTATE_IDDB} && !-f $opt{ANNOTATE_IDDB}) { say "ERROR: $opt{ANNOTATE_IDDB} does not exist"; $checkFailed = 1; }
        }
    }

    ## KINSHIP
    if ($opt{KINSHIP} && $opt{KINSHIP} eq "yes") {
        if (!$opt{KINSHIP_QUEUE}) { say "ERROR: No KINSHIP_QUEUE found in .ini file"; $checkFailed = 1; }
        if (!$opt{KINSHIP_THREADS}) { say "ERROR: No KINSHIP_THREADS found in .ini file"; $checkFailed = 1; }
        if (!$opt{KINSHIP_MEM}) { say "ERROR: No KINSHIP_MEM found in .ini file"; $checkFailed = 1; }
        if (!$opt{KINSHIP_TIME}) { say "ERROR: No KINSHIP_TIME option found in config files."; $checkFailed = 1; }
        if (!$opt{PLINK_PATH}) { say "ERROR: No PLINK_PATH found in .ini file"; $checkFailed = 1; }
        if (!$opt{KING_PATH}) { say "ERROR: No KING_PATH found in .ini file"; $checkFailed = 1; }
        if (!$opt{VCFTOOLS_PATH}) { say "ERROR: No VCFTOOLS_PATH found in .ini file"; $checkFailed = 1; }
    }

    ## FINALIZE
    if ($opt{FINALIZE} && $opt{FINALIZE} eq "yes") {
        if (!$opt{FINALIZE_QUEUE}) { say "ERROR: No FINALIZE_QUEUE found in .ini file"; $checkFailed = 1; }
        if (!$opt{FINALIZE_THREADS}) { say "ERROR: No FINALIZE_THREADS found in .ini file"; $checkFailed = 1; }
        if (!$opt{FINALIZE_MEM}) { say "ERROR: No FINALIZE_MEM found in .ini file"; $checkFailed = 1; }
        if (!$opt{FINALIZE_TIME}) { say "ERROR: No FINALIZE_TIME found in .ini file"; $checkFailed = 1; }
    }

    die "One or more options not found or invalid in config files" if $checkFailed;
    $opt->{RUN_NAME} = basename($opt{OUTPUT_DIR});
    return;
}

sub lockRun {
    my ($dir) = @_;
    my $lock_file = catfile($dir, "run.lock");
    sysopen my $fh, $lock_file, O_WRONLY | O_CREAT | O_EXCL or die "Couldn't obtain lock file, are you *sure* there are no more jobs running? (error: $!)";
    close $fh;
    return;
}

sub recordGitVersion {
    my ($opt) = @_;

    my $git_dir = catfile(dirname(abs_path($0)), ".git");
    $opt->{VERSION} = `git --git-dir $git_dir describe --tags`;
    chomp $opt->{VERSION};
    return;
}

sub copyConfigAndScripts {
    my ($opt) = @_;

    my $pipeline_path = dirname(abs_path($0));
    my $slice_dir = catfile($pipeline_path, "settings", "slicing");
    my $strelka_dir = catfile($pipeline_path, "settings", "strelka");
    my $script_dir = catfile($pipeline_path, "scripts");
    my $qscript_dir = catfile($pipeline_path, "QScripts");

    rcopy $slice_dir, catfile($opt->{OUTPUT_DIR}, "settings", "slicing") or die "Failed to copy slice settings $slice_dir: $!";
    rcopy $strelka_dir, catfile($opt->{OUTPUT_DIR}, "settings", "strelka") or die "Failed to copy Strelka settings $strelka_dir: $!";
    rcopy $script_dir, catfile($opt->{OUTPUT_DIR}, "scripts") or die "Failed to copy scripts $script_dir: $!";
    rcopy $qscript_dir, catfile($opt->{OUTPUT_DIR}, "QScripts") or die "Failed to copy QScripts $qscript_dir: $!";
    foreach my $ini_file (@{$opt->{INIFILE}}) {
        rcopy $ini_file, catfile($opt->{OUTPUT_DIR}, "logs") or die "Failed to copy INI file $ini_file: $!";
        rcopy $ini_file, catfile($opt->{OUTPUT_DIR}, "settings") or die "Failed to copy INI file $ini_file: $!";
    }

    my $final_ini = catfile($opt->{OUTPUT_DIR}, "logs", "final.ini");
    open my $fh, ">", $final_ini or die "Couldn't open $final_ini: $!";
    say $fh join "\n", map { "$_\t$opt->{$_}" if defined $opt->{$_} } keys %{$opt};
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
