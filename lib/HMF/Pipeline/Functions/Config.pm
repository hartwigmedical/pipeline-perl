package HMF::Pipeline::Functions::Config;

use FindBin::libs;
use discipline;

use Carp;
use File::Basename;
use File::Copy::Recursive qw(rcopy);
use File::Path qw(make_path);
use File::Spec::Functions;
use FindBin;
use Getopt::Long;
use IO::Pipe;
use List::MoreUtils qw(uniq);
use POSIX qw(strftime);
use Time::HiRes qw(gettimeofday);

use HMF::Pipeline::Functions::Validate qw(parseFastqName verifyConfig verifyBam);
use HMF::Pipeline::Functions::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(
    parse
    validate
    createDirs
    addSubDir
    setupLogging
    addSamples
    sampleBamAndJobs
    refSampleBamAndJobs
    sampleBamsAndJobs
    sampleControlBamsAndJobs
    allRunningJobs
    recordGitVersion
    copyConfigAndScripts
);

sub parse {
    my ($configurationFile, $opt) = @_;

    open my $fh, "<", $configurationFile or die "Couldn't open $configurationFile: $!";
    parseFile($fh, $opt);
    close $fh;
    return;
}

sub parseFile {
    my ($fh, $opt) = @_;

    while (<$fh>) {
        chomp;
        next if m/^#/ or not $_;
        my ($key, $val) = split "\t", $_, 2;
        die "Key '$key' is missing a value - is it badly formatted?" if not defined $val;

        if ($key eq 'INIFILE') {
            $val = catfile(pipelinePath(), $val) unless file_name_is_absolute($val);
            push @{$opt->{$key}}, $val;
            parse($val, $opt);
        } elsif ($key eq 'FASTQ' or $key eq 'BAM') {
            $opt->{$key}->{$val} = 1;
        } else {
            $opt->{$key} = $val;
        }
    }
    return;
}

sub validate {
    my ($opt) = @_;

    my @errors = @{verifyConfig($opt)};
    if (@errors) {
        foreach my $error (@errors) {
            warn "ERROR: $error";
        }
        die "One or more options not found or invalid in config files";
    }
    $opt->{RUN_NAME} = basename($opt->{OUTPUT_DIR});
    return;
}

sub makePaths {
    my (@dirs) = @_;

    make_path(@dirs, {error => \my $errors});
    my $messages = join ", ", map { join ": ", each %{$_} } @{$errors};
    die "Couldn't create directories: $messages" if $messages;
    return;
}

sub createDirs {
    my ($output_dir, %extra_dirs) = @_;

    foreach (@extra_dirs{keys %extra_dirs}) { $_ = catfile($output_dir, $_) }

    my $dirs = {
        out => $output_dir,
        tmp => catfile($output_dir, "tmp"),
        log => catfile($output_dir, "logs"),
        job => catfile($output_dir, "jobs"),
        %extra_dirs,
    };
    makePaths(values %{$dirs});
    return $dirs;
}

sub addSubDir {
    my ($dirs, $dir) = @_;
    my $output_dir = catfile($dirs->{out}, $dir);
    makePaths($output_dir);
    return $output_dir;
}

sub setupLogging {
    my ($output_dir) = @_;

    my ($seconds, $microseconds) = gettimeofday;
    my $datetime = strftime('%Y%m%d_%H%M%S_', localtime $seconds) . sprintf('%.6d', $microseconds);
    my $output_file = catfile($output_dir, "logs", "submitlog_${datetime}.out");
    my $error_file = catfile($output_dir, "logs", "submitlog_${datetime}.err");
    my $output_fh = IO::Pipe->new()->writer("tee $output_file") or die "Couldn't tee to $output_file: $!";
    my $error_fh = IO::Pipe->new()->writer("tee $error_file >&2") or die "Couldn't tee to $error_file: $!";
    open STDOUT, ">&", $output_fh or die "STDOUT redirection failed: $!";
    open STDERR, ">&", $error_fh or die "STDERR redirection failed: $!";
    ## no critic (Modules::RequireExplicitInclusion)
    STDOUT->autoflush(1);
    STDERR->autoflush(1);
    $output_fh->autoflush(1);
    ## use critic
    return;
}

sub addSamples {
    my ($opt) = @_;

    $opt->{SAMPLES} = {};

    if ($opt->{FASTQ}) {
        foreach my $input_path (sort keys %{$opt->{FASTQ}}) {
            my $sample_name = parseFastqName($input_path)->{sampleName};
            $opt->{SAMPLES}->{$sample_name} = [] if not exists $opt->{SAMPLES}->{$sample_name};
            push @{$opt->{SAMPLES}->{$sample_name}}, $input_path;
            @{$opt->{RUNNING_JOBS}->{$sample_name}} = ();
        }
    }

    if ($opt->{BAM}) {
        foreach my $input_path (sort keys %{$opt->{BAM}}) {
            my $sample_name = verifyBam($input_path, $opt);
            not exists $opt->{SAMPLES}->{$sample_name} or die "sample '$sample_name' from $input_path already used by $opt->{SAMPLES}->{$sample_name}[0]";
            $opt->{SAMPLES}->{$sample_name} = [$input_path];
            @{$opt->{RUNNING_JOBS}->{$sample_name}} = ();
        }
    }
    return;
}

sub sampleBamAndJobs {
    my ($sample, $opt) = @_;

    my $bam = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $opt->{BAM_FILES}->{$sample});

    return ($bam, $opt->{RUNNING_JOBS}->{$sample});
}

sub sampleBamsAndJobs {
    my ($opt) = @_;

    my $all_bams = {};
    my @all_jobs;
    foreach my $sample (sort keys %{$opt->{SAMPLES}}) {
        my ($bam, $jobs) = sampleBamAndJobs($sample, $opt);
        $all_bams->{$sample} = $bam;
        push @all_jobs, @{$jobs};
    }
    return ($all_bams, [ uniq @all_jobs ]);
}

sub refSampleBamAndJobs {
    my ($opt) = @_;
    my ($ref_sample) = HMF::Pipeline::Functions::Metadata::refSampleName($opt);
    $opt->{BAM_FILES}->{$ref_sample} or die "metadata ref_sample $ref_sample not in BAM file list: " . join ", ", keys %{$opt->{BAM_FILES}};
    my ($ref_sample_bam, $ref_sample_jobs) = sampleBamAndJobs($ref_sample, $opt);
    return ($ref_sample, $ref_sample_bam, [ uniq @{$ref_sample_jobs} ]);
}

sub sampleControlBamsAndJobs {
    my ($opt) = @_;

    my ($ref_sample, $tumor_sample, $joint_name) = HMF::Pipeline::Functions::Metadata::sampleControlNames($opt);

    $opt->{BAM_FILES}->{$ref_sample} or die "metadata ref_sample $ref_sample not in BAM file list: " . join ", ", keys %{$opt->{BAM_FILES}};
    $opt->{BAM_FILES}->{$tumor_sample} or die "metadata tumor_sample $tumor_sample not in BAM file list: " . join ", ", keys %{$opt->{BAM_FILES}};

    my ($ref_sample_bam, $ref_sample_jobs) = sampleBamAndJobs($ref_sample, $opt);
    my ($tumor_sample_bam, $tumor_sample_jobs) = sampleBamAndJobs($tumor_sample, $opt);

    return ($ref_sample, $tumor_sample, $ref_sample_bam, $tumor_sample_bam, $joint_name, [ uniq @{$ref_sample_jobs}, @{$tumor_sample_jobs} ]);
}

sub allRunningJobs {
    my ($opt) = @_;

    my @running_job_chains = grep { defined } values %{$opt->{RUNNING_JOBS}};
    my @running_jobs = map { @$_ } @running_job_chains;
    my @unique_jobs = uniq sort grep { defined } @running_jobs;
    return \@unique_jobs;
}

sub recordGitVersion {
    my ($opt) = @_;

    my $git_dir = catfile(pipelinePath(), ".git");
    $opt->{VERSION} = qx(git --git-dir $git_dir describe --tags);
    chomp $opt->{VERSION};
    return;
}

sub copyConfigAndScripts {
    my ($opt) = @_;

    my $pipeline_path = pipelinePath();
    my $slice_dir = catfile($pipeline_path, "settings", "slicing");
    my $strelka_dir = catfile($pipeline_path, "settings", "strelka");
    my $gridss_dir = catfile($pipeline_path, "settings", "gridss");
    my $scripts_dir = catfile($pipeline_path, "scripts");
    my $qscript_dir = catfile($pipeline_path, "QScripts");

    rcopy $slice_dir, catfile($opt->{OUTPUT_DIR}, "settings", "slicing") or die "Failed to copy slice settings $slice_dir: $!";
    rcopy $strelka_dir, catfile($opt->{OUTPUT_DIR}, "settings", "strelka") or die "Failed to copy Strelka settings $strelka_dir: $!";
    rcopy $gridss_dir, catfile($opt->{OUTPUT_DIR}, "settings", "gridss") or die "Failed to copy Gridss settings $gridss_dir: $!";
    rcopy $scripts_dir, catfile($opt->{OUTPUT_DIR}, "scripts") or die "Failed to copy scripts directory $scripts_dir: $!";
    rcopy $qscript_dir, catfile($opt->{OUTPUT_DIR}, "QScripts") or die "Failed to copy QScripts $qscript_dir: $!";
    foreach my $ini_file (@{$opt->{INIFILE}}) {
        rcopy $ini_file, catfile($opt->{OUTPUT_DIR}, "logs") or die "Failed to copy INI file $ini_file: $!";
        rcopy $ini_file, catfile($opt->{OUTPUT_DIR}, "settings") or die "Failed to copy INI file $ini_file: $!";
    }

    my $final_ini = catfile($opt->{OUTPUT_DIR}, "logs", "final.ini");
    open my $fh, ">", $final_ini or die "Couldn't open $final_ini: $!";
    say $fh join "\n", map { "$_\t$opt->{$_}" } grep { defined $opt->{$_} and not ref $opt->{$_} } sort keys %{$opt};
    close $fh;

    return;
}

# SABR: do NOT depend on this from jobs
sub pipelinePath {
    return catfile($FindBin::Bin, updir());
}

1;
