package HMF::Pipeline::Config;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Copy::Recursive qw(rcopy);
use File::Path qw(make_path);
use File::Spec::Functions;
use FindBin;
use Getopt::Long;
use IO::Pipe;
use POSIX qw(strftime);
use Time::HiRes qw(gettimeofday);

use HMF::Pipeline::Config::Validate qw(verifyConfig verifyBam);

use parent qw(Exporter);
our @EXPORT_OK = qw(
    parse
    validate
    createDirs
    addSubDir
    setupLogging
    addSamples
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
    my $out_dir = catfile($dirs->{out}, $dir);
    makePaths($out_dir);
    return $out_dir;
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
    ## no critic (Modules::RequireExplicitInclusion)
    STDOUT->autoflush(1);
    STDERR->autoflush(1);
    $out_fh->autoflush(1);
    ## use critic
    return;
}

sub addSamples {
    my ($opt) = @_;

    $opt->{SAMPLES} = {};

    if ($opt->{FASTQ}) {
        foreach my $input_file (sort keys %{$opt->{FASTQ}}) {
            my $fastqFile = fileparse($input_file);
            my ($sampleName) = split "_", $fastqFile;
            $opt->{SAMPLES}->{$sampleName} = $input_file;
            @{$opt->{RUNNING_JOBS}->{$sampleName}} = ();
        }
    }

    if ($opt->{BAM}) {
        foreach my $input_file (sort keys %{$opt->{BAM}}) {
            my $sampleName = verifyBam($input_file, $opt);
            not exists $opt->{SAMPLES}->{$sampleName} or die "sample '$sampleName' from $input_file already used by $opt->{SAMPLES}->{$sampleName}";
            $opt->{SAMPLES}->{$sampleName} = $input_file;
            @{$opt->{RUNNING_JOBS}->{$sampleName}} = ();
        }
    }
    return;
}

sub recordGitVersion {
    my ($opt) = @_;

    my $git_dir = catfile(pipelinePath(), ".git");
    $opt->{VERSION} = `git --git-dir $git_dir describe --tags`;
    chomp $opt->{VERSION};
    return;
}

sub copyConfigAndScripts {
    my ($opt) = @_;

    my $pipeline_path = pipelinePath();
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
    say $fh join "\n", map { "$_\t$opt->{$_}" } grep { defined $opt->{$_} and not ref $opt->{$_} } sort keys %{$opt};
    close $fh;

    return;
}

# do NOT depend on this from jobs
sub pipelinePath {
    return catfile($FindBin::Bin, updir());
}

1;
