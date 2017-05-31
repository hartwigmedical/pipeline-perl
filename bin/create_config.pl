#!/usr/bin/env perl

use FindBin::libs;
use discipline;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
use Getopt::Long;
use File::Path qw(make_path);
use File::Spec::Functions;
use File::Find::Rule;

use HMF::Pipeline::Template qw(writeFromTemplate);


my $settingsDir = catfile(dirname(abs_path($0)), updir(), "settings");
exit interactive() if @ARGV == 0;

GetOptions(
    "iniFile|i=s" => \my $iniFile,
    "iniPath|ip=s" => \my $iniPath,
    "outputDir|o=s" => \my $outputDir,
    "fastqDir|f=s" => \my @fastqDirs,
    "bamDir|b=s" => \my @bamDirs,
    "vcfFile|v=s" => \my $vcfFile,
    "mail|m=s" => \my $mail,
    "help|h" => \my $help,
    "run" => \my $run,
    )
    or warn usage()
    and exit 1;

say usage() and exit if $help;
warn usage() and exit 1 if not($iniFile or $iniPath) or not $outputDir or not(@fastqDirs or @bamDirs or $vcfFile) or not $mail;

my $ini;
if ($iniFile) {
    $ini = catfile($settingsDir, $iniFile);
} elsif ($iniPath) {
    $ini = $iniPath;
}
say "ERROR: $ini does not exist." if !-f $ini;

exit createConfig($ini, $outputDir, \@fastqDirs, \@bamDirs, $vcfFile, $mail, $run);


sub getIniFiles {
    my ($iniDir) = @_;

    -d $iniDir or die "Can't get INI files from $iniDir: $!";
    my @iniFiles = File::Find::Rule->file() #
        ->name("*.ini")                     #
        ->in($iniDir);                      #
    return \@iniFiles;
}

sub createConfig {
    my ($iniFile, $outputDir, $fastqDirs, $bamDirs, $vcfFile, $mail, $run) = @_;

    make_path($outputDir);

    map { die "$_ does not exist" if not -d } @{$fastqDirs};
    map { die "$_ does not exist" if not -d } @{$bamDirs};
    die "$vcfFile does not exist" if $vcfFile and not -f $vcfFile;

    my @fastqFiles = File::Find::Rule->file() #
        ->extras({follow => 1})               #
        ->name("*.fastq.gz")                  #
        ->in(@{$fastqDirs});                  #
    my @bamFiles = File::Find::Rule->file()   #
        ->extras({follow => 1})               #
        ->name("*.bam")                       #
        ->in(@{$bamDirs});                    #

    my $configFile = catfile($outputDir, "settings.config");
    writeFromTemplate(
        "Config.tt", $configFile,
        iniFile => $iniFile,
        outputDir => $outputDir,
        mail => $mail,
        fastqFiles => \@fastqFiles,
        bamFiles => \@bamFiles,
        vcfFiles => [$vcfFile],
    );

    if ($run) {
        my $pipeline = catfile(dirname(abs_path($0)), "pipeline.pl");
        return system "$pipeline $configFile";
    }
    return 0;
}

sub usage {
    my @iniFiles = @{getIniFiles($settingsDir)};
    my @usage = (
        "Usage: $0",
        "",
        "Advanced usage:",
        "$0 REQUIRED_ARGUMENTS [-run]",
        "",
        "Required INI file:",
        "",
        "\t-i, --inifile settings.ini",
        "\t-ip, --inipath /path/to/settings.ini",
        "",
        "Required input data:",
        "",
        "\t-f, --fastqdir /fastqFolder",
        "\t-b, --bamdir /bamFolder",
        "\t-v, --vcfFile vcfFile.vcf",
        "",
        "Required output config:",
        "",
        "\t-o, --outputdir /path/to/outputDir",
        "\t-m, --mail example\@mail.nl",
        "",
        "Available INI files: (use -i)",

        map { "\t$_:\t$iniFiles[$_]" } 0 .. $#iniFiles,
    );
    return join "\n", @usage;
}

sub interactive {
    say "Using interactive mode";
    say "Available INI files:";

    my @iniFiles = @{getIniFiles($settingsDir)};
    while (my ($iniIndex, $iniFile) = each @iniFiles) {
        say "\t${iniIndex}:\t${iniFile}";
    }

    print "Choose setting file [index]: ";
    chomp(my $iniIndex = <>);
    die "Please provide a correct INI index number" unless $iniIndex and $iniFiles[$iniIndex];
    my $iniFile = catfile($settingsDir, "$iniFiles[$iniIndex]");

    print "Output dir: ";
    chomp($outputDir = <>);
    die "Please provide a correct output directory" unless $outputDir;

    print "Input FASTQ data dir: ";
    chomp(my $rawDataDir = <>);
    push @fastqDirs, $rawDataDir if $rawDataDir;

    print "Input BAM data dir: ";
    chomp($rawDataDir = <>);
    push @bamDirs, $rawDataDir if $rawDataDir;

    die "Please provide a correct input data directory" unless @fastqDirs or @bamDirs;

    print "Input VCF: ";
    chomp($vcfFile = <>);
    $vcfFile = undef unless $vcfFile;

    print "Mail address: ";
    chomp($mail = <>);
    die "Please provide a correct mail address" unless $mail;

    return createConfig($iniFile, $outputDir, \@fastqDirs, \@bamDirs, $vcfFile, $mail, undef);
}
