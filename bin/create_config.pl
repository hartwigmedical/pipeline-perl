#!/usr/bin/env perl

use FindBin::libs;
use discipline;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
use Getopt::Long;
use File::Path qw(make_path);
use File::Spec::Functions;
use File::Find::Rule;

use HMF::Pipeline::Functions::Template qw(writeFromTemplate);

my $settingsDir = catfile(dirname(abs_path($0)), updir(), "settings");

GetOptions(
    "iniFile|i=s" => \my $iniFile,
    "outputDir|o=s" => \my $outputDir,
    "fastqDir|f=s" => \my @fastqDirs,
    "bamDir|b=s" => \my @bamDirs,
    "help|h" => \my $help,
    )
    or warn usage()
    and exit 1;

say usage() and exit if $help;
warn usage() and exit 1 if not $iniFile or not $outputDir or not(@fastqDirs or @bamDirs);

my $ini = catfile($settingsDir, $iniFile);
say "ERROR: $ini does not exist." if !-f $ini;

exit createConfig($ini, $outputDir, \@fastqDirs, \@bamDirs);

sub createConfig {
    my ($iniFile, $outputDir, $fastqDirs, $bamDirs) = @_;

    make_path($outputDir);

    map { die "$_ does not exist" if not -d } @{$fastqDirs};
    map { die "$_ does not exist" if not -d } @{$bamDirs};

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
        fastqFiles => \@fastqFiles,
        bamFiles => \@bamFiles
    );

    return 0;
}

sub usage {
    my @iniFiles = @{getIniFiles($settingsDir)};
    my @usage = (
        "Usage: $0",
        "",
        "Required INI file:",
        "",
        "\t-i, --inifile settings.ini",
        "",
        "Required input data:",
        "",
        "\t-f, --fastqdir /fastqFolder",
        "\t-b, --bamdir /bamFolder",
        "",
        "Required output config:",
        "",
        "\t-o, --outputdir /path/to/outputDir",
        "",
        "Available INI files: (use -i)",

        map { "\t$_:\t$iniFiles[$_]" } 0 .. $#iniFiles,
    );
    return join "\n", @usage;
}

sub getIniFiles {
    my ($iniDir) = @_;

    -d $iniDir or die "Can't get INI files from $iniDir: $!";
    my @iniFiles = File::Find::Rule->file() #
        ->name("*.ini")                     #
        ->maxdepth(1)                       #
        ->in($iniDir);                      #

    return \@iniFiles;
}