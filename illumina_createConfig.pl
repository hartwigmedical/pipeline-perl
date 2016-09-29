#!/usr/bin/env perl

use 5.16.0;
use strict;
use warnings;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
use Getopt::Long;
use File::Path qw(make_path);
use File::Spec::Functions;
use File::Find::Rule;

use FindBin;
use lib "$FindBin::Bin";

use illumina_template;


my $settingsDir = catfile(dirname(abs_path($0)), "settings");

GetOptions("iniFile|i=s" => \my $iniFile,
           "iniPath|ip=s" => \my $iniPath,
           "outputDir|o=s" => \my $outputDir,
           "fastqDir|f=s" => \my @fastqDirs,
           "bamDir|b=s" => \my @bamDirs,
           "vcfFile|v=s" => \my $vcfFile,
           "mail|m=s" => \my $mail,
           "help|h" => \my $help,
           "run" => \my $run
          )
or die usage();

usage() if $help || !($iniFile || $iniPath) || !$outputDir || !(@fastqDirs || @bamDirs || $vcfFile) || !$mail;

my $ini;
if ($iniFile) {
    $ini = catfile($settingsDir, $iniFile);
} elsif ($iniPath) {
    $ini = $iniPath;
}
say "ERROR: $ini does not exist." if !-f $ini;

createConfig($ini, $outputDir, \@fastqDirs, \@bamDirs, $vcfFile, $mail, $run);


sub getIniFiles {
    my ($iniDir) = @_;

    -d $iniDir or die "Can't get INI files from $iniDir: $!";
    my @iniFiles = File::Find::Rule->file()
        ->name("*.ini")
        ->in($iniDir);
    while (my ($iniIndex, $iniFile) = each @iniFiles) {
        say "\t${iniIndex}:\t${iniFile}";
    }
    return \@iniFiles;
}

sub createConfig {
    my ($iniFile, $outputDir, $fastqDirs, $bamDirs, $vcfFile, $mail, $run) = @_;

    if (!-d $outputDir) {
        make_path($outputDir) or die "Couldn't create directory $outputDir: $!";
    }

    map { die "$_ does not exist" if !-d } @{$fastqDirs};
    map { die "$_ does not exist" if !-d } @{$bamDirs};
    die "$vcfFile does not exist" if $vcfFile and !-f $vcfFile;

    my @fastqFiles = File::Find::Rule->file()
        ->name("*.fastq.gz")
        ->in(@{$fastqDirs});
    my @bamFiles = File::Find::Rule->file()
        ->name("*.bam")
        ->in(@{$bamDirs});

    my $configFile = catfile($outputDir, "settings.config");
    from_template("Config.tt", $configFile,
                  iniFile => $iniFile,
                  outputDir => $outputDir,
                  mail => $mail,
                  fastqFiles => \@fastqFiles,
                  bamFiles => \@bamFiles,
                  vcfFile => $vcfFile);

    if ($run) {
        my $pipeline = catfile(dirname(abs_path($0)), "illumina_pipeline.pl");
        system "$pipeline $configFile";
    }
}

sub usage {
    say "Usage: perl illumina_createConfig.pl";
    say "";
    say "Advanced usage:";
    say "illumina_createConfig.pl REQUIRED_ARGUMENTS [-run]";
    say "";
    say "Required INI file:";
    say "";
    say "	-i, --inifile settings.ini";
    say "	-ip, --inipath /path/to/settings.ini";
    say "";
    say "Required input data:";
    say "";
    say "	-f, --fastqdir /fastqFolder";
    say "	-b, --bamdir /bamFolder";
    say "	-v, --vcfFile vcfFile.vcf";
    say "";
    say "Required output config:";
    say "";
    say "	-o, --outputdir /path/to/outputDir";
    say "	-m, --mail example\@mail.nl";
    say "";
    say "Available INI files: (use -i)";
    getIniFiles($settingsDir);
    exit;
}

sub interactive {
    say "Using interactive mode";
    say "Available INI files:";
    my @iniFiles = getIniFiles($settingsDir);

    print "Choose setting file [index]: ";
    chomp(my $iniIndex = <STDIN>);
    die "Please provide a correct INI index number" unless $iniIndex and $iniFiles[$iniIndex];
    my $iniFile = catfile($settingsDir, "$iniFiles[$iniIndex]");

    print "Output dir: ";
    chomp($outputDir = <STDIN>);
    die "Please provide a correct output directory" unless $outputDir;

    print "Raw data project dir: ";
    chomp(my $rawDataDir = <STDIN>);
    die "Please provide a correct raw data directory" unless $rawDataDir;
    push @fastqDirs, $rawDataDir;

    print "Mail address: ";
    chomp($mail = <STDIN>);
    die "Please provide a correct mail address" unless $mail;

    createConfig($iniFile, $outputDir, \@fastqDirs, $mail, $run);
}
