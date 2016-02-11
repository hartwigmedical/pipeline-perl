#!/usr/bin/perl -w
###################################################
### illumina_createConfig.pl
### - Creates a config file based on user input
### Author: R.F.Ernst
###################################################

use strict;
use Cwd            qw( abs_path );
use File::Basename qw( dirname );
use Getopt::Long;
use File::Path qw(make_path);

### Variables ###
my $settingsDir = dirname(abs_path($0))."/settings";

### Check usage ###
#interactive() if @ARGV == 0;

### Parse options ###
my $iniFile;
my $iniPath;
my $ini;
my $outputDir;
my @fastqDirs;
my @bamDirs;
my $vcfFile;
my $mail;
my $help;
my $run;

GetOptions ("iniFile|i=s" => \$iniFile,
	    "iniPath|ip=s" => \$iniPath,
	    "outputDir|o=s" => \$outputDir,
	    "fastqDir|f=s" => \@fastqDirs,
	    "bamDir|b=s" => \@bamDirs,
	    "vcfFile|v=s" => \$vcfFile,
	    "mail|m=s" => \$mail,
	    "help|h" => \$help,
	    "run" => \$run)
or die usage();

if ($help || ! ( $iniFile || $iniPath ) || ! $outputDir || ! (@fastqDirs || @bamDirs || $vcfFile ) || ! $mail ) { usage() }

### Non interactive mode ###
if ( $iniFile) {
    $ini = $settingsDir."/".$iniFile;
} elsif ($iniPath) {
    $ini = $iniPath;
}

if ( ! -e $ini ) {
    print "ERROR: $ini does not exist.\n";
}
createConfig($ini,$outputDir,\@fastqDirs,\@bamDirs,$vcfFile,$mail,$run);

### Parse and print available ini files ###
sub getIniFiles{
    my $iniDir = shift;
    my @iniFiles;
    my $iniIndex = -1;

    opendir (INIDIR, $iniDir) or die "Can't open $iniDir";
    while (my $iniFile = readdir(INIDIR)) {
	next unless ($iniFile =~ /\.ini$/); #skip non .ini files
	push(@iniFiles, $iniFile);
	$iniIndex ++;
	print "\t$iniIndex: \t $iniFile\n";
    }
    closedir(INIDIR);

    return(@iniFiles);
}

### Create config file ###
sub createConfig {
    my $iniFile = $_[0];
    my $outputDir = $_[1];
    my @fastqDirs = @{$_[2]};
    my @bamDirs = @{$_[3]};
    my $vcfFile = $_[4];
    my $mail = $_[5];
    my $run = $_[6];

    my $configFile = $outputDir."/settings.config";

    if(! -e $outputDir){
	make_path($outputDir) or die "Couldn't create directory: $outputDir\n";
    }
    # Create settings.config file in outputDir
    open CONFIG, ">$configFile" or die "cannot open file $configFile \n";
    print CONFIG "### SETTINGS ###\n";
    print CONFIG "INIFILE\t$iniFile\n";
    print CONFIG "OUTPUT_DIR\t$outputDir\n";
    print CONFIG "MAIL\t$mail\n";

    if(@fastqDirs){
	print CONFIG "\n### FASTQ FILES ###\n";
	#Find fastq files for each rawDataDir
	foreach my $fastqDir (@fastqDirs){
	    if(! -e $fastqDir) { die "$fastqDir does not exist." }
	    print CONFIG "# $fastqDir\n";
	    my $query_one = $fastqDir."/*.fastq.gz "; # look in fastq folder
	    my $query_two = $fastqDir."/*/*.fastq.gz "; # look one folder deep
	    my $query_three = $fastqDir."/*/*/*.fastq.gz "; # look two folders deep
	    my @fastqFiles = glob($query_one.$query_two.$query_three);
	    foreach my $fastqFile (@fastqFiles){ print CONFIG "FASTQ\t$fastqFile\n" }
	}
    }

    if(@bamDirs){
	print CONFIG "\n### BAM FILES###\n";
	foreach my $bamDir (@bamDirs){
	    if(! -e $bamDir) { die "$bamDir does not exist." }
	    print CONFIG "# $bamDir\n";
	    my $query_one = $bamDir."/*.bam "; # look in bam folder
	    my $query_two = $bamDir."/*/*.bam "; # look one folder deep
	    my @bamFiles = glob($query_one.$query_two);
	    foreach my $bamFile (@bamFiles){
		my $baiFile = $bamFile;
		$baiFile =~ s/\.bam/.bai/;
		if (! (-e $baiFile || -e "$bamFile.bai")) { die "ERROR: $bamFile.bai or $baiFile does not exist please create an index for $bamFile" }
		print CONFIG "BAM\t$bamFile\n"
	    }
	}
    }

    if($vcfFile){
	print CONFIG "\n### VCF FILE###\n";
	if(! -e $vcfFile) { die "$vcfFile does not exist." }
	print CONFIG "VCF\t$vcfFile\n"
    }

    close CONFIG;

    ###Run pipeline if -run is specified 
    if($run) {
	### run pipeline
	my $pipeline = dirname(abs_path($0))."/illumina_pipeline.pl";
	system "perl $pipeline $configFile";
    }
}

### Help information ###
sub usage{
    print "Usage: perl illumina_createConfig.pl\n\n";
    print "Advanced usage: \n";
    print "illumina_createConfig.pl (-i|-iniFile settings.ini OR -ip|-iniPath /path/to/settings.ini) -o|-outputDir /path/to/outputDir (-f|-fastqDir /fastqFolder OR -b|-bamDir /bamFolder OR -v|-vcfFile vcfFile.vcf) -m|-mail example\@mail.nl [-run]\n\n";
    print "Available ini files: (use -i)\n";
    getIniFiles($settingsDir);
    exit;
}

### Interactive mode ###
sub interactive{
    print "Using interactive mode \n";
    print "Avaible setting files:\n";
    my @iniFiles = getIniFiles($settingsDir);

    # Settings file
    print "Choose setting file [index]: ";
    chomp(my $iniIndex = <STDIN>);
    if ($iniIndex eq "" || ! $iniFiles[$iniIndex] ) { die "Please provide a correct ini index number." }
    my $iniFile = $settingsDir ."/". $iniFiles[$iniIndex];

    #Output dir
    print "Output dir: ";
    chomp($outputDir = <STDIN>); # no tab completion
    if ($outputDir eq "") { die "Please provide a correct output directory." }

    #Raw data dir -> add while loop to allow for multiple raw data dirs.
    print "Raw data project dir: ";
    chomp(my $rawDataDir = <STDIN>); # no tab completion
    if($rawDataDir eq "") { die "Please provide a correct raw data directory." } #check for existence
    push(@fastqDirs, $rawDataDir);

    #Output dir
    print "Mail address: ";
    chomp($mail = <STDIN>);
    if ($mail eq "") { die "Please provide a correct mail address." }

    #Create config
    createConfig($iniFile,$outputDir,\@fastqDirs,$mail,$run);
}