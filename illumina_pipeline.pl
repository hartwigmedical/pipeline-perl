#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use File::Copy qw(copy);
use Cwd qw( abs_path );
use File::Basename qw( dirname );

use FindBin;
use lib "$FindBin::Bin";

use illumina_prestats;
use illumina_mapping;
use illumina_poststats;
use illumina_realign;
use illumina_germlineCalling;
use illumina_germlineFiltering;
use illumina_germlineAnnotation;
use illumina_somaticVariants;
use illumina_copyNumber;
use illumina_baf;
use illumina_kinship;
use illumina_finalize;

die usage() if @ARGV == 0;

my %opt;
my $configurationFile;

# KODU - Do not use in .conf or .ini
%opt = (
    'RUNNING_JOBS'		=> {},
    'BAM_FILES'			=> {},
    'SAMPLES'			=> undef,
    'PIPELINE_PATH'		=> dirname(abs_path($0)),
    'SINGLE_END'        => undef
);

############ READ RUN SETTINGS FORM .conf FILE ############
$configurationFile = $ARGV[0];

open (CONFIGURATION, "<$configurationFile") or die "Couldn't open .conf file: $configurationFile\n";
while(<CONFIGURATION>) {
    chomp;
    next if m/^#/ or ! $_;
    my ($key, $val) = split("\t",$_,2);

    if($key eq 'INIFILE') {
        $opt{$key} = $val;
        open (INI, "<$val") or die "Couldn't open .ini file $val\n";
        while(<INI>) {
            chomp;
            next if m/^#/ or ! $_;
            my ($key, $val) = split("\t",$_,2);
            $opt{$key} = $val;
        }
        close INI;
    } elsif ($key eq 'FASTQ' or $key eq 'BAM') {
        $opt{$key}->{$val} = 1;
    } else {
        $opt{$key} = $val;
    }
}
close CONFIGURATION;

############ START PIPELINE  ############

checkConfig();
getSamples();
createOutputDirs();

system "cp $opt{INIFILE} $opt{OUTPUT_DIR}/logs";

my $opt_ref;

if ($opt{FASTQ}) {
    if ($opt{PRESTATS} eq "yes") {
        print "### SCHEDULING PRESTATS ###\n";
        illumina_prestats::runPreStats(\%opt);
    }

    if ($opt{MAPPING} eq "yes") {
        print "\n### SCHEDULING MAPPING ###\n";
        $opt_ref = illumina_mapping::runMapping(\%opt);
        %opt = %$opt_ref;
    }
} elsif ($opt{BAM}) {
    $opt{MAPPING} = "no";
    $opt{PRESTATS} = "no";

    print "\n###SCHEDULING BAM PREP###\n";
    $opt_ref = illumina_mapping::runBamPrep(\%opt);
    %opt = %$opt_ref;
}

if ($opt{FASTQ} or $opt{BAM}) {
    if ($opt{POSTSTATS} eq "yes") {
        print "\n### SCHEDULING POSTSTATS ###\n";
        my $postStatsJob = illumina_poststats::runPostStats(\%opt);
        $opt{RUNNING_JOBS}->{'postStats'} = $postStatsJob;
    }

    if ($opt{INDELREALIGNMENT} eq "yes") {
        print "\n### SCHEDULING INDELREALIGNMENT ###\n";
        $opt_ref = illumina_realign::runRealignment(\%opt);
        %opt = %$opt_ref;
    }

    if ($opt{SOMATIC_VARIANTS} eq "yes") {
        print "\n### SCHEDULING SOMATIC VARIANT CALLERS ####\n";
        my $somVar_jobs = illumina_somaticVariants::runSomaticVariantCallers(\%opt);
        $opt{RUNNING_JOBS}->{'somVar'} = $somVar_jobs;
    }

    if ($opt{COPY_NUMBER} eq "yes") {
        print "\n### SCHEDULING COPY NUMBER TOOLS ####\n";
        my $cnv_jobs = illumina_copyNumber::runCopyNumberTools(\%opt);
        $opt{RUNNING_JOBS}->{'CNV'} = $cnv_jobs;
    }

    if ($opt{BAF} eq "yes") {
        print "\n### SCHEDULING BAF ANALYSIS ###\n";
        my $baf_jobs = illumina_baf::runBAF(\%opt);
        $opt{RUNNING_JOBS}->{'baf'} = $baf_jobs;
    }

    if ($opt{VARIANT_CALLING} eq "yes") {
        print "\n### SCHEDULING VARIANT CALLING ####\n";
        $opt_ref = illumina_germlineCalling::runVariantCalling(\%opt);
        %opt = %$opt_ref;
    }

    if ($opt{FILTER_VARIANTS} eq "yes") {
        print "\n### SCHEDULING VARIANT FILTRATION ####\n";
        my $FVJob = illumina_germlineFiltering::runFilterVariants(\%opt);

        foreach my $sample (@{$opt{SAMPLES}}){
            push (@{$opt{RUNNING_JOBS}->{$sample}} , $FVJob);
        }
    }

    if ($opt{ANNOTATE_VARIANTS} eq "yes") {
        print "\n### SCHEDULING VARIANT ANNOTATION ####\n";
        my $AVJob = illumina_germlineAnnotation::runAnnotateVariants(\%opt);

        foreach my $sample (@{$opt{SAMPLES}}){
            push (@{$opt{RUNNING_JOBS}->{$sample}} , $AVJob);
        }
    }

    if ($opt{KINSHIP} eq "yes") {
        print "\n### SCHEDULING KINSHIP ####\n";
        my $kinship_job = illumina_kinship::runKinship(\%opt);
        $opt{RUNNING_JOBS}->{'Kinship'} = $kinship_job;
    }

    if ($opt{FINALIZE} eq "yes") {
        print "\n### SCHEDULING PIPELINE FINALIZE ####\n";
        illumina_finalize::runFinalize(\%opt);
    }
}

############ SUBROUTINES  ############
sub getSamples{
    my %samples;

    if ($opt{FASTQ}) {
        foreach my $input (keys %{$opt{FASTQ}}){
            my $fastqFile = (split("/", $input))[-1];
            my $sampleName = (split("_", $fastqFile))[0];
            $samples{$sampleName}++;
            @{$opt{RUNNING_JOBS}->{$sampleName}} = ();
        }
    }

    if ($opt{BAM}) {
        foreach my $input (keys %{$opt{BAM}}){
            my $bamFile = (split("/", $input))[-1];
            my $sampleName = $bamFile;
            $sampleName =~ s/\.bam//g;
            $samples{$sampleName}++;
            @{$opt{RUNNING_JOBS}->{$sampleName}} = ();
        }
    }

    @{$opt{SAMPLES}} = keys(%samples);
}

sub createOutputDirs{
    if (! -e $opt{OUTPUT_DIR}) {
	    make_path($opt{OUTPUT_DIR}) or die "Couldn't create directory: $opt{OUTPUT_DIR}\n";
    }

    if (! -e "$opt{OUTPUT_DIR}/QCStats") {
	    mkdir("$opt{OUTPUT_DIR}/QCStats") or die "Couldn't create directory: $opt{OUTPUT_DIR}/QCStats\n";
    }

    if (! -e "$opt{OUTPUT_DIR}/jobs") {
	    mkdir("$opt{OUTPUT_DIR}/jobs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/jobs\n";
    }

    if (! -e "$opt{OUTPUT_DIR}/logs") {
	    mkdir("$opt{OUTPUT_DIR}/logs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/logs\n";
    }

    if (! -e "$opt{OUTPUT_DIR}/tmp") {
	    mkdir("$opt{OUTPUT_DIR}/tmp") or die "Couldn't create directory: $opt{OUTPUT_DIR}/tmp\n";
    }

    foreach my $sample (@{$opt{SAMPLES}}) {
        if (! -e "$opt{OUTPUT_DIR}/$sample"){
            mkdir("$opt{OUTPUT_DIR}/$sample") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample\n";
        }

        if (! -e "$opt{OUTPUT_DIR}/$sample/mapping") {
            mkdir("$opt{OUTPUT_DIR}/$sample/mapping") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/mapping\n";
        }

        if (! -e "$opt{OUTPUT_DIR}/$sample/QCStats") {
            mkdir("$opt{OUTPUT_DIR}/$sample/QCStats") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/QCStats\n";
        }

        if (! -e "$opt{OUTPUT_DIR}/$sample/jobs") {
            mkdir("$opt{OUTPUT_DIR}/$sample/jobs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/jobs\n";
        }

        if (! -e "$opt{OUTPUT_DIR}/$sample/logs") {
            mkdir("$opt{OUTPUT_DIR}/$sample/logs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/logs\n";
        }

        if (! -e "$opt{OUTPUT_DIR}/$sample/tmp") {
            mkdir("$opt{OUTPUT_DIR}/$sample/tmp") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/tmp\n";
        }
    }
}

sub usage{
    warn <<END;
    Usage: perl illumina_pipeline.pl configurationFile.conf
END
    exit;
}

sub checkConfig {
    my $checkFailed = 0;
    my $runName = "";
    ### Input and Output
    if(! $opt{INIFILE}){ print "ERROR: No INIFILE option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{OUTPUT_DIR}){ print "ERROR: No OUTPUT_DIR found in config files.\n"; $checkFailed = 1; } else { $runName = (split("/", $opt{OUTPUT_DIR}))[-1];}
    if(! ($opt{FASTQ} || $opt{BAM})){ print "ERROR: No FASTQ or BAM files found in config files.\n"; $checkFailed = 1; }
    if(! $opt{MAIL}){ print "ERROR: No MAIL address specified in config files.\n"; $checkFailed = 1; }

    ### Cluster settings
    if(! $opt{CLUSTER_PATH}){ print "ERROR: No CLUSTER_PATH option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{CLUSTER_TMP}){ print "ERROR: No CLUSTER_TMP option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{CLUSTER_RESERVATION}){ print "ERROR: No CLUSTER_RESERVATION option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{CLUSTER_PROJECT}){ print "ERROR: No CLUSTER_PROJECT option found in config files.\n"; $checkFailed = 1; }

    ### Module Yes or No
    if(! $opt{PRESTATS}){ print "ERROR: No PRESTATS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{MAPPING}){ print "ERROR: No MAPPING option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{POSTSTATS}){ print "ERROR: No POSTSTATS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{INDELREALIGNMENT}){ print "ERROR: No INDELREALIGNMENT option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{VARIANT_CALLING}){ print "ERROR: No VARIANT_CALLING option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{FILTER_VARIANTS}){ print "ERROR: No FILTER_VARIANTS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{SOMATIC_VARIANTS}){ print "ERROR: No SOMATIC_VARIANTS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{COPY_NUMBER}){ print "ERROR: No COPY_NUMBER option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{BAF}){ print "ERROR: No BAF option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{ANNOTATE_VARIANTS}){ print "ERROR: No ANNOTATE_VARIANTS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{KINSHIP}){ print "ERROR: No KINSHIP option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{FINALIZE}){ print "ERROR: No FINALIZE option found in config files. \n"; $checkFailed = 1; }

    ### Module Settings / tools
    if(! $opt{GENOME}){ print "ERROR: No GENOME option found in config files.\n"; $checkFailed = 1; }
    elsif(! -e $opt{GENOME}){ print"ERROR: $opt{GENOME} does not exist\n"}
    if(! $opt{SAMBAMBA_PATH}){ print "ERROR: No SAMBAMBA_PATH option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{QUEUE_PATH}){ print "ERROR: No QUEUE_PATH option found in config files.\n"; $checkFailed = 1; }

    ## PRESTATS
    if($opt{PRESTATS} && $opt{PRESTATS} eq "yes"){
        if(! $opt{FASTQC_PATH}){ print "ERROR: No FASTQC_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{PRESTATS_THREADS}){ print "ERROR: No PRESTATS_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{PRESTATS_MEM}){ print "ERROR: No PRESTATS_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{PRESTATS_QUEUE}){ print "ERROR: No PRESTATS_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{PRESTATS_TIME}){ print "ERROR: No PRESTATS_TIME option found in config files.\n"; $checkFailed = 1; }
    }

    ## MAPPING
    if ($opt{MAPPING} && $opt{MAPPING} eq "yes") {
        if(! $opt{BWA_PATH}){ print "ERROR: No BWA_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{MAPPING_THREADS}){ print "ERROR: No MAPPING_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{MAPPING_MEM}){ print "ERROR: No MAPPING_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{MAPPING_QUEUE}){ print "ERROR: No MAPPING_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{MAPPING_TIME}){ print "ERROR: No MAPPING_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{MAPPING_SETTINGS}){ print "ERROR: No MAPPING_SETTINGS option found in config files.\n"; $checkFailed = 1; }

        if(! $opt{MARKDUP_QUEUE}){ print "ERROR: No MARKDUP_QUEUE option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{MARKDUP_TIME}){ print "ERROR: No MARKDUP_TIME option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{MARKDUP_THREADS}){ print "ERROR: No MARKDUP_THREADS option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{MARKDUP_MEM}){ print "ERROR: No MARKDUP_MEM option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{MARKDUP_OVERFLOW_LIST_SIZE}){ print "ERROR: No MARKDUP_OVERFLOW_LIST_SIZE option found in config files.\n"; $checkFailed = 1; }
	}

	if(! $opt{FLAGSTAT_QUEUE}){ print "ERROR: No FLAGSTAT_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FLAGSTAT_THREADS}){ print "ERROR: No FLAGSTAT_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FLAGSTAT_MEM}){ print "ERROR: No FLAGSTAT_MEM option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FLAGSTAT_TIME}){ print "ERROR: No FLAGSTAT_TIME option found in config files.\n"; $checkFailed = 1; }

    ## POSTSTATS
    if($opt{POSTSTATS} && $opt{POSTSTATS} eq "yes"){
        if(! $opt{BAMMETRICS_PATH}){ print "ERROR: No BAMMETRICS_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{PICARD_PATH}){ print "ERROR: No PICARD_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{POSTSTATS_THREADS}){ print "ERROR: No POSTSTATS_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{POSTSTATS_MEM}){ print "ERROR: No POSTSTATS_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{POSTSTATS_QUEUE}){ print "ERROR: No POSTSTATS_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{POSTSTATS_TIME}){ print "ERROR: No POSTSTATS_TIME option found in config files.\n"; $checkFailed = 1; }
    }

    ## INDELREALIGNMENT
    if($opt{INDELREALIGNMENT} && $opt{INDELREALIGNMENT} eq "yes"){
        if(! $opt{BAMUTIL_PATH}){ print "ERROR: No BAMUTIL_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_MASTER_QUEUE}){ print "ERROR: No REALIGNMENT_MASTER_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_MASTER_THREADS}){ print "ERROR: No REALIGNMENT_MASTER_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_MASTER_TIME}){ print "ERROR: No REALIGNMENT_MASTER_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_MASTER_MEM}){ print "ERROR: No REALIGNMENT_MASTER_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_QUEUE}){ print "ERROR: No REALIGNMENT_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_THREADS}){ print "ERROR: No REALIGNMENT_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_MEM}){ print "ERROR: No REALIGNMENT_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_TIME}){ print "ERROR: No REALIGNMENT_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_MERGETHREADS}){ print "ERROR: No REALIGNMENT_MERGETHREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_SCALA}){ print "ERROR: No REALIGNMENT_SCALA option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{REALIGNMENT_SCATTER}){ print "ERROR: No REALIGNMENT_SCATTER option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FLAGSTAT_QUEUE}){ print "ERROR: No FLAGSTAT_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FLAGSTAT_THREADS}){ print "ERROR: No FLAGSTAT_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FLAGSTAT_MEM}){ print "ERROR: No FLAGSTAT_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FLAGSTAT_TIME}){ print "ERROR: No FLAGSTAT_TIME option found in config files.\n"; $checkFailed = 1; }
    }

    ## VARIANT_CALLING
    if($opt{VARIANT_CALLING} && $opt{VARIANT_CALLING} eq "yes"){
        if(! $opt{CALLING_MASTER_QUEUE}){ print "ERROR: No CALLING_MASTER_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_MASTER_TIME}){ print "ERROR: No CALLING_MASTER_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_MASTER_THREADS}){ print "ERROR: No CALLING_MASTER_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_MASTER_MEM}){ print "ERROR: No CALLING_MASTER_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_QUEUE}){ print "ERROR: No CALLING_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_THREADS}){ print "ERROR: No CALLING_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_MEM}){ print "ERROR: No CALLING_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_TIME}){ print "ERROR: No CALLING_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_SCATTER}){ print "ERROR: No CALLING_SCATTER option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_GVCF}){ print "ERROR: No CALLING_GVCF option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_SCALA}){ print "ERROR: No CALLING_SCALA option found in config files.\n"; $checkFailed = 1; }
        if($opt{CALLING_UGMODE}){
            if($opt{CALLING_UGMODE} ne "SNP" and $opt{CALLING_UGMODE} ne "INDEL" and $opt{CALLING_UGMODE} ne "BOTH"){ print "ERROR: UGMODE: $opt{CALLING_UGMODE} does not exist use SNP, INDEL or BOTH\n"; $checkFailed = 1; }
        }
        if(! $opt{CALLING_STANDCALLCONF}){ print "ERROR: No CALLING_STANDCALLCONF option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{CALLING_STANDEMITCONF}){ print "ERROR: No CALLING_STANDEMITCONF option found in config files.\n"; $checkFailed = 1; }
        if( $opt{CALLING_TARGETS} && ! -e $opt{CALLING_TARGETS}) { print"ERROR: $opt{CALLING_TARGETS} does not exist\n"; $checkFailed = 1; }
        if( $opt{CALLING_DBSNP} && ! -e $opt{CALLING_DBSNP}) { print"ERROR: $opt{CALLING_DBSNP} does not exist\n"; $checkFailed = 1; }
    }

    ## FILTER_VARIANTS
    if($opt{FILTER_VARIANTS} && $opt{FILTER_VARIANTS} eq "yes"){
        if(! $opt{FILTER_MASTER_QUEUE}){ print "ERROR: No FILTER_MASTER_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_MASTER_TIME}){ print "ERROR: No FILTER_MASTER_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_MASTER_THREADS}){ print "ERROR: No FILTER_MASTER_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_MASTER_MEM}){ print "ERROR: No FILTER_MASTER_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_QUEUE}){ print "ERROR: No FILTER_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_THREADS}){ print "ERROR: No FILTER_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_MEM}){ print "ERROR: No FILTER_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_TIME}){ print "ERROR: No FILTER_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_SCATTER}){ print "ERROR: No FILTER_SCATTER option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_SCALA}){ print "ERROR: No FILTER_SCALA option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_SNPTYPES}){ print "ERROR: No FILTER_SNPTYPES option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_SNPNAME}){ print "ERROR: No FILTER_SNPNAME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_SNPEXPR}){ print "ERROR: No FILTER_SNPEXPR  option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_INDELTYPES}){ print "ERROR: No FILTER_INDELTYPES option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_INDELNAME}){ print "ERROR: No FILTER_INDELNAME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{FILTER_INDELEXPR}){ print "ERROR: No FILTER_INDELEXPR option found in config files.\n"; $checkFailed = 1; }
    }

    ## SOMATIC_VARIANTS
    if($opt{SOMATIC_VARIANTS} && $opt{SOMATIC_VARIANTS} eq "yes") {
        if(! $opt{VCFTOOLS_PATH}){ print "ERROR: No VCFTOOLS_PATH found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{SAMTOOLS_PATH}){ print "ERROR: No SAMTOOLS_PATH option found in config files.\n"; $checkFailed = 1; }
        if( $opt{SOMVAR_TARGETS} && ! -e $opt{SOMVAR_TARGETS}) { print"ERROR: $opt{SOMVAR_TARGETS} does not exist\n"; $checkFailed = 1; }
        if(! $opt{SOMVAR_STRELKA}){ print "ERROR: No SOMVAR_STRELKA option found in config files.\n"; $checkFailed = 1; }
        if($opt{SOMVAR_STRELKA} && $opt{SOMVAR_STRELKA} eq "yes") {
            if(! $opt{STRELKA_PATH}){ print "ERROR: No STRELKA_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{STRELKA_INI}){ print "ERROR: No STRELKA_INI option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{STRELKA_QUEUE}){ print "ERROR: No STRELKA_QUEUE option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{STRELKA_THREADS}){ print "ERROR: No STRELKA_THREADS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{STRELKA_MEM}){ print "ERROR: No STRELKA_MEM option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{STRELKA_TIME}){ print "ERROR: No STRELKA_TIME option found in config files.\n"; $checkFailed = 1; }
        }
        if(! $opt{SOMVAR_VARSCAN}){ print "ERROR: No SOMVAR_VARSCAN option found in config files.\n"; $checkFailed = 1; }
        if($opt{SOMVAR_VARSCAN} && $opt{SOMVAR_VARSCAN} eq "yes") {
            if(! $opt{VARSCAN_PATH}){ print "ERROR: No VARSCAN_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{PBGZIP_PATH}){ print "ERROR: No PBGZIP_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{TABIX_PATH}){ print "ERROR: No TABIX_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{VARSCAN_QUEUE}){ print "ERROR: No VARSCAN_QUEUE option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{VARSCAN_THREADS}){ print "ERROR: No VARSCAN_THREADS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{VARSCAN_TIME}){ print "ERROR: No VARSCAN_TIME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{VARSCAN_MEM}){ print "ERROR: No VARSCAN_MEM option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{VARSCAN_SETTINGS}){ print "ERROR: No VARSCAN_SETTINGS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{VARSCAN_POSTSETTINGS}){ print "ERROR: No VARSCAN_POSTSETTINGS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{PILEUP_QUEUE}){ print "ERROR: No PILEUP_QUEUE option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{PILEUP_DIVISOR}){ print "ERROR: No PILEUP_DIVISOR option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{PILEUP_THREADS}){ print "ERROR: No PILEUP_THREADS option found in config files.\n"; $checkFailed = 1; }
            elsif( $opt{PILEUP_THREADS} < $opt{PILEUP_DIVISOR}){ print "ERROR: PILEUP_THREADS ($opt{PILEUP_THREADS}) must be at least PILEUP_DIVISOR ($opt{PILEUP_DIVISOR}).\n", $checkFailed = 1; }
            if(! $opt{PILEUP_MEM}){ print "ERROR: No PILEUP_MEM option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{PILEUP_TIME}){ print "ERROR: No PILEUP_TIME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FINALIZE_KEEP_PILEUP}){ print "ERROR: No FINALIZE_KEEP_PILEUP found in .ini file\n"; $checkFailed = 1; }
        }
        if(! $opt{SOMVAR_FREEBAYES}){ print "ERROR: No SOMVAR_FREEBAYES option found in config files.\n"; $checkFailed = 1; }
        if($opt{SOMVAR_FREEBAYES} && $opt{SOMVAR_FREEBAYES} eq "yes") {
            if(! $opt{FREEBAYES_PATH}){ print "ERROR: No FREEBAYES_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{VCFLIB_PATH}){ print "ERROR: No VCFLIB_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEBAYES_QUEUE}){ print "ERROR: No FREEBAYES_QUEUE option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEBAYES_THREADS}){ print "ERROR: No FREEBAYES_THREADS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEBAYES_MEM}){ print "ERROR: No FREEBAYES_MEM option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEBAYES_TIME}){ print "ERROR: No FREEBAYES_TIME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEBAYES_SETTINGS}){ print "ERROR: No FREEBAYES_SETTINGS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEBAYES_SOMATICFILTER}){ print "ERROR: No FREEBAYES_SOMATICFILTER option found in config files.\n"; $checkFailed = 1; }
        }
        if(! $opt{SOMVAR_MUTECT}){ print "ERROR: No SOMVAR_MUTECT option found in config files.\n"; $checkFailed = 1; }
        if($opt{SOMVAR_MUTECT} && $opt{SOMVAR_MUTECT} eq "yes") {
            if(! $opt{MUTECT_PATH}){ print "ERROR: No MUTECT_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{MUTECT_QUEUE}){ print "ERROR: No MUTECT_QUEUE option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{MUTECT_THREADS}){ print "ERROR: No MUTECT_THREADS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{MUTECT_MEM}){ print "ERROR: No MUTECT_MEM option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{MUTECT_TIME}){ print "ERROR: No MUTECT_TIME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{MUTECT_COSMIC}){ print "ERROR: No MUTECT_COSMIC option found in config files.\n"; $checkFailed = 1; }
        }
        if(! $opt{SOMVARMERGE_QUEUE}){ print "ERROR: No SOMVARMERGE_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{SOMVARMERGE_THREADS}){ print "ERROR: No SOMVARMERGE_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{SOMVARMERGE_MEM}){ print "ERROR: No SOMVARMERGE_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{SOMVARMERGE_TIME}){ print "ERROR: No SOMVARMERGE_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{SOMVAR_ANNOTATE}){ print "ERROR: No SOMVAR_ANNOTATE option found in config files.\n"; $checkFailed = 1; }
        if($opt{SOMVAR_ANNOTATE} && $opt{SOMVAR_ANNOTATE} eq "yes") {
            if(! $opt{ANNOTATE_DB}){ print "ERROR: No ANNOTATE_DB option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_FLAGS}){ print "ERROR: No ANNOTATE_FLAGS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_IDNAME}){ print "ERROR: No ANNOTATE_IDNAME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_IDDB}){ print "ERROR: No ANNOTATE_IDDB option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{CALLING_DBSNP}){ print "ERROR: No CALLING_DBSNP option found in config files.\n"; $checkFailed = 1; }
        }
    }

    ## COPY_NUMBER
    if($opt{COPY_NUMBER} && $opt{COPY_NUMBER} eq "yes"){
        if(! $opt{CNVCHECK_QUEUE} ) { print "ERROR: No CNVCHECK_QUEUE in config files.\n"; $checkFailed = 1; }
        if(! $opt{CNVCHECK_THREADS} ) { print "ERROR: No CNVCHECK_THREADS  in config files.\n"; $checkFailed = 1; }
        if(! $opt{CNVCHECK_MEM} ) { print "ERROR: No CNVCHECK_MEM in config files.\n"; $checkFailed = 1; }
        if(! $opt{CNVCHECK_TIME} ) { print "ERROR: No CNVCHECK_TIME in config files.\n"; $checkFailed = 1; }
        if(! $opt{CNV_MODE}){ print "ERROR: No CNV_MODE in config files. \n"; $checkFailed = 1; }
        if(! $opt{CNV_FREEC}){ print "ERROR: No CNV_FREEC  in config files.\n"; $checkFailed = 1; }
        if($opt{CNV_FREEC} eq "yes"){
            if(! $opt{FREEC_PATH}){ print "ERROR: No FREEC_PATH option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_QUEUE}){ print "ERROR: No FREEC_QUEUE option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_THREADS}){ print "ERROR: No FREEC_THREADS option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_MEM}){ print "ERROR: No FREEC_MEM option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_TIME}){ print "ERROR: No FREEC_TIME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_CHRLENFILE}){ print "ERROR: No FREEC_CHRLENFILE option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_CHRFILES}){ print "ERROR: No FREEC_CHRFILES option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_PLOIDY}){ print "ERROR: No FREEC_PLOIDY option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_WINDOW}){ print "ERROR: No FREEC_WINDOW option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{FREEC_TELOCENTROMERIC}){ print "ERROR: No FREEC_TELOCENTROMERIC option found in config files.\n"; $checkFailed = 1; }
        }
    }

    ##BAF Analysis
    if($opt{BAF} && $opt{BAF} eq "yes"){
        if(! $opt{BAF_QUEUE}){ print "ERROR: No BAF_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{BAF_THREADS}){ print "ERROR: No BAF_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{BAF_MEM}){ print "ERROR: No BAF_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{BAF_TIME}){ print "ERROR: No BAF_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{BIOVCF_PATH}){ print "ERROR: No BIOVCF_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{BAF_SNPS}){ print "ERROR: No BAF_SNPS option found in config files.\n"; $checkFailed = 1; }
	}

    ## ANNOTATE_VARIANTS
    if($opt{ANNOTATE_VARIANTS} && $opt{ANNOTATE_VARIANTS} eq "yes"){
        if(! $opt{SNPEFF_PATH}){ print "ERROR: No SNPEFF_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{IGVTOOLS_PATH}){ print "ERROR: No IGVTOOLS_PATH option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{ANNOTATE_QUEUE}){ print "ERROR: No ANNOTATE_QUEUE option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{ANNOTATE_THREADS}){ print "ERROR: No ANNOTATE_THREADS option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{ANNOTATE_MEM}){ print "ERROR: No ANNOTATE_MEM option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{ANNOTATE_TIME}){ print "ERROR: No ANNOTATE_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{ANNOTATE_SNPEFF}){ print "ERROR: No ANNOTATE_SNPEFF option found in config files.\n"; $checkFailed = 1; }
        if($opt{ANNOTATE_SNPEFF} eq "yes"){
            if(! $opt{ANNOTATE_DB}){ print "ERROR: No ANNOTATE_DB option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_FLAGS}){ print "ERROR: No ANNOTATE_FLAGS option found in config files.\n"; $checkFailed = 1; }
        }
        if(! $opt{ANNOTATE_SNPSIFT}){ print "ERROR: No ANNOTATE_SNPSIFT option found in config files.\n"; $checkFailed = 1; }
        if($opt{ANNOTATE_SNPSIFT} eq "yes"){
            if(! $opt{ANNOTATE_DBNSFP}){ print "ERROR: No ANNOTATE_DBNSFP option found in config files.\n"; $checkFailed = 1; }
            elsif( $opt{ANNOTATE_DBNSFP} && ! -e $opt{ANNOTATE_DBNSFP}) { print"ERROR: $opt{ANNOTATE_DBNSFP} does not exist\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_FIELDS}){ print "ERROR: No ANNOTATE_FIELDS option found in config files.\n"; $checkFailed = 1; }
        }
        if(! $opt{ANNOTATE_FREQUENCIES}){ print "ERROR: No ANNOTATE_FREQUENCIES option found in config files.\n"; $checkFailed = 1; }
        if($opt{ANNOTATE_FREQUENCIES} eq "yes"){
            if(! $opt{ANNOTATE_FREQNAME}){ print "ERROR: No ANNOTATE_FREQNAME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_FREQDB}){ print "ERROR: No ANNOTATE_FREQDB option found in config files.\n"; $checkFailed = 1; }
            elsif( $opt{ANNOTATE_FREQDB} && ! -e $opt{ANNOTATE_FREQDB}) { print"ERROR: $opt{ANNOTATE_FREQDB} does not exist\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_FREQINFO}){ print "ERROR: No ANNOTATE_FREQINFO option found in config files.\n"; $checkFailed = 1; }
        }
        if(! $opt{ANNOTATE_IDFIELD}){ print "ERROR: No ANNOTATE_IDFIELD option found in config files.\n"; $checkFailed = 1; }
            if($opt{ANNOTATE_IDFIELD} eq "yes"){
            if(! $opt{ANNOTATE_IDNAME}){ print "ERROR: No ANNOTATE_IDNAME option found in config files.\n"; $checkFailed = 1; }
            if(! $opt{ANNOTATE_IDDB}){ print "ERROR: No ANNOTATE_IDDB option found in config files.\n"; $checkFailed = 1; }
        }
    }

    ## KINSHIP
    if($opt{KINSHIP} && $opt{KINSHIP} eq "yes") {
        if(! $opt{KINSHIP_QUEUE}){ print "ERROR: No KINSHIP_QUEUE found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{KINSHIP_THREADS}){ print "ERROR: No KINSHIP_THREADS found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{KINSHIP_MEM}){ print "ERROR: No KINSHIP_MEM found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{KINSHIP_TIME}){ print "ERROR: No KINSHIP_TIME option found in config files.\n"; $checkFailed = 1; }
        if(! $opt{PLINK_PATH}){ print "ERROR: No PLINK_PATH found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{KING_PATH}){ print "ERROR: No KING_PATH found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{VCFTOOLS_PATH}){ print "ERROR: No VCFTOOLS_PATH found in .ini file\n"; $checkFailed = 1; }
    }

    ## FINALIZE
    if($opt{FINALIZE} && $opt{FINALIZE} eq "yes") {
        if(! $opt{FINALIZE_QUEUE}){ print "ERROR: No FINALIZE_QUEUE found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{FINALIZE_THREADS}){ print "ERROR: No FINALIZE_THREADS found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{FINALIZE_MEM}){ print "ERROR: No FINALIZE_MEM found in .ini file\n"; $checkFailed = 1; }
        if(! $opt{FINALIZE_TIME}){ print "ERROR: No FINALIZE_TIME found in .ini file\n"; $checkFailed = 1; }
    }

    if ($checkFailed) {
        die "One or more options not found in config files";
    }
}

1;
