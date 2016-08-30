#!/usr/bin/perl -w

package illumina_copyNumber;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);
use lib "$FindBin::Bin";
use illumina_sge;
use illumina_metadataParser;

sub runCopyNumberTools {
    my $configuration = shift;
    my %opt = %{$configuration};
    my @check_cnv_jobs;

    if ($opt{CNV_MODE} eq "sample_control") {
        my @cnv_jobs;
        my $metadata = metadataParse($opt{OUTPUT_DIR});
        my $sample_ref = $metadata->{'ref_sample'};
        my $sample_tumor = $metadata->{'tumor_sample'};

        my $sample_tumor_name = "$sample_ref\_$sample_tumor";
        my $sample_tumor_out_dir = "$opt{OUTPUT_DIR}/copyNumber/$sample_tumor_name";
        my $sample_tumor_log_dir = "$sample_tumor_out_dir/logs/";
        my $sample_tumor_job_dir = "$sample_tumor_out_dir/jobs/";

        if (!-e $sample_tumor_out_dir) {
            make_path($sample_tumor_out_dir) or die "Couldn't create directory:  $sample_tumor_out_dir\n";
        }
        if (!-e $sample_tumor_job_dir) {
            make_path($sample_tumor_job_dir) or die "Couldn't create directory: $sample_tumor_job_dir\n";
        }
        if (!-e $sample_tumor_log_dir) {
            make_path($sample_tumor_log_dir) or die "Couldn't create directory: $sample_tumor_log_dir\n";
        }

        my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
        my @running_jobs;
        if (@{$opt{RUNNING_JOBS}->{$sample_tumor}}) {
            push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
        }
        my $sample_ref_bam = "$opt{OUTPUT_DIR}/$sample_ref/mapping/$opt{BAM_FILES}->{$sample_ref}";
        if (@{$opt{RUNNING_JOBS}->{$sample_ref}}) {
            push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_ref}});
        }

        print "\n$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

        if (-e "$sample_tumor_log_dir/$sample_tumor_name.done") {
            print "WARNING: $sample_tumor_log_dir/$sample_tumor_name.done exists, skipping \n";
            next;
        }

        if ($opt{CNV_FREEC} eq "yes") {
            print "\n###SCHEDULING FREEC####\n";
            my $freec_job = runFreec($sample_tumor, $sample_tumor_out_dir, $sample_tumor_job_dir,
                $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
            if ($freec_job) {push(@cnv_jobs, $freec_job)};
        }

        my $job_id = "CHECK_".$sample_tumor."_".get_job_id();
        my $bash_file = $sample_tumor_job_dir."/".$job_id.".sh";

        open CHECK_SH, ">$bash_file" or die "cannot open file $bash_file \n";
        print CHECK_SH "#!/bin/bash\n\n";
        print CHECK_SH "echo \"Start Check\t\" `date` `uname -n` >> $sample_tumor_log_dir/check.log\n\n";
        print CHECK_SH "if [[ ";
        if ($opt{CNV_FREEC} eq "yes") {print CHECK_SH "-f $sample_tumor_log_dir/freec.done"}
        print CHECK_SH " ]]\n";
        print CHECK_SH "then\n";
        print CHECK_SH "\ttouch $sample_tumor_log_dir/$sample_tumor_name.done\n";
        print CHECK_SH "fi\n\n";
        print CHECK_SH "echo \"End Check\t\" `date` `uname -n` >> $sample_tumor_log_dir/check.log\n";
        close CHECK_SH;
        my $qsub = &qsubTemplate(\%opt, "CNVCHECK");
        if (@cnv_jobs) {
            system "$qsub -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id -hold_jid ".join(","
                    , @cnv_jobs)." $bash_file";
        } else {
            system "$qsub -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id $bash_file";
        }
        push(@check_cnv_jobs, $job_id);
    } elsif ($opt{CNV_MODE} eq "sample") {
        foreach my $sample (@{$opt{SAMPLES}}) {
            my @cnv_jobs;

            my $sample_out_dir = "$opt{OUTPUT_DIR}/copyNumber/$sample";
            my $sample_log_dir = "$sample_out_dir/logs/";
            my $sample_job_dir = "$sample_out_dir/jobs/";

            if (!-e $sample_out_dir) {
                make_path($sample_out_dir) or die "Couldn't create directory:  $sample_out_dir\n";
            }
            if (!-e $sample_job_dir) {
                make_path($sample_job_dir) or die "Couldn't create directory: $sample_job_dir\n";
            }
            if (!-e $sample_log_dir) {
                make_path($sample_log_dir) or die "Couldn't create directory: $sample_log_dir\n";
            }

            my $sample_bam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
            my @running_jobs;
            if (@{$opt{RUNNING_JOBS}->{$sample}}) {
                push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample}});
            }
            print "\n$sample \t $sample_bam \n";

            if (-e "$sample_log_dir/$sample.done") {
                print "WARNING: $sample_log_dir/$sample.done exists, skipping \n";
                next;
            }

            if ($opt{CNV_FREEC} eq "yes") {
                print "\n###SCHEDULING FREEC####\n";
                my $freec_job = runFreec($sample, $sample_out_dir, $sample_job_dir, $sample_log_dir, $sample_bam, "",
                    \@running_jobs, \%opt);
                if ($freec_job) {push(@cnv_jobs, $freec_job)};
            }

            my $job_id = "CHECK_".$sample."_".get_job_id();
            my $bash_file = $sample_job_dir."/".$job_id.".sh";

            open CHECK_SH, ">$bash_file" or die "cannot open file $bash_file \n";
            print CHECK_SH "#!/bin/bash\n\n";
            print CHECK_SH "echo \"Start Check\t\" `date` `uname -n` >> $sample_log_dir/check.log\n\n";
            print CHECK_SH "if [[ -f $sample_log_dir/freec.done ]]\n";
            print CHECK_SH "then\n";
            print CHECK_SH "\ttouch $sample_log_dir/$sample.done\n";
            print CHECK_SH "fi\n\n";
            print CHECK_SH "echo \"End Check\t\" `date` `uname -n` >> $sample_log_dir/check.log\n";
            close CHECK_SH;
            my $qsub = &qsubTemplate(\%opt, "CNVCHECK");
            if (@cnv_jobs) {
                system "$qsub -o $sample_log_dir -e $sample_log_dir -N $job_id -hold_jid ".join(",",
                        @cnv_jobs)." $bash_file";
            } else {
                system "$qsub -o $sample_log_dir -e $sample_log_dir -N $job_id $bash_file";
            }
            push(@check_cnv_jobs, $job_id);
        }
    }
    return \@check_cnv_jobs;
}

sub runFreec {
    my ($sample_name, $out_dir, $job_dir, $log_dir, $sample_bam, $control_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};

    if (-e "$log_dir/freec.done") {
        print "WARNING: $log_dir/freec.done exists, skipping \n";
        return;
    }

    my $freec_out_dir = "$out_dir/freec";
    if (!-e $freec_out_dir) {
        make_path($freec_out_dir) or die "Couldn't create directory: $freec_out_dir\n";
    }

    my @mappabilityTracks;
    if ($opt{FREEC_MAPPABILITY_TRACKS}) {
        @mappabilityTracks = split('\t', $opt{FREEC_MAPPABILITY_TRACKS});
    }

    my $freec_config = $freec_out_dir."/freec_config.txt";
    open FREEC_CONFIG, ">$freec_config" or die "cannot open file $freec_config \n";

    print FREEC_CONFIG "[general]\n";
    print FREEC_CONFIG "chrLenFile= $opt{FREEC_CHRLENFILE}\n";
    print FREEC_CONFIG "ploidy=2\n";
    print FREEC_CONFIG "samtools=$opt{SAMTOOLS_PATH}/samtools\n";
    print FREEC_CONFIG "chrFiles= $opt{FREEC_CHRFILES}\n";
    print FREEC_CONFIG "window=$opt{FREEC_WINDOW}\n";
    print FREEC_CONFIG "maxThreads=$opt{FREEC_THREADS}\n";
    print FREEC_CONFIG "telocentromeric=$opt{FREEC_TELOCENTROMERIC}\n";
    print FREEC_CONFIG "BedGraphOutput=TRUE\n";
    print FREEC_CONFIG "outputDir=$freec_out_dir\n";

    foreach my $mappabilityTrack (@mappabilityTracks) {
        print FREEC_CONFIG "gemMappabilityFile=$mappabilityTrack\n";
    }

    print FREEC_CONFIG "[sample]\n";
    print FREEC_CONFIG "mateFile=$sample_bam\n";
    print FREEC_CONFIG "inputFormat=BAM\n";
    print FREEC_CONFIG "mateOrientation=FR\n";
    if ($control_bam) {
        print FREEC_CONFIG "[control]\n";
        print FREEC_CONFIG "mateFile=$control_bam\n";
        print FREEC_CONFIG "inputFormat=BAM\n";
        print FREEC_CONFIG "mateOrientation=FR\n";
    }

    if ($opt{CNV_TARGETS}) {
        print FREEC_CONFIG "[target]\n";
        print FREEC_CONFIG "captureRegions=$opt{CNV_TARGETS}\n";
    }

    close FREEC_CONFIG;

    my $job_id = "FREEC_".$sample_name."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";
    my $sample_bam_name = (split('/', $sample_bam))[-1];

    open FREEC_SH, ">$bash_file" or die "cannot open file $bash_file \n";

    print FREEC_SH "#!/bin/bash\n\n";
    if ($control_bam) {
        print FREEC_SH "if [ -s $sample_bam -a -s $control_bam ]\n";
    } else {
        print FREEC_SH "if [ -s $sample_bam ]\n";
    }
    print FREEC_SH "then\n";
    print FREEC_SH "\techo \"Start FREEC\t\" `date` \"\t $sample_bam \t $control_bam\t\" `uname -n` >> $log_dir/freec.log\n\n";

    print FREEC_SH "\t$opt{FREEC_PATH}/freec -conf $freec_config\n";
    print FREEC_SH "\tcd $freec_out_dir\n";
    print FREEC_SH "\tcat $opt{FREEC_PATH}/assess_significance.R | R --slave --args ".$sample_bam_name."_CNVs ".$sample_bam_name."_ratio.txt\n";
    print FREEC_SH "\tcat $opt{FREEC_PATH}/makeGraph.R | R --slave --args 2 ".$sample_bam_name."_ratio.txt\n";
    print FREEC_SH "\tcat $opt{PIPELINE_PATH}/scripts/makeKaryotype.R | R --slave --args 2 24 4 500000 ".$sample_bam_name."_ratio.txt\n";
    print FREEC_SH "\ttouch $log_dir/freec.done\n";
    print FREEC_SH "\techo \"End FREEC\t\" `date` \"\t $sample_bam \t $control_bam\t\" `uname -n` >> $log_dir/freec.log\n\n";
    print FREEC_SH "else\n";
    print FREEC_SH "\techo \"ERROR: $sample_bam or $control_bam does not exist.\" >> $log_dir/freec.log\n";
    print FREEC_SH "fi\n";

    close FREEC_SH;
    my $qsub = &qsubTemplate(\%opt, "FREEC");

    if (@running_jobs) {
        system "$qsub -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",", @running_jobs)." $bash_file";
    } else {
        system "$qsub -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id =~ s/\/tmp\/file//;
    return $id;
}
############

1;
