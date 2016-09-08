package illumina_mapping;

use strict;
use warnings;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin";
use File::Basename;
use File::Spec::Functions;

use illumina_sge;
use illumina_template;

sub validateFastQName {
    my ($input) = @_;
    my $name = fileparse($input);
    (my $R1 = $input) =~ s/_R2/_R1/;
    (my $R2 = $input) =~ s/_R1/_R2/;

    my $fastQPattern = qr/^(?<sampleName>[^_]+)_(?<flowcellID>[^_]+)_(?<index>[^_]+)_(?<lane>[^_]+)_(?<tag>R1|R2)_(?<suffix>[^\.]+)(?<ext>\.fastq\.gz)$/;
    $name =~ $fastQPattern or die "ERROR: FASTQ filename '$name' must match regex '$fastQPattern' (for example: SAMPLENAME_FLOWCELLID_S1_L001_R1_001.fastq.gz)\n";

    return {
        R1 => $R1,
        R2 => $R2,
        coreName => "$+{sampleName}_$+{flowcellID}_$+{index}_$+{lane}_$+{suffix}",
        %+,
    };
}

sub runMapping {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];

    my $FAI = "$opt{GENOME}\.fai";
    die "GENOME: $opt{GENOME} does not exists!\n" if !-e "$opt{GENOME}";
    die "GENOME BWT: $opt{GENOME}.bwt does not exists!\n" if !-e "$opt{GENOME}.bwt";
    die "GENOME FAI: $FAI does not exists!\n" if !-e $FAI;

    my $mainJobID = "$opt{OUTPUT_DIR}/jobs/MapMainJob_".get_job_id().".sh";

    open (my $QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print $QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";

    my $samples = {};

    foreach my $input (keys %{$opt{FASTQ}}) {
        my $metadata = validateFastQName($input);
        print "Skipping R2 sample $input\n" and next if $input eq $metadata->{R2};
        if (exists $opt{FASTQ}->{$metadata->{R2}}) {
            print "Switching to paired end mode!\n";
        } else {
            print "Switching to fragment mode!\n";
            $opt{SINGLE_END} = 1;
        }
        print "Creating $opt{OUTPUT_DIR}/$metadata->{sampleName}/mapping/$metadata->{coreName}\_sorted.bam with:\n";
        createIndividualMappingJobs(\%opt, $QSUB, $samples, $metadata->{sampleName}, $metadata->{coreName}, $metadata->{R1}, $metadata->{R2}, $metadata->{flowcellID});
    }

    print "\n";
    foreach my $sample (keys %{$samples}) {
        my @bamList = ();
        my @jobIds = ();

        foreach my $chunk (@{$samples->{$sample}}) {
            push( @bamList, $chunk->{'file'} );
            push( @jobIds, $chunk->{'jobId'} );
        }

        $opt{BAM_FILES}->{$sample} = "$sample\_dedup.bam";
        print "Creating $opt{BAM_FILES}->{$sample}\n";

        if (-e "$opt{OUTPUT_DIR}/$sample/logs/Mapping_$sample.done") {
            print "\tWARNING: $opt{OUTPUT_DIR}/$sample/logs/Mapping_$sample.done exists, skipping\n";
            next;
        }

        my $jobId = "Merge_$sample\_".get_job_id();
        push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobId);
        my $bams = join(" ", @bamList);

        from_template("Merge.sh.tt", "$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh", sample => $sample, bamList => \@bamList, bams => $bams, runName => $runName, opt => \%opt);

        my $qsub = &qsubTemplate(\%opt, "MARKDUP");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sample,"/logs/Merge_",$sample,".out -e ",$opt{OUTPUT_DIR},"/",$sample,"/logs/Merge_",$sample,".err -N ",$jobId,
        " -hold_jid ",join(",",@jobIds)," ",$opt{OUTPUT_DIR},"/",$sample,"/jobs/",$jobId,".sh\n\n";
    }

    close $QSUB;

    system("sh $mainJobID");
    return \%opt;
}

sub createIndividualMappingJobs {
    my ($opt, $QSUB, $samples, $sampleName, $coreName, $R1, $R2, $flowcellID) = @_;
    my %opt = %$opt;
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM, $RG_PU) = ('ILLUMINA', $coreName, $sampleName, $sampleName, $flowcellID);

    my $mappingJobId = "Map_$coreName\_".get_job_id();
    my $mappingFSJobId = "MapFS_$coreName\_".get_job_id();
    my $sortJobId = "Sort_$coreName\_".get_job_id();
    my $sortFSJobId = "SortFS_$coreName\_".get_job_id();
    my $indexJobId = "Index_$coreName\_".get_job_id();
    my $cleanAndCheckJobId = "CheckAndClean_$coreName\_".get_job_id();

    push(@{$samples->{$sampleName}}, {'jobId'=>$cleanAndCheckJobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam"});

    if (-e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done") {
        print "\tWARNING: $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done exists, skipping\n";
        return;
    }

    if (! -e "$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa.done"){
        print $R2 ? "\t$R1\n\t$R2\n" : "\t$R1\n";
        from_template("PerLaneMap.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh", coreName => $coreName, sampleName => $sampleName, R1 => $R1, R2 => $R2,
            RG_ID => $RG_ID, RG_SM => $RG_SM, RG_PL => $RG_PL, RG_LB => $RG_LB, RG_PU => $RG_PU, runName => $runName, opt => \%opt);

        my $qsub = &qsubTemplate(\%opt,"MAPPING");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
        $mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$mappingJobId,".sh\n";
    } else {
        print "\t$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa.done exists, skipping bwa\n";
    }

    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat")){
        from_template("PerLaneMapFS.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh", sampleName => $sampleName, coreName => $coreName, runName => $runName, opt => \%opt);

        my $qsub = &qsubTemplate(\%opt,"FLAGSTAT");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",
        $coreName,".err -N ",$mappingFSJobId," -hold_jid ",$mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$mappingFSJobId,".sh\n";
    } else {
        print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat exist and is not empty, skipping bwa flagstat\n";
    }

    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam")) {
        from_template("PerLaneSort.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh", coreName => $coreName, sampleName => $sampleName, runName => $runName, opt => \%opt);

        my $qsub = &qsubTemplate(\%opt,"MAPPING");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",$sortJobId,
        " -hold_jid ",$mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$sortJobId,".sh\n";
    } else {
        print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam exist and is not empty, skipping sort\n";
    }

    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat")){
        from_template("PerLaneSortFS.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh", sampleName => $sampleName, coreName => $coreName, runName => $runName, opt => \%opt);

        my $qsub = &qsubTemplate(\%opt,"FLAGSTAT");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
        $sortFSJobId," -hold_jid ",$sortJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$sortFSJobId,".sh\n";
    } else {
        print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat exist and is not empty, skipping sorted bam flagstat\n";
    }

    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai")) {
        from_template("PerLaneIndex.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh", sampleName => $sampleName, coreName => $coreName, runName => $runName, opt => \%opt);
        my $qsub = &qsubTemplate(\%opt,"MAPPING");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
        $indexJobId," -hold_jid ",$sortJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$indexJobId,".sh\n";
    } else {
        print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai exist and is not empty, skipping sorted bam index\n";
    }

    my $checkAndCleanSh = "$opt{OUTPUT_DIR}/$sampleName/jobs/$cleanAndCheckJobId.sh";
    from_template("PerLaneCheckAndClean.sh.tt", $checkAndCleanSh, sampleName => $sampleName, coreName => $coreName, opt => \%opt);

    my $qsub = &qsubTemplate(\%opt,"MAPPING");
    print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",$cleanAndCheckJobId,
    " -hold_jid ",$mappingFSJobId,",",$sortFSJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$cleanAndCheckJobId,".sh\n\n";
}

sub runBamPrep {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $jobIds = {};

    foreach my $input (keys %{$opt{BAM}}) {
        my $bam_file = fileparse($input);
        (my $input_bai = $input) =~ s/\.bam/\.bai/g;
        (my $input_flagstat = $input) =~ s/\.bam/\.flagstat/g;

        (my $sample = $bam_file) =~ s/\.bam//g;
        $opt{BAM_FILES}->{$sample} = $bam_file;
        my $sample_bam = catfile($opt{OUTPUT_DIR}, $sample, "mapping", $bam_file);
        my $sample_bai = catfile($opt{OUTPUT_DIR}, $sample, "mapping", "${sample}.bai");
        my $sample_flagstat = catfile($opt{OUTPUT_DIR}, $sample, "mapping", "${sample}.flagstat");

        -e $input or die "ERROR: $input does not exist.";
        symlink($input, $sample_bam);
        -e $input_bai and symlink($input_bai, $sample_bai);
        -e $input_flagstat and symlink($input_flagstat, $sample_flagstat);

        next if -e $sample_bai && -e $sample_flagstat;

        my $job_id = "PrepBam_${sample}_" . get_job_id();
        my $log_dir = catfile($opt{OUTPUT_DIR}, $sample, "logs");
        my $bash_file = catfile($opt{OUTPUT_DIR}, $sample, "jobs", "$job_id.sh");

        from_template("PrepBam.sh.tt", $bash_file,
                      sample => $sample,
                      sample_bam => $sample_bam,
                      sample_bai => $sample_bai,
                      sample_flagstat => $sample_flagstat,
                      log_dir => $log_dir,
                      opt => \%opt);

        my $qsub = &qsubTemplate(\%opt, "MAPPING");
        system "$qsub -o $log_dir/PrepBam_${sample}.out -e $log_dir/PrepBam_${sample}.err -N $job_id $bash_file";
        push(@{$opt{RUNNING_JOBS}->{$sample}}, $job_id);
    }
    return \%opt;
}

############
sub get_job_id {
   my $id = tmpnam();
      $id=~s/\/tmp\/file//;
   return $id;
}
############

1;
