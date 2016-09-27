package illumina_mapping;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;
use Carp;
use List::MoreUtils qw(zip);

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge;
use illumina_template;

sub validateFastQName {
    my ($input_file) = @_;
    my $name = fileparse($input_file);
    (my $R1 = $input_file) =~ s/_R2/_R1/;
    (my $R2 = $input_file) =~ s/_R1/_R2/;

    my $fastQPattern = qr/^(?<sampleName>[^_]+)_(?<flowcellID>[^_]+)_(?<index>[^_]+)_(?<lane>[^_]+)_(?<tag>R1|R2)_(?<suffix>[^\.]+)(?<ext>\.fastq\.gz)$/;
    $name =~ $fastQPattern or confess "ERROR: FASTQ filename '$name' must match regex '$fastQPattern' (for example: SAMPLENAME_FLOWCELLID_S1_L001_R1_001.fastq.gz)\n";

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
    my $runName = basename($opt{OUTPUT_DIR});

    my $FAI = "$opt{GENOME}.fai";
    die "GENOME: $opt{GENOME} does not exist!" if !-e $opt{GENOME};
    die "GENOME BWT: $opt{GENOME}.bwt does not exist!" if !-e "$opt{GENOME}.bwt";
    die "GENOME FAI: $FAI does not exist!" if !-e $FAI;

    my $mainJobID = catfile($opt{OUTPUT_DIR}, "jobs", "MapMainJob_" . getJobId() . ".sh");
    open (my $QSUB, ">", $mainJobID) or die "ERROR: Couldn't create $mainJobID: $!";
    print $QSUB "\#!/usr/bin/env bash\n\nsource $opt{CLUSTER_PATH}/settings.sh\n\n";

    my $samples = {};
    foreach my $input_fastq (keys %{$opt{FASTQ}}) {
        my $metadata = validateFastQName($input_fastq);
        say "Skipping R2 sample $input_fastq" and next if $input_fastq eq $metadata->{R2};
        if (exists $opt{FASTQ}->{$metadata->{R2}}) {
            say "Switching to paired end mode!";
        } else {
            say "Switching to fragment mode!";
            $opt{SINGLE_END} = 1;
        }
        say "Creating $opt{OUTPUT_DIR}/$metadata->{sampleName}/mapping/$metadata->{coreName}_sorted.bam with:";
        createIndividualMappingJobs(\%opt, $QSUB, $samples, $metadata->{sampleName}, $metadata->{coreName}, $metadata->{R1}, $metadata->{R2}, $metadata->{flowcellID});
    }

    say "";
    foreach my $sample (keys %{$samples}) {
        my @bamList = ();
        my @jobIds = ();

        foreach my $chunk (@{$samples->{$sample}}) {
            push @bamList, $chunk->{'file'};
            push @jobIds, $chunk->{'jobId'};
        }

        $opt{BAM_FILES}->{$sample} = "${sample}_dedup.bam";
        say "Creating $opt{BAM_FILES}->{$sample}";

        my $done_file = catfile($opt{OUTPUT_DIR}, $sample, "logs", "Mapping_${sample}.done");
        if (-e $done_file) {
            say "WARNING: $done_file exists, skipping";
            next;
        }

        my $jobId = "Merge_${sample}_" . getJobId();
        my $bash_file = catfile($opt{OUTPUT_DIR}, $sample, "jobs", "${jobId}.sh");
        push @{$opt{RUNNING_JOBS}->{$sample}}, $jobId;
        my $bams = join " ", @bamList;

        from_template("Merge.sh.tt", $bash_file, sample => $sample, bams => $bams, runName => $runName, opt => \%opt);

        my $qsub = qsubTemplate(\%opt, "MARKDUP");
        my $stdout = catfile($opt{OUTPUT_DIR}, $sample, "logs", "Merge_${sample}.out");
        my $stderr = catfile($opt{OUTPUT_DIR}, $sample, "logs", "Merge_${sample}.err");
        my $hold_jids = join ",", @jobIds;
        print $QSUB "$qsub -o $stdout -e $stderr -N $jobId -hold_jid $hold_jids $bash_file\n\n";
    }

    close $QSUB;

    system("sh $mainJobID");
    return \%opt;
}

sub createIndividualMappingJobs {
    my ($opt, $QSUB, $samples, $sampleName, $coreName, $R1, $R2, $flowcellID) = @_;
    my %opt = %$opt;
    my $runName = basename($opt{OUTPUT_DIR});
    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM, $RG_PU) = ('ILLUMINA', $coreName, $sampleName, $sampleName, $flowcellID);

    my $mappingJobId = "Map_${coreName}_" . getJobId();
    my $mappingFSJobId = "MapFS_${coreName}_" . getJobId();
    my $sortJobId = "Sort_${coreName}_" . getJobId();
    my $sortFSJobId = "SortFS_${coreName}_" . getJobId();
    my $indexJobId = "Index_${coreName}_" . getJobId();
    my $cleanAndCheckJobId = "CheckAndClean_${coreName}_" . getJobId();

    push @{$samples->{$sampleName}}, {
                                      'jobId' => $cleanAndCheckJobId,
                                      'file' => catfile($opt{OUTPUT_DIR}, $sampleName, "mapping", "${coreName}_sorted.bam")
                                     };

    my $done_file = catfile($opt{OUTPUT_DIR}, $sampleName, "mapping", "${coreName}.done");
    if (-e $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    $done_file = catfile($opt{OUTPUT_DIR}, $sampleName, "logs", "${coreName}_bwa.done");
    if (!-e $done_file) {
        print $R2 ? "\t$R1\n\t$R2\n" : "\t$R1\n";
        from_template("PerLaneMap.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh", coreName => $coreName, sampleName => $sampleName, R1 => $R1, R2 => $R2,
            RG_ID => $RG_ID, RG_SM => $RG_SM, RG_PL => $RG_PL, RG_LB => $RG_LB, RG_PU => $RG_PU, runName => $runName, opt => \%opt);

        my $qsub = qsubTemplate(\%opt, "MAPPING");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
        $mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$mappingJobId,".sh\n";
    } else {
        say "WARNING: $done_file exists, skipping BWA";
    }

    if (!-s "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat") {
        from_template("PerLaneMapFS.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh", sampleName => $sampleName, coreName => $coreName, runName => $runName, opt => \%opt);

        my $qsub = qsubTemplate(\%opt, "FLAGSTAT");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",
        $coreName,".err -N ",$mappingFSJobId," -hold_jid ",$mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$mappingFSJobId,".sh\n";
    } else {
        say "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat exist and is not empty, skipping bwa flagstat";
    }

    if (!-s "$opt{OUTPUT_DIR}/$sampleName/mapping/${coreName}_sorted.bam") {
        from_template("PerLaneSort.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh", coreName => $coreName, sampleName => $sampleName, runName => $runName, opt => \%opt);

        my $qsub = qsubTemplate(\%opt, "MAPPING");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",$sortJobId,
        " -hold_jid ",$mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$sortJobId,".sh\n";
    } else {
        say "\t$opt{OUTPUT_DIR}/$sampleName/mapping/${coreName}_sorted.bam exist and is not empty, skipping sort";
    }

    if (!-s "$opt{OUTPUT_DIR}/$sampleName/mapping/${coreName}_sorted.flagstat") {
        from_template("PerLaneSortFS.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh", sampleName => $sampleName, coreName => $coreName, runName => $runName, opt => \%opt);

        my $qsub = qsubTemplate(\%opt, "FLAGSTAT");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
        $sortFSJobId," -hold_jid ",$sortJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$sortFSJobId,".sh\n";
    } else {
        say "\t$opt{OUTPUT_DIR}/$sampleName/mapping/${coreName}_sortedlagstat exist and is not empty, skipping sorted bam flagstat";
    }

    if (!-s "$opt{OUTPUT_DIR}/$sampleName/mapping/${coreName}_sorted.bai") {
        from_template("PerLaneIndex.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh", sampleName => $sampleName, coreName => $coreName, runName => $runName, opt => \%opt);
        my $qsub = qsubTemplate(\%opt, "MAPPING");
        print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
        $indexJobId," -hold_jid ",$sortJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$indexJobId,".sh\n";
    } else {
        say "\t$opt{OUTPUT_DIR}/$sampleName/mapping/${coreName}_sorted.bai exist and is not empty, skipping sorted bam index";
    }

    my $checkAndCleanSh = "$opt{OUTPUT_DIR}/$sampleName/jobs/$cleanAndCheckJobId.sh";
    from_template("PerLaneCheckAndClean.sh.tt", $checkAndCleanSh, sampleName => $sampleName, coreName => $coreName, opt => \%opt);

    my $qsub = qsubTemplate(\%opt, "MAPPING");
    print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",$cleanAndCheckJobId,
    " -hold_jid ",$mappingFSJobId,",",$sortFSJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$cleanAndCheckJobId,".sh\n\n";
}

sub runBamPrep {
    my ($opt) = @_;
    my %opt = %{$opt};

    while (my ($sample, $input_bam) = each $opt{SAMPLES}) {
        (my $input_bai = $input_bam) =~ s/\.bam$/\.bai/;
        (my $input_flagstat = $input_bam) =~ s/\.bam$/\.flagstat/;

        my $bam_file = "${sample}.bam";
        $opt{BAM_FILES}->{$sample} = $bam_file;
        my $sample_bam = catfile($opt{OUTPUT_DIR}, $sample, "mapping", $bam_file);
        my $sample_bai = catfile($opt{OUTPUT_DIR}, $sample, "mapping", "${sample}.bai");
        my $sample_flagstat = catfile($opt{OUTPUT_DIR}, $sample, "mapping", "${sample}.flagstat");

        my $bai_good = verifyBai($input_bai, $input_bam, $opt);
        my $flagstat_good = verifyFlagstat($input_flagstat, $input_bam);

        symlink($input_bam, $sample_bam);
        $bai_good and symlink($input_bai, $sample_bai);
        $flagstat_good and symlink($input_flagstat, $sample_flagstat);

        next if $bai_good and $flagstat_good;

        my $job_id = "PrepBam_${sample}_" . getJobId();
        my $log_dir = catfile($opt{OUTPUT_DIR}, $sample, "logs");
        my $bash_file = catfile($opt{OUTPUT_DIR}, $sample, "jobs", "$job_id.sh");

        from_template("PrepBam.sh.tt", $bash_file,
                      sample => $sample,
                      sample_bam => $sample_bam,
                      sample_bai => $sample_bai,
                      sample_flagstat => $sample_flagstat,
                      log_dir => $log_dir,
                      opt => $opt);

        my $qsub = qsubTemplate($opt, "MAPPING");
        system "$qsub -o $log_dir/PrepBam_${sample}.out -e $log_dir/PrepBam_${sample}.err -N $job_id $bash_file";
        push @{$opt{RUNNING_JOBS}->{$sample}}, $job_id;
    }
    return $opt;
}

sub verifyBam {
    my ($bam_file, $opt) = @_;

    -e $bam_file or confess "ERROR: $bam_file does not exist";
    (my $sample = fileparse($bam_file)) =~ s/\.bam$//;

    my $headers = bamHeaders($bam_file, $opt);
    my @read_groups = grep $_->{name} eq '@RG', @$headers;

    my @sample_names = map $_->{tags}{SM}, @read_groups;
    confess "too many samples in BAM $bam_file: @sample_names" unless keys { map { $_, 1 } @sample_names } < 2;
    warn "missing sample name in BAM $bam_file, using file name" unless @sample_names;
    $sample = $sample_names[0] if @sample_names;

    my %header_contigs = map { $_->{tags}{SN}, $_->{tags}{LN} } grep $_->{name} eq '@SQ', @$headers;
    verifyContigs(\%header_contigs, refGenomeContigs($opt));
    verifyReadGroups(\@read_groups, bamReads($bam_file, 1000, $opt));

    return $sample;
}

sub verifyBai {
    my ($bai_file, $bam_file, $opt) = @_;

    -e $bai_file && -M $bam_file > -M $bai_file or return 0;

    # this check does not happen if the .bai is missing/old, so re-implemented in the job :(
    my $headers = bamHeaders($bam_file, $opt);
    my %header_contigs = map { $_->{tags}{SN}, $_->{tags}{LN} } grep $_->{name} eq '@SQ', @$headers;
    verifyContigs(indexContigs($bam_file, $opt), \%header_contigs);
    return 1;
}

sub verifyFlagstat {
    my ($flagstat_file, $bam_file) = @_;

    return -e $flagstat_file && -M $bam_file > -M $flagstat_file;
}

sub bamHeaders {
    my ($bam_file, $opt) = @_;

    my $samtools = catfile($opt->{SAMTOOLS_PATH}, "samtools");

    my @lines = qx($samtools view -H $bam_file);
    $? == 0 or confess "could not read BAM headers from $bam_file";

    chomp @lines;
    my @fields = map [ split qr/[\t:]/ ], @lines;
    my @headers = map {
        name => shift $_,
        tags => { @$_ },
    }, @fields;
    return \@headers;
}

sub bamReads {
    my ($bam_file, $num_lines, $opt) = @_;

    my $samtools = catfile($opt->{SAMTOOLS_PATH}, "samtools");

    my @lines = qx($samtools view $bam_file | head -$num_lines);
    chomp @lines;

    my @field_names = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual tags);
    my @fields = map [ split "\t", $_, @field_names ], @lines;
    my @reads = map +{ zip @field_names, @$_ }, @fields;
    map { $_->{tags} =
          {
           map {
               shift @$_,
               {
                type => shift @$_,
                value => shift @$_,
               },
           } map [ split ":" ], split "\t", $_->{tags}
          },
      } @reads;
    return \@reads;
}

sub refGenomeContigs {
    my ($opt) = @_;

    my $fasta_file = "$opt->{GENOME}.fai";
    my @lines = qx(cat $fasta_file);
    $? == 0 or confess "could not read from $fasta_file";

    return contigLengths(\@lines);
}

sub indexContigs {
    my ($bam_file, $opt) = @_;

    my $samtools = catfile($opt->{SAMTOOLS_PATH}, "samtools");
    my @lines = qx($samtools idxstats $bam_file);
    $? == 0 or confess "could not read index stats from $bam_file";

    my $contigs = contigLengths(\@lines);
    delete $contigs->{'*'};
    return $contigs;
}

sub contigLengths {
    my ($lines) = @_;

    chomp @$lines;
    my %contigs = map { splice [ split "\t" ], 0, 2 } @$lines;
    return \%contigs;
}

sub verifyContigs {
    my ($contigs, $ref_contigs) = @_;

    my @warnings;
    while (my ($key, $value) = each $contigs) {
        push @warnings, "contig $key missing" and next if not exists $ref_contigs->{$key};
        push @warnings, "contig $key length $value does not match $ref_contigs->{$key}" if $value != $ref_contigs->{$key};
    }
    warn $_ foreach @warnings;
    confess "contigs do not match" if @warnings;
}

sub verifyReadGroups {
    my ($read_groups, $bam_reads) = @_;

    my %header_rgids = map { $_->{tags}{ID} => 1 } @$read_groups;
    my %reads_rgids = map { $_->{tags}{RG}{value} => 1 } @$bam_reads;
    my @unknown_rgids = grep { not exists $header_rgids{$_} } keys %reads_rgids;
    confess "read group IDs from read tags not in BAM header:\n\t" . join "\n\t", @unknown_rgids if @unknown_rgids;
}

1;
