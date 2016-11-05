package illumina_mapping;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;
use Carp;
use List::MoreUtils qw(zip);

use illumina_sge qw(qsubTemplate);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);

use parent qw(Exporter);
our @EXPORT_OK = qw(runMapping runBamPrep verifyBam);


sub runMapping {
    my ($opt) = @_;

    say "\n### SCHEDULING MAPPING ###";

    die "GENOME: $opt->{GENOME} does not exist!" if !-f $opt->{GENOME};
    die "GENOME BWT: $opt->{GENOME}.bwt does not exist!" if !-f "$opt->{GENOME}.bwt";
    die "GENOME FAI: $opt->{GENOME}.fai does not exist!" if !-f "$opt->{GENOME}.fai";

    my $samples = {};
    foreach my $input_fastq (keys %{$opt->{FASTQ}}) {
        my $metadata = validateFastQName($input_fastq);
        say "Skipping R2 sample $input_fastq" and next if $input_fastq eq $metadata->{R2};
        if (exists $opt->{FASTQ}->{$metadata->{R2}}) {
            say "Switching to paired end mode!";
        } else {
            say "Switching to fragment mode!";
            $opt->{SINGLE_END} = 1;
        }

        addDirectories($samples, $opt->{OUTPUT_DIR}, $metadata->{sampleName});
        say "Creating $samples->{$metadata->{sampleName}}{dirs}{mapping}/$metadata->{coreName}_sorted.bam with:";
        createIndividualMappingJobs($opt, $metadata, $samples);
    }

    say "";
    while (my ($sample, $info) = each %{$samples}) {
        $opt->{BAM_FILES}->{$sample} = "${sample}_dedup.bam";
        say "Creating $opt->{BAM_FILES}->{$sample}";

        my $done_file = catfile($info->{dirs}->{log}, "Mapping_${sample}.done");
        if (-f $done_file) {
            say "WARNING: $done_file exists, skipping";
            next;
        }

        my $bams = join " ", map { $_->{file} } @{$info->{jobs}};
        my $hold_jids = join ",", map { $_->{job_id} } @{$info->{jobs}};
        my $job_id = "MergeMarkdup_${sample}_" . getJobId();
        my $bash_file = catfile($info->{dirs}->{job}, "${job_id}.sh");

        from_template("MergeMarkdup.sh.tt", $bash_file, sample => $sample, bams => $bams, opt => $opt);
        my $qsub = qsubTemplate($opt, "MARKDUP");
        my $stdout = catfile($info->{dirs}->{log}, "MergeMarkdup_${sample}.out");
        my $stderr = catfile($info->{dirs}->{log}, "MergeMarkdup_${sample}.err");
        system("$qsub -o $stdout -e $stderr -N $job_id -hold_jid $hold_jids $bash_file");
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }
    return;
}

sub addDirectories {
    my ($samples, $run_dir, $sample_name) = @_;

    my $out_dir = catfile($run_dir, $sample_name);
    $samples->{$sample_name}->{dirs} = {
        out => $out_dir,
        log => catfile($out_dir, "logs"),
        job => catfile($out_dir, "jobs"),
        mapping => catfile($out_dir, "mapping"),
    };
    return;
}

sub createIndividualMappingJobs {
    my ($opt, $metadata, $samples) = @_;

    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM, $RG_PU) = ('ILLUMINA', $metadata->{coreName}, $metadata->{sampleName}, $metadata->{sampleName}, $metadata->{flowcellID});

    my %jid = (
        mapping => "Map_$metadata->{coreName}_" . getJobId(),
        mapping_flagstat => "MapFS_$metadata->{coreName}_" . getJobId(),
        sort => "Sort_$metadata->{coreName}_" . getJobId(),
        sort_flagstat => "SortFS_$metadata->{coreName}_" . getJobId(),
        check_clean => "CheckAndClean_$metadata->{coreName}_" . getJobId(),
    );

    my $dirs = $samples->{$metadata->{sampleName}}->{dirs};
    push @{$samples->{$metadata->{sampleName}}->{jobs}}, {
        job_id => $jid{check_clean},
        file => catfile($dirs->{mapping}, "$metadata->{coreName}_sorted.bam"),
        };

    my $done_file = catfile($dirs->{mapping}, "$metadata->{coreName}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my $stdout = catfile($dirs->{log}, "Mapping_$metadata->{coreName}.out");
    my $stderr = catfile($dirs->{log}, "Mapping_$metadata->{coreName}.err");
    my $core_file = catfile($dirs->{mapping}, $metadata->{coreName});

    $done_file = catfile($dirs->{log}, "$metadata->{coreName}_bwa.done");
    if (!-f $done_file) {
        say $metadata->{R2} ? "\t$metadata->{R1}\n\t$metadata->{R2}" : "\t$metadata->{R1}";

        my $bash_file = catfile($dirs->{job}, "$jid{mapping}.sh");
        from_template(
            "PerLaneMap.sh.tt", $bash_file,
            coreName => $metadata->{coreName},
            sampleName => $metadata->{sampleName},
            R1 => $metadata->{R1},
            R2 => $metadata->{R2},
            RG_ID => $RG_ID,
            RG_SM => $RG_SM,
            RG_PL => $RG_PL,
            RG_LB => $RG_LB,
            RG_PU => $RG_PU,
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "MAPPING");
        system("$qsub -o $stdout -e $stderr -N $jid{mapping} $bash_file");
    } else {
        say "WARNING: $done_file exists, skipping BWA";
    }

    if (!-s "${core_file}.flagstat") {
        my $bash_file = catfile($dirs->{job}, "$jid{mapping_flagstat}.sh");
        from_template(
            "PerLaneMapFS.sh.tt", $bash_file,
            sampleName => $metadata->{sampleName},
            coreName => $metadata->{coreName},
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "FLAGSTAT");
        system("$qsub -o $stdout -e $stderr -N $jid{mapping_flagstat} -hold_jid $jid{mapping} $bash_file");
    } else {
        say "\t${core_file}.flagstat exists and is not empty, skipping BWA flagstat";
    }

    if (!-s "${core_file}_sorted.bam") {
        my $bash_file = catfile($dirs->{job}, "$jid{sort}.sh");
        from_template(
            "PerLaneSort.sh.tt", $bash_file,
            coreName => $metadata->{coreName},
            sampleName => $metadata->{sampleName},
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "MAPPING");
        system("$qsub -o $stdout -e $stderr -N $jid{sort} -hold_jid $jid{mapping} $bash_file");
    } else {
        say "\t${core_file}_sorted.bam exists and is not empty, skipping sort";
    }

    if (!-s "${core_file}_sorted.flagstat") {
        my $bash_file = catfile($dirs->{job}, "$jid{sort_flagstat}.sh");
        from_template(
            "PerLaneSortFS.sh.tt", $bash_file,
            sampleName => $metadata->{sampleName},
            coreName => $metadata->{coreName},
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "FLAGSTAT");
        system("$qsub -o $stdout -e $stderr -N $jid{sort_flagstat} -hold_jid $jid{sort} $bash_file");
    } else {
        say "\t${core_file}_sorted.flagstat exists and is not empty, skipping sorted BAM flagstat";
    }

    my $bash_file = catfile($dirs->{job}, "$jid{check_clean}.sh");
    from_template(
        "PerLaneCheckAndClean.sh.tt", $bash_file,
        sampleName => $metadata->{sampleName},
        coreName => $metadata->{coreName},
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "MAPPING");
    system("$qsub -o $stdout -e $stderr -N $jid{check_clean} -hold_jid $jid{mapping_flagstat},$jid{sort_flagstat} $bash_file");
    return;
}

sub validateFastQName {
    my ($input_file) = @_;
    my $name = fileparse($input_file);
    (my $R1 = $input_file) =~ s/_R2/_R1/;
    (my $R2 = $input_file) =~ s/_R1/_R2/;

    my $fastQPattern = qr/^(?<sampleName>[^_]+)
                          _(?<flowcellID>[^_]+)
                          _(?<index>[^_]+)
                          _(?<lane>[^_]+)
                          _(?<tag>R1|R2)
                          _(?<suffix>[^\.]+)
                          (?<ext>\.fastq\.gz)$
                         /x;
    $name =~ $fastQPattern or confess "ERROR: FASTQ filename '$name' must match regex '$fastQPattern' (for example: SAMPLENAME_FLOWCELLID_S1_L001_R1_001.fastq.gz)\n";

    return {
        R1 => $R1,
        R2 => $R2,
        coreName => "$+{sampleName}_$+{flowcellID}_$+{index}_$+{lane}_$+{suffix}",
        %+,
    };
}

sub runBamPrep {
    my ($opt) = @_;

    $opt->{MAPPING} = "no";
    $opt->{PRESTATS} = "no";
    say "\n### SCHEDULING BAM PREP ###";

    while (my ($sample, $input_bam) = each %{$opt->{SAMPLES}}) {
        (my $input_bai = $input_bam) =~ s/\.bam$/.bam.bai/;
        (my $input_flagstat = $input_bam) =~ s/\.bam$/.flagstat/;

        my $bam_file = "${sample}.bam";
        $opt->{BAM_FILES}->{$sample} = $bam_file;
        my $sample_bam = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $bam_file);
        my $sample_bai = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", "${sample}.bam.bai");
        my $sample_flagstat = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", "${sample}.flagstat");

        my $bai_good = verifyBai($input_bai, $input_bam, $opt);
        my $flagstat_good = verifyFlagstat($input_flagstat, $input_bam);

        symlink($input_bam, $sample_bam);
        $bai_good and symlink($input_bai, $sample_bai);
        $flagstat_good and symlink($input_flagstat, $sample_flagstat);

        next if $bai_good and $flagstat_good;

        my $job_id = "PrepBam_${sample}_" . getJobId();
        my $log_dir = catfile($opt->{OUTPUT_DIR}, $sample, "logs");
        my $bash_file = catfile($opt->{OUTPUT_DIR}, $sample, "jobs", "$job_id.sh");

        from_template(
            "PrepBam.sh.tt", $bash_file,
            sample => $sample,
            sample_bam => $sample_bam,
            sample_bai => $sample_bai,
            sample_flagstat => $sample_flagstat,
            log_dir => $log_dir,
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "MAPPING");
        system "$qsub -o $log_dir/PrepBam_${sample}.out -e $log_dir/PrepBam_${sample}.err -N $job_id $bash_file";
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }
    return;
}

sub verifyBam {
    my ($bam_file, $opt) = @_;

    -f $bam_file or confess "ERROR: $bam_file does not exist";
    (my $sample = fileparse($bam_file)) =~ s/\.bam$//;

    my $headers = bamHeaders($bam_file, $opt);
    my @read_groups = grep { $_->{name} eq '@RG' } @$headers;

    my @sample_names = map { $_->{tags}{SM} } @read_groups;
    confess "too many samples in BAM $bam_file: @sample_names" if keys {map { $_ => 1 } @sample_names} > 1;
    warn "missing sample name in BAM $bam_file, using file name" unless @sample_names;
    $sample = $sample_names[0] if @sample_names;

    my %header_contigs = map { $_->{tags}{SN} => $_->{tags}{LN} } grep { $_->{name} eq '@SQ' } @$headers;
    verifyContigs(\%header_contigs, refGenomeContigs($opt));
    verifyReadGroups(\@read_groups, bamReads($bam_file, 1000, $opt));

    return $sample;
}

sub verifyBai {
    my ($bai_file, $bam_file, $opt) = @_;

    (-f $bai_file && -M $bam_file > -M $bai_file) or return 0;

    # this check does not happen if the .bai is missing/old, so re-implemented in the job :(
    my $headers = bamHeaders($bam_file, $opt);
    my %header_contigs = map { $_->{tags}{SN} => $_->{tags}{LN} } grep { $_->{name} eq '@SQ' } @$headers;
    verifyContigs(indexContigs($bam_file, $opt), \%header_contigs);
    return 1;
}

sub verifyFlagstat {
    my ($flagstat_file, $bam_file) = @_;

    return -f $flagstat_file && -M $bam_file > -M $flagstat_file;
}

sub bamHeaders {
    my ($bam_file, $opt) = @_;

    my $samtools = catfile($opt->{SAMTOOLS_PATH}, "samtools");

    my @lines = qx($samtools view -H $bam_file);
    $? == 0 or confess "could not read BAM headers from $bam_file";

    chomp @lines;
    my @fields = grep { @$_[0] ne '@CO' } map { [ split /\t/ ] } @lines;
    #<<< no perltidy
    my @headers = map {
        name => shift $_,
        tags => {map { split /:/, $_, 2 } @$_},
    }, @fields;
    #>>> no perltidy
    return \@headers;
}

sub bamReads {
    my ($bam_file, $num_lines, $opt) = @_;

    my $samtools = catfile($opt->{SAMTOOLS_PATH}, "samtools");

    my @lines = qx($samtools view $bam_file | head -$num_lines);
    chomp @lines;

    my @field_names = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual tags);
    my @fields = map { [ split "\t", $_, @field_names ] } @lines;
    my @reads = map { +{zip @field_names, @$_} } @fields;
    #<<< no perltidy
    map {
        $_->{tags} = {
            map {
                shift @$_ => {
                    type => shift @$_,
                    value => shift @$_,}
            } map { [ split ":" ] } split "\t", $_->{tags}}
    } @reads;
    #>>> no perltidy
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
    return;
}

sub verifyReadGroups {
    my ($read_groups, $bam_reads) = @_;

    my %header_rgids = map { $_->{tags}{ID} => 1 } @$read_groups;
    my %reads_rgids = map { $_->{tags}{RG}{value} => 1 } @$bam_reads;
    my @unknown_rgids = grep { not exists $header_rgids{$_} } keys %reads_rgids;
    confess "read group IDs from read tags not in BAM header:\n\t" . join "\n\t", @unknown_rgids if @unknown_rgids;
    return;
}

1;
