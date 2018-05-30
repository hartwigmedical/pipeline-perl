package HMF::Pipeline::Functions::Validate;

use FindBin::libs;
use discipline;

use Carp;
use File::Basename;
use File::Spec::Functions;
use List::MoreUtils qw(zip);

use parent qw(Exporter);
our @EXPORT_OK = qw(parseFastqName verifyConfig verifyBam verifyBai verifyFlagstat);

sub parseFastqName {
    my ($input_file) = @_;
    my $name = fileparse($input_file);
    (my $R1 = $input_file) =~ s/_R2/_R1/;
    (my $R2 = $input_file) =~ s/_R1/_R2/;

    my $fastQPattern = qr/^(?<sampleName>[^-_][^_]*)
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

sub verifyBam {
    my ($bam_file, $opt) = @_;

    -f $bam_file or confess "ERROR: $bam_file does not exist";
    (my $sample = fileparse($bam_file)) =~ s/\.bam$//;

    my $headers = bamHeaders($bam_file, $opt);
    my @read_groups = grep { $_->{name} eq '@RG' } @$headers;

    my @sample_names = map { $_->{tags}{SM} } grep { $_->{tags}->{SM} } @read_groups;
    confess "too many samples in BAM $bam_file: @sample_names" if keys %{{map { $_ => 1 } @sample_names}} > 1;
    warn "missing sample name (\@RG SM tag) in BAM $bam_file, using file name" unless @sample_names;
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
    $? == 0 or confess "could not parse BAM headers from $bam_file";

    chomp @lines;
    my @fields = grep { @$_[0] ne '@CO' } map { [ split /\t/ ] } @lines;
    #<<< no perltidy
    my @headers = map {
        name => shift @$_,
        tags => {map { split /:/, $_, 2 } @$_},
    }, @fields;
    #>>> no perltidy
    return \@headers;
}

sub bamReads {
    my ($bam_file, $num_lines, $opt) = @_;

    my $samtools = catfile($opt->{SAMTOOLS_PATH}, "samtools");

    my @lines = qx(bash -c "set -o pipefail; $samtools view $bam_file | head -$num_lines");
    $? == 0 or ($? >> 8) == 141 or confess "could not parse BAM reads from $bam_file";
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

    my $fai_file = "$opt->{REF_GENOME}.fai";
    my @lines = qx(cat $fai_file);
    $? == 0 or confess "could not read from $fai_file";

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
    my %contigs = map { splice @{[ split "\t" ]}, 0, 2 } @$lines;
    return \%contigs;
}

sub verifyContigs {
    my ($contigs, $ref_contigs) = @_;

    my @warnings;
    foreach my $key (sort keys %{$contigs}) {
        my $value = $contigs->{$key};
        push @warnings, "contig $key missing" and next if not exists $ref_contigs->{$key};
        push @warnings, "contig $key length $value does not match $ref_contigs->{$key}" if $value != $ref_contigs->{$key};
    }
    warn $_ foreach @warnings;
    confess "contigs do not match" if @warnings;
    return 1;
}

sub verifyReadGroups {
    my ($read_groups, $bam_reads) = @_;

    my %header_rgids = map { $_->{tags}{ID} => 1 } @$read_groups;
    my %reads_rgids = map { $_->{tags}{RG}{value} => 1 } @$bam_reads;
    my @unknown_rgids = grep { not exists $header_rgids{$_} } keys %reads_rgids;
    confess "read group IDs from read tags not in BAM header:\n\t" . join "\n\t", @unknown_rgids if @unknown_rgids;
    return 1;
}

sub verifyConfig {
    my ($opt) = @_;

    return applyChecks(configChecks(), $opt);
}

sub applyChecks {
    my ($checks, $opt) = @_;

    my @errors;
    while (my ($key, $checker) = each %{$checks}) {
        my $error = &$checker($key, $opt->{$key}, $opt);
        if (ref $error eq "ARRAY") {
            push @errors, @{$error};
        } else {
            push @errors, $error if $error;
        }
    }
    return \@errors;
}

sub configChecks {
    return {
        INIFILE => \&key_not_present,
        OUTPUT_DIR => \&key_not_present,
        FASTQ => key_not_present_and_not_present("BAM"),
        CLUSTER_PATH => \&missing_directory,
        CLUSTER_TMP => \&key_not_present,
        CLUSTER_RESERVATION => \&key_not_present,
        CLUSTER_PROJECT => \&key_not_present,

        REF_GENOME => \&missing_genome_files,

        # KODU: GATK and SAM tools are used all over the place so check them regardless of configuration.
        GATK_PATH => \&missing_directory,
        QUEUE_PATH => \&missing_directory,
        SAMTOOLS_PATH => \&missing_directory,

        # KODU: Sambamba creates bam indices regardless of how pipeline is run so we always need below params.
        SAMBAMBA_PATH => \&missing_directory,
        INDEXING_THREADS => \&key_not_present,

        PRESTATS => if_enabled({
                FASTQC_PATH => \&missing_directory,
                PRESTATS_MEM => \&key_not_present,
                PRESTATS_QUEUE => \&key_not_present,
                PRESTATS_THREADS => \&key_not_present,
                PRESTATS_TIME => \&key_not_present,
            }
        ),
        MAPPING => if_enabled({
                BWA_PATH => \&missing_directory,
                MAPPING_QUEUE => \&key_not_present,
                MAPPING_TIME => \&key_not_present,
                MAPPING_THREADS => \&key_not_present,
                MAPPING_MEM => \&key_not_present,
                MAPPING_SETTINGS => \&key_not_present,

                MARKDUP_QUEUE => \&key_not_present,
                MARKDUP_TIME => \&key_not_present,
                MARKDUP_THREADS => \&key_not_present,
                MARKDUP_MEM => \&key_not_present,
                MARKDUP_OVERFLOW_LIST_SIZE => \&key_not_present,

                FLAGSTAT_QUEUE => \&key_not_present,
                FLAGSTAT_THREADS => \&key_not_present,
                FLAGSTAT_MEM => \&key_not_present,
                FLAGSTAT_TIME => \&key_not_present,
            }
        ),
        POSTSTATS => if_enabled({
                BAMMETRICS_PATH => \&missing_directory,
                PICARD_PATH => \&missing_directory,

                POSTSTATS_QUEUE => \&key_not_present,
                POSTSTATS_TIME => \&key_not_present,
                POSTSTATS_THREADS => \&key_not_present,
                POSTSTATS_MEM => \&key_not_present,
            }
        ),
        DAMAGE_ESTIMATE => if_enabled({
                DAMAGE_ESTIMATOR_PATH => \&missing_directory,
                DAMAGE_ESTIMATE_QUEUE => \&key_not_present,
                DAMAGE_ESTIMATE_TIME => \&key_not_present,
                DAMAGE_ESTIMATE_THREADS => \&key_not_present,
                DAMAGE_ESTIMATE_MEM => \&key_not_present,
                DAMAGE_ESTIMATE_DOWNSAMPLE_THREADS => \&key_not_present,
                DAMAGE_ESTIMATE_DOWNSAMPLE_BAM_SIZE => \&key_not_present,
                DAMAGE_ESTIMATE_MIN_COVERAGE_LIMIT => \&key_not_present,
                DAMAGE_ESTIMATE_MAX_COVERAGE_LIMIT => \&key_not_present,
            }
        ),
        INDEL_REALIGNMENT => if_enabled({
                BAMUTIL_PATH => \&missing_directory,

                REALIGNMENT_MASTER_QUEUE => \&key_not_present,
                REALIGNMENT_MASTER_THREADS => \&key_not_present,
                REALIGNMENT_MASTER_TIME => \&key_not_present,
                REALIGNMENT_MASTER_MEM => \&key_not_present,
                REALIGNMENT_QUEUE => \&key_not_present,
                REALIGNMENT_THREADS => \&key_not_present,
                REALIGNMENT_MEM => \&key_not_present,
                REALIGNMENT_TIME => \&key_not_present,

                REALIGNMENT_MERGETHREADS => \&key_not_present,
                REALIGNMENT_SCALA => \&key_not_present,
                REALIGNMENT_SCATTER => \&key_not_present,
                REALIGNMENT_KNOWN => \&missing_optional_files,

                FLAGSTAT_QUEUE => \&key_not_present,
                FLAGSTAT_THREADS => \&key_not_present,
                FLAGSTAT_MEM => \&key_not_present,
                FLAGSTAT_TIME => \&key_not_present,
            }
        ),
        GERMLINE_CALLING => if_enabled({
                GERMLINE_CALLING_MASTER_QUEUE => \&key_not_present,
                GERMLINE_CALLING_MASTER_TIME => \&key_not_present,
                GERMLINE_CALLING_MASTER_THREADS => \&key_not_present,
                GERMLINE_CALLING_MASTER_MEM => \&key_not_present,

                GERMLINE_CALLING_QUEUE => \&key_not_present,
                GERMLINE_CALLING_TIME => \&key_not_present,
                GERMLINE_CALLING_THREADS => \&key_not_present,
                GERMLINE_CALLING_MEM => \&key_not_present,

                GERMLINE_CALLING_SCATTER => \&key_not_present,
                GERMLINE_CALLING_SCALA => \&key_not_present,
                GERMLINE_CALLING_UGMODE => invalid_choice([ "SNP", "INDEL", "BOTH" ]),
                GERMLINE_CALLING_STANDCALLCONF => \&key_not_present,
                DBSNP_DB => \&missing_optional_file,
                GERMLINE_FILTER_MASTER_QUEUE => \&key_not_present,
                GERMLINE_FILTER_MASTER_TIME => \&key_not_present,
                GERMLINE_FILTER_MASTER_THREADS => \&key_not_present,
                GERMLINE_FILTER_MASTER_MEM => \&key_not_present,
                GERMLINE_FILTER_QUEUE => \&key_not_present,
                GERMLINE_FILTER_THREADS => \&key_not_present,
                GERMLINE_FILTER_MEM => \&key_not_present,
                GERMLINE_FILTER_TIME => \&key_not_present,
                GERMLINE_FILTER_SCATTER => \&key_not_present,
                GERMLINE_FILTER_SCALA => \&key_not_present,
                GERMLINE_FILTER_SNPTYPES => \&key_not_present,
                GERMLINE_FILTER_SNPNAME => \&key_not_present,
                GERMLINE_FILTER_SNPEXPR => \&key_not_present,
                GERMLINE_FILTER_INDELTYPES => \&key_not_present,
                GERMLINE_FILTER_INDELNAME => \&key_not_present,
                GERMLINE_FILTER_INDELEXPR => \&key_not_present,
                SNPEFF_PATH => \&missing_directory,
                TABIX_PATH => \&missing_directory,
                ANNOTATE_QUEUE => \&key_not_present,
                ANNOTATE_THREADS => \&key_not_present,
                ANNOTATE_MEM => \&key_not_present,
                ANNOTATE_TIME => \&key_not_present,
                ANNOTATE_SNPEFF_DB => \&key_not_present,
                ANNOTATE_SNPEFF_FLAGS => \&key_not_present,
                ANNOTATE_DBNSFP_DB => \&missing_file,
                ANNOTATE_DBNSFP_FIELDS => \&key_not_present,
                ANNOTATE_COSMIC_DB => \&missing_file,
                ANNOTATE_COSMIC_ID => \&key_not_present,
                ANNOTATE_FREQ_DB => \&missing_file,
                ANNOTATE_FREQ_NAME => \&key_not_present,
                ANNOTATE_FREQ_INFO => \&key_not_present,
            }
        ),
        AMBER => if_enabled({
                AMBER_PATH => \&missing_directory,
                AMBER_QUEUE => \&key_not_present,
                AMBER_TIME => \&key_not_present,
                AMBER_THREADS => \&key_not_present,
                AMBER_MEM => \&key_not_present,
                BAF_SNPS => \&missing_file,
            }
        ),
        COBALT => if_enabled({
                COBALT_PATH => \&missing_directory,
                COBALT_QUEUE => \&key_not_present,
                COBALT_TIME => \&key_not_present,
                COBALT_THREADS => \&key_not_present,
                COBALT_MEM => \&key_not_present,
            }
        ),
        PURPLE => if_enabled({
                PURPLE_PATH => \&missing_directory,
                PURPLE_QUEUE => \&key_not_present,
                PURPLE_TIME => \&key_not_present,
                PURPLE_THREADS => \&key_not_present,
                PURPLE_MEM => \&key_not_present,
                CIRCOS_PATH => \&missing_directory,
                GC_PROFILE => \&missing_file,
            }
        ),
        STRELKA => if_enabled({
                QUEUE_LOW_GZIP_COMPRESSION_PATH => \&missing_directory,
                BAMUTIL_PATH => \&missing_directory,

                BASERECALIBRATION_MASTER_QUEUE => \&key_not_present,
                BASERECALIBRATION_MASTER_TIME => \&key_not_present,
                BASERECALIBRATION_MASTER_THREADS => \&key_not_present,
                BASERECALIBRATION_MASTER_MEM => \&key_not_present,
                BASERECALIBRATION_QUEUE => \&key_not_present,
                BASERECALIBRATION_THREADS => \&key_not_present,
                BASERECALIBRATION_MEM => \&key_not_present,
                BASERECALIBRATION_TIME => \&key_not_present,
                BASERECALIBRATION_SCALA => \&key_not_present,
                BASERECALIBRATION_SCATTER => \&key_not_present,
                BASERECALIBRATION_KNOWN => \&missing_optional_files,

                FLAGSTAT_QUEUE => \&key_not_present,
                FLAGSTAT_THREADS => \&key_not_present,
                FLAGSTAT_MEM => \&key_not_present,
                FLAGSTAT_TIME => \&key_not_present,

                STRELKA_PATH => \&missing_directory,
                STRELKA_POST_PROCESS_PATH => \&missing_directory,
                STRELKA_INI => \&key_not_present,
                STRELKA_QUEUE => \&key_not_present,
                STRELKA_THREADS => \&key_not_present,
                STRELKA_MEM => \&key_not_present,
                STRELKA_TIME => \&key_not_present,

                HIGH_CONFIDENCE_BED => \&missing_file,
                STRELKAPOSTPROCESS_QUEUE => \&key_not_present,
                STRELKAPOSTPROCESS_THREADS => \&key_not_present,
                STRELKAPOSTPROCESS_MEM => \&key_not_present,
                STRELKAPOSTPROCESS_TIME => \&key_not_present,

                ANNOTATE_SNPEFF_DB => \&key_not_present,
                ANNOTATE_SNPEFF_FLAGS => \&key_not_present,
                ANNOTATE_COSMIC_DB => \&missing_file,
                ANNOTATE_COSMIC_ID => \&key_not_present,
                DBSNP_DB => \&missing_file,
                HMF_GERMLINE_PON => \&missing_file,
                HMF_SOMATIC_PON => \&missing_file,
                MAPPABILITY_BED => \&missing_file,
                MAPPABILITY_HDR => \&missing_file,
                HMF_SOMATIC_HOTSPOTS => \&missing_file,

                TABIX_PATH => \&missing_directory,
            }
        ),
        STRUCTURAL_VARIANT_CALLING => if_enabled({
                MANTA => if_enabled({
                        MANTA_PATH => \&missing_directory,
                        MANTA_QUEUE => \&key_not_present,
                        MANTA_THREADS => \&key_not_present,
                        MANTA_MEM => \&key_not_present,
                        MANTA_TIME => \&key_not_present,

                        BPI_PATH => \&missing_directory,
                        BPI_QUEUE => \&key_not_present,
                        BPI_THREADS => \&key_not_present,
                        BPI_MEM => \&key_not_present,
                        BPI_TIME => \&key_not_present,

                        TABIX_PATH => \&missing_directory,
                    }
                ),
            }
        ),
        HEALTHCHECK => if_enabled({
                HEALTHCHECK_QUEUE => \&key_not_present,
                HEALTHCHECK_THREADS => \&key_not_present,
                HEALTHCHECK_MEM => \&key_not_present,
                HEALTHCHECK_TIME => \&key_not_present,
            }
        ),
        FINALIZE_THREADS => \&key_not_present,
        FINALIZE_QUEUE => \&key_not_present,
        FINALIZE_MEM => \&key_not_present,
        FINALIZE_TIME => \&key_not_present,
    };
}

sub if_enabled {
    my ($more_checks) = @_;

    return sub {
        my ($key, $value, $opt) = @_;
        key_not_present($key, $value) or $value eq "yes" and return applyChecks($more_checks, $opt);
    };
}

sub key_not_present {
    my ($key, $value) = @_;

    not defined $value and return "No $key option found in config files";
    return;
}

sub missing_file {
    my ($key, $value) = @_;

    return (key_not_present($key, $value) or not -f $value and "$key file $value does not exist");
}

sub missing_files {
    my ($key, $value) = @_;

    return (key_not_present($key, $value) or join " and\n", grep { $_ } map { missing_file($key, $_) } split /\t/, $value);
}

sub missing_optional_file {
    my ($key, $value) = @_;

    return (not key_not_present($key, $value) and missing_file($key, $value));
}

sub missing_optional_files {
    my ($key, $value) = @_;

    return (not key_not_present($key, $value) and join " and\n", grep { $_ } map { missing_optional_file($key, $_) } split /\t/, $value);
}

sub missing_directory {
    my ($key, $value) = @_;

    return (key_not_present($key, $value) or not -d $value and "$key directory $value does not exist");
}

sub missing_genome_files {
    my ($key, $value) = @_;

    return (key_not_present($key, $value) or join " and ", grep { $_ } map { missing_file($key, $_) } ("${value}", "${value}.fai", "${value}.bwt"));
}

sub invalid_choice {
    my ($choices) = @_;

    return sub {
        my ($key, $choice) = @_;
        defined $choice and not grep { /^$choice$/ } @{$choices} and return "$key must be one of " . join ", ", @{$choices};
    };
}

sub key_not_present_and_not_present {
    my (@alt_keys) = @_;

    return sub {
        my ($key, $value, $opt) = @_;
        not defined $value and not grep { defined $opt->{$_} } @alt_keys and return "No $key or " . join(", ", @alt_keys) . " found in config files";
    };
}

sub compare_to {
    my ($other_key, $func, $name) = @_;

    return sub {
        my ($key, $value, $opt) = @_;
        my $other_value = $opt->{$other_key};
        key_not_present($key, $value)
            or key_not_present($other_key, $other_value)
            or not &$func($value, $other_value)
            and return "$key ($value) must be $name $other_key ($other_value)";
    };
}

1;
