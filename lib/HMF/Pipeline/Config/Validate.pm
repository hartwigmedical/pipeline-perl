package HMF::Pipeline::Config::Validate;

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

    # pseudo-DSL to de-duplicate checks
    my $key_not_present = sub { not defined $_[1] and return "No $_[0] option found in config files" };
    my $missing_file = sub { &$key_not_present or not -f $_[1] and return "$_[0] file $_[1] does not exist" };
    my $missing_optional_file = sub { not &$key_not_present and &$missing_file };
    my $missing_optional_files = sub {
        join " and\n", grep { $_ } map { &$missing_optional_file($_[0], $_) } split /\t/, $_[1];
    };
    my $invalid_choice = sub {
        my ($choices) = @_;
        return sub {
            my ($key, $choice) = @_;
            defined $choice and not grep { /^$choice$/ } @{$choices} and return "$key must be one of " . join ", ", @{$choices};
        };
    };
    my $key_not_present_and_not_present = sub {
        my @alt_keys = @_;
        return sub {
            my ($key, $value, $opt) = @_;
            not defined $value and not grep { defined $opt->{$_} } @alt_keys and return "No $key or " . join(", ", @alt_keys) . " found in config files";
        };
    };
    my $if_enabled = sub {
        my ($more_checks) = @_;
        return sub {
            my ($key, $value, $opt) = @_;
            defined $value and $value eq "yes" and return applyChecks($more_checks, $opt);
        };
    };
    my $compare = sub {
        my ($other_key, $func, $name) = @_;
        return sub {
            my ($key, $value, $opt) = @_;
            my $other_value = $opt->{$other_key};
            &$key_not_present or not &$func($value, $other_value) and return "$key ($value) must be $name $other_key ($other_value)";
        };
    };

    return {
        INIFILE => $key_not_present,
        OUTPUT_DIR => $key_not_present,
        FASTQ => &$key_not_present_and_not_present('BAM'),
        MAIL => $key_not_present,
        CLUSTER_PATH => $key_not_present,
        CLUSTER_TMP => $key_not_present,
        CLUSTER_RESERVATION => $key_not_present,
        CLUSTER_PROJECT => $key_not_present,
        PRESTATS => $key_not_present,
        MAPPING => $key_not_present,
        POSTSTATS => $key_not_present,
        INDELREALIGNMENT => $key_not_present,
        BASEQUALITYRECAL => $key_not_present,
        VARIANT_CALLING => $key_not_present,
        FILTER_VARIANTS => $key_not_present,
        SOMATIC_VARIANTS => $key_not_present,
        COPY_NUMBER => $key_not_present,
        BAF => $key_not_present,
        ANNOTATE_VARIANTS => $key_not_present,
        KINSHIP => $key_not_present,
        FINALIZE => $key_not_present,
        GENOME => $missing_file,
        SAMBAMBA_PATH => $key_not_present,
        QUEUE_PATH => $key_not_present,
        PRESTATS => &$if_enabled({
                FASTQC_PATH => $key_not_present,
                PRESTATS_THREADS => $key_not_present,
                PRESTATS_MEM => $key_not_present,
                PRESTATS_QUEUE => $key_not_present,
                PRESTATS_TIME => $key_not_present,
            }
        ),
        MAPPING => &$if_enabled({
                BWA_PATH => $key_not_present,
                MAPPING_THREADS => $key_not_present,
                MAPPING_MEM => $key_not_present,
                MAPPING_QUEUE => $key_not_present,
                MAPPING_TIME => $key_not_present,
                MAPPING_SETTINGS => $key_not_present,

                MARKDUP_QUEUE => $key_not_present,
                MARKDUP_TIME => $key_not_present,
                MARKDUP_THREADS => $key_not_present,
                MARKDUP_MEM => $key_not_present,
                MARKDUP_OVERFLOW_LIST_SIZE => $key_not_present,
            }
        ),
        FLAGSTAT_QUEUE => $key_not_present,
        FLAGSTAT_THREADS => $key_not_present,
        FLAGSTAT_MEM => $key_not_present,
        FLAGSTAT_TIME => $key_not_present,
        POSTSTATS => &$if_enabled({
                BAMMETRICS_PATH => $key_not_present,
                PICARD_PATH => $key_not_present,
                POSTSTATS_THREADS => $key_not_present,
                POSTSTATS_MEM => $key_not_present,
                POSTSTATS_QUEUE => $key_not_present,
                POSTSTATS_TIME => $key_not_present,
                EXONCALLCOV => $key_not_present,
                EXONCALLCOV => &$if_enabled({
                        EXONCALLCOV_QUEUE => $key_not_present,
                        EXONCALLCOV_TIME => $key_not_present,
                        EXONCALLCOV_MEM => $key_not_present,
                        EXONCALLCOV_PATH => $key_not_present,
                        EXONCALLCOV_BED => $missing_file,
                        EXONCALLCOV_PREF => $missing_file,
                        EXONCALLCOV_PANEL => $missing_file,
                        EXONCALLCOV_ENS => $missing_file,
                    }
                ),
            }
        ),
        INDELREALIGNMENT => &$if_enabled({
                BAMUTIL_PATH => $key_not_present,
                REALIGNMENT_MASTER_QUEUE => $key_not_present,
                REALIGNMENT_MASTER_THREADS => $key_not_present,
                REALIGNMENT_MASTER_TIME => $key_not_present,
                REALIGNMENT_MASTER_MEM => $key_not_present,
                REALIGNMENT_QUEUE => $key_not_present,
                REALIGNMENT_THREADS => $key_not_present,
                REALIGNMENT_MEM => $key_not_present,
                REALIGNMENT_TIME => $key_not_present,
                REALIGNMENT_MERGETHREADS => $key_not_present,
                REALIGNMENT_SCALA => $key_not_present,
                REALIGNMENT_SCATTER => $key_not_present,
                REALIGNMENT_KNOWN => $missing_optional_files,
                FLAGSTAT_QUEUE => $key_not_present,
                FLAGSTAT_THREADS => $key_not_present,
                FLAGSTAT_MEM => $key_not_present,
                FLAGSTAT_TIME => $key_not_present,
            }
        ),
        BASEQUALITYRECAL => &$if_enabled({
                BASERECALIBRATION_MASTER_QUEUE => $key_not_present,
                BASERECALIBRATION_MASTER_TIME => $key_not_present,
                BASERECALIBRATION_MASTER_THREADS => $key_not_present,
                BASERECALIBRATION_MASTER_MEM => $key_not_present,
                BASERECALIBRATION_QUEUE => $key_not_present,
                BASERECALIBRATION_THREADS => $key_not_present,
                BASERECALIBRATION_MEM => $key_not_present,
                BASERECALIBRATION_TIME => $key_not_present,
                BASERECALIBRATION_SCALA => $key_not_present,
                BASERECALIBRATION_SCATTER => $key_not_present,
                BASERECALIBRATION_KNOWN => $missing_optional_files,
                FLAGSTAT_QUEUE => $key_not_present,
                FLAGSTAT_THREADS => $key_not_present,
                FLAGSTAT_MEM => $key_not_present,
                FLAGSTAT_TIME => $key_not_present,
            }
        ),
        VARIANT_CALLING => &$if_enabled({
                CALLING_MASTER_QUEUE => $key_not_present,
                CALLING_MASTER_TIME => $key_not_present,
                CALLING_MASTER_THREADS => $key_not_present,
                CALLING_MASTER_MEM => $key_not_present,
                CALLING_QUEUE => $key_not_present,
                CALLING_THREADS => $key_not_present,
                CALLING_MEM => $key_not_present,
                CALLING_TIME => $key_not_present,
                CALLING_SCATTER => $key_not_present,
                CALLING_GVCF => $key_not_present,
                CALLING_SCALA => $key_not_present,
                CALLING_UGMODE => &$invalid_choice([ "SNP", "INDEL", "BOTH" ]),
                CALLING_STANDCALLCONF => $key_not_present,
                CALLING_STANDEMITCONF => $key_not_present,
                CALLING_TARGETS => $missing_optional_file,
                CALLING_DBSNP => $missing_optional_file,
            }
        ),
        FILTER_VARIANTS => &$if_enabled({
                FILTER_MASTER_QUEUE => $key_not_present,
                FILTER_MASTER_TIME => $key_not_present,
                FILTER_MASTER_THREADS => $key_not_present,
                FILTER_MASTER_MEM => $key_not_present,
                FILTER_QUEUE => $key_not_present,
                FILTER_THREADS => $key_not_present,
                FILTER_MEM => $key_not_present,
                FILTER_TIME => $key_not_present,
                FILTER_SCATTER => $key_not_present,
                FILTER_SCALA => $key_not_present,
                FILTER_SNPTYPES => $key_not_present,
                FILTER_SNPNAME => $key_not_present,
                FILTER_SNPEXPR => $key_not_present,
                FILTER_INDELTYPES => $key_not_present,
                FILTER_INDELNAME => $key_not_present,
                FILTER_INDELEXPR => $key_not_present,
            }
        ),
        SOMATIC_VARIANTS => &$if_enabled({
                VCFTOOLS_PATH => $key_not_present,
                SAMTOOLS_PATH => $key_not_present,
                SOMVAR_TARGETS => $missing_optional_file,
                SOMVAR_STRELKA => $key_not_present,
                SOMVAR_STRELKA => &$if_enabled({
                        STRELKA_PATH => $key_not_present,
                        STRELKA_INI => $key_not_present,
                        STRELKA_QUEUE => $key_not_present,
                        STRELKA_THREADS => $key_not_present,
                        STRELKA_MEM => $key_not_present,
                        STRELKA_TIME => $key_not_present,
                    }
                ),
                SOMVAR_VARSCAN => $key_not_present,
                SOMVAR_VARSCAN => &$if_enabled({
                        VARSCAN_PATH => $key_not_present,
                        PBGZIP_PATH => $key_not_present,
                        TABIX_PATH => $key_not_present,
                        VARSCAN_QUEUE => $key_not_present,
                        VARSCAN_THREADS => $key_not_present,
                        VARSCAN_TIME => $key_not_present,
                        VARSCAN_MEM => $key_not_present,
                        VARSCAN_SETTINGS => $key_not_present,
                        VARSCAN_POSTSETTINGS => $key_not_present,
                        PILEUP_QUEUE => $key_not_present,
                        PILEUP_DIVISOR => $key_not_present,
                        PILEUP_THREADS => &$compare('PILEUP_DIVISOR', sub { $_[0] >= $_[1] }, "greater than or equal to"),
                        PILEUP_MEM => $key_not_present,
                        PILEUP_TIME => $key_not_present,
                        FINALIZE_KEEP_PILEUP => $key_not_present,
                    }
                ),
                SOMVAR_FREEBAYES => $key_not_present,
                SOMVAR_FREEBAYES => &$if_enabled({
                        FREEBAYES_PATH => $key_not_present,
                        VCFLIB_PATH => $key_not_present,
                        FREEBAYES_QUEUE => $key_not_present,
                        FREEBAYES_THREADS => $key_not_present,
                        FREEBAYES_MEM => $key_not_present,
                        FREEBAYES_TIME => $key_not_present,
                        FREEBAYES_SETTINGS => $key_not_present,
                        FREEBAYES_SOMATICFILTER => $key_not_present,
                    }
                ),
                SOMVAR_MUTECT => $key_not_present,
                SOMVAR_MUTECT => &$if_enabled({
                        MUTECT_PATH => $key_not_present,
                        MUTECT_QUEUE => $key_not_present,
                        MUTECT_THREADS => $key_not_present,
                        MUTECT_MEM => $key_not_present,
                        MUTECT_TIME => $key_not_present,
                        MUTECT_COSMIC => $key_not_present,
                    }
                ),
                SOMVARMERGE_QUEUE => $key_not_present,
                SOMVARMERGE_THREADS => $key_not_present,
                SOMVARMERGE_MEM => $key_not_present,
                SOMVARMERGE_TIME => $key_not_present,
                SOMVAR_ANNOTATE => $key_not_present,
                SOMVAR_ANNOTATE => &$if_enabled({
                        ANNOTATE_DB => $key_not_present,
                        ANNOTATE_FLAGS => $key_not_present,
                        ANNOTATE_IDNAME => $key_not_present,
                        ANNOTATE_IDDB => $missing_file,
                        CALLING_DBSNP => $missing_file,
                    }
                ),
            }
        ),
        COPY_NUMBER => &$if_enabled({
                CNVCHECK_QUEUE => $key_not_present,
                CNVCHECK_THREADS => $key_not_present,
                CNVCHECK_MEM => $key_not_present,
                CNVCHECK_TIME => $key_not_present,
                CNV_MODE => $key_not_present,
                CNV_FREEC => $key_not_present,
                CNV_FREEC => &$if_enabled({
                        FREEC_PATH => $key_not_present,
                        FREEC_QUEUE => $key_not_present,
                        FREEC_THREADS => $key_not_present,
                        FREEC_MEM => $key_not_present,
                        FREEC_TIME => $key_not_present,
                        FREEC_CHRLENFILE => $key_not_present,
                        FREEC_CHRFILES => $key_not_present,
                        FREEC_PLOIDY => $key_not_present,
                        FREEC_WINDOW => $key_not_present,
                        FREEC_TELOCENTROMERIC => $key_not_present,
                    }
                ),
                CNV_QDNASEQ => $key_not_present,
                CNV_QDNASEQ => &$if_enabled({
                        QDNASEQ_PATH => $key_not_present,
                        QDNASEQ_QUEUE => $key_not_present,
                        QDNASEQ_THREADS => $key_not_present,
                        QDNASEQ_MEM => $key_not_present,
                        QDNASEQ_TIME => $key_not_present,
                    }
                ),
            }
        ),
        BAF => &$if_enabled({
                BAF_QUEUE => $key_not_present,
                BAF_THREADS => $key_not_present,
                BAF_MEM => $key_not_present,
                BAF_TIME => $key_not_present,
                BIOVCF_PATH => $key_not_present,
                BAF_SNPS => $key_not_present,
            }
        ),
        ANNOTATE_VARIANTS => &$if_enabled({
                SNPEFF_PATH => $key_not_present,
                IGVTOOLS_PATH => $key_not_present,
                ANNOTATE_QUEUE => $key_not_present,
                ANNOTATE_THREADS => $key_not_present,
                ANNOTATE_MEM => $key_not_present,
                ANNOTATE_TIME => $key_not_present,
                ANNOTATE_SNPEFF => $key_not_present,
                ANNOTATE_SNPEFF => &$if_enabled({
                        ANNOTATE_DB => $key_not_present,
                        ANNOTATE_FLAGS => $key_not_present,
                    }
                ),
                ANNOTATE_SNPSIFT => $key_not_present,
                ANNOTATE_SNPSIFT => &$if_enabled({
                        ANNOTATE_DBNSFP => $missing_file,
                        ANNOTATE_FIELDS => $key_not_present,
                    }
                ),
                ANNOTATE_FREQUENCIES => $key_not_present,
                ANNOTATE_FREQUENCIES => &$if_enabled({
                        ANNOTATE_FREQNAME => $key_not_present,
                        ANNOTATE_FREQDB => $missing_file,
                        ANNOTATE_FREQINFO => $key_not_present,
                    }
                ),
                ANNOTATE_IDFIELD => $key_not_present,
                ANNOTATE_IDFIELD => &$if_enabled({
                        ANNOTATE_IDNAME => $key_not_present,
                        ANNOTATE_IDDB => $missing_file,
                    }
                ),
            }
        ),
        KINSHIP => &$if_enabled({
                KINSHIP_QUEUE => $key_not_present,
                KINSHIP_THREADS => $key_not_present,
                KINSHIP_MEM => $key_not_present,
                KINSHIP_TIME => $key_not_present,
                PLINK_PATH => $key_not_present,
                KING_PATH => $key_not_present,
                VCFTOOLS_PATH => $key_not_present,
            }
        ),
        FINALIZE => &$if_enabled({
                FINALIZE_QUEUE => $key_not_present,
                FINALIZE_THREADS => $key_not_present,
                FINALIZE_MEM => $key_not_present,
                FINALIZE_TIME => $key_not_present,
            }
        ),
    };
}

1;
