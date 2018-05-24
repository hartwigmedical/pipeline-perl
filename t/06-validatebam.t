#!/usr/bin/env perl

use discipline;

use File::Copy qw(move);
use File::Spec::Functions;
use File::Temp;
use File::Touch;
use Test::Fatal;
use Test::More;
use Test::Output;
use Test::Warn;

use lib "t";
# bug in Test::Prereq 1.x needs filename for test dependencies
require "Util.pm"; ## no critic (Modules::RequireBarewordIncludes)

use HMF::Pipeline::Functions::Validate;


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

SKIP: {
    skip "no SAMTOOLS_PATH set", 64 if not $ENV{SAMTOOLS_PATH};

    my ($exception, $bam_path, $headers, $reads, $sample, $result);

    $bam_path = Util::convertToBam(catfile("t", "data", "empty.sam"));
    $sample = HMF::Pipeline::Functions::Validate::verifyBam(
        $bam_path, {
            SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
            REF_GENOME => catfile("t", "data", "empty.fasta"),
        }
    );
    is($sample, "empty", "validates empty bam");

    $bam_path = Util::convertToBam(catfile("t", "data", "sample.sam"));
    $sample = HMF::Pipeline::Functions::Validate::verifyBam(
        $bam_path, {
            SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
            REF_GENOME => catfile("t", "data", "empty.fasta"),
        }
    );
    is($sample, "samplename", "sample name from bam");

    $bam_path = Util::convertToBam(catfile("t", "data", "nosamplename.sam"));
    warning_is {
        $sample = HMF::Pipeline::Functions::Validate::verifyBam(
            $bam_path, {
                SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
                REF_GENOME => catfile("t", "data", "empty.fasta"),
            }
            )
    }
    "missing sample name (\@RG SM tag) in BAM t/data/nosamplename.bam, using file name", "warns for missing SM tag";
    is($sample, "nosamplename", "uses file name for missing SM tag");

    $bam_path = Util::convertToBam(catfile("t", "data", "twosamplenames.sam"));
    $exception = exception {
        $sample = HMF::Pipeline::Functions::Validate::verifyBam(
            $bam_path, {
                SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
                REF_GENOME => catfile("t", "data", "empty.fasta"),
            }
            )
    };
    like($exception, qr(^too many samples in BAM t/data/twosamplenames.bam: sampleone sampletwo), "detects more than one SM tag");

    $bam_path = Util::convertToBam(catfile("t", "data", "extracontig.sam"));
    warning_is {
        $exception = exception {
            $sample = HMF::Pipeline::Functions::Validate::verifyBam(
                $bam_path, {
                    SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
                    REF_GENOME => catfile("t", "data", "empty.fasta"),
                }
                )
        }
    }
    "contig 3 missing", "warning for failed contig verification";
    like($exception, qr/^contigs do not match/, "detects contig not in ref genome");

    $bam_path = Util::convertToBam(catfile("t", "data", "extrareadgroup.sam"));
    $exception = exception {
        $sample = HMF::Pipeline::Functions::Validate::verifyBam(
            $bam_path, {
                SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
                REF_GENOME => catfile("t", "data", "empty.fasta"),
            }
            )
    };
    like($exception, qr/^read group IDs from read tags not in BAM header:\n\tmissing/, "detects read group not in header");

    $bam_path = Util::convertToBam(catfile("t", "data", "sample.sam"), 1);
    is(
        HMF::Pipeline::Functions::Validate::verifyBai(
            "${bam_path}.bai",
            $bam_path, {
                SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
                REF_GENOME => catfile("t", "data", "empty.fasta"),
            }
        ),
        1,
        "validates BAI"
    );

    # required for stupid filesystem timestamp granularity
    sleep 1;
    touch($bam_path);
    is(
        HMF::Pipeline::Functions::Validate::verifyBai(
            "${bam_path}.bai",
            $bam_path, {
                SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
                REF_GENOME => catfile("t", "data", "empty.fasta"),
            }
        ),
        0,
        "detects BAI too old"
    );

    unlink("${bam_path}.bai");
    is(
        HMF::Pipeline::Functions::Validate::verifyBai(
            "${bam_path}.bai",
            $bam_path, {
                SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
                REF_GENOME => catfile("t", "data", "empty.fasta"),
            }
        ),
        0,
        "detects BAI missing"
    );

    # TODO: {
    #        local $TODO = "samtools idxstats re-maps indices and cannot fail";

    #        my $basic_bam_path = Util::convertToBam(catfile("t", "data", "sample.sam"));
    #        $bam_path = Util::convertToBam(catfile("t", "data", "extracontig.sam"), 1);
    #        move "${bam_path}.bai", "${basic_bam_path}.bai";
    #        warning_is {
    #            $exception = exception {
    #                HMF::Pipeline::Functions::Validate::verifyBai(
    #                    "${basic_bam_path}.bai",
    #                    $basic_bam_path, {
    #                        SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
    #                        REF_GENOME => catfile("t", "data", "empty.fasta"),
    #                    }
    #                )
    #                }
    #        }
    #            "contig 3 missing", "warning for failed contig verification";
    #        like($exception, qr/^contigs do not match/, "detects BAI index contig not in header");
    #    }

    $bam_path = Util::convertToBam(catfile("t", "data", "sample.sam"), 1);
    (my $flagstat_path = $bam_path) =~ s/\.bam$/.flagstat/;
    ok(HMF::Pipeline::Functions::Validate::verifyFlagstat($flagstat_path, $bam_path), "validates flagstat");

    # required for stupid filesystem timestamp granularity
    sleep 1;
    touch($bam_path);
    ok(!HMF::Pipeline::Functions::Validate::verifyFlagstat($flagstat_path, $bam_path), "detects flagstat too old");

    unlink($flagstat_path);
    ok(!HMF::Pipeline::Functions::Validate::verifyFlagstat($flagstat_path, $bam_path), "detects flagstat missing");

    $bam_path = Util::convertToBam(catfile("t", "data", "sample.sam"));
    $headers = HMF::Pipeline::Functions::Validate::bamHeaders($bam_path, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}});
    is_deeply(
        $headers, [ {
                name => '@SQ',
                tags => {SN => 1, LN => 0}
            }, {
                name => '@RG',
                tags => {ID => "id", SM => "samplename"}
            },
        ],
        "parses BAM headers"
    );

    $bam_path = Util::convertToBam(catfile("t", "data", "reads.sam"));
    $reads = HMF::Pipeline::Functions::Validate::bamReads($bam_path, 100, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}});
    is_deeply(
        $reads, [ {
                qname => "seq",
                flag => "0",
                rname => "1",
                pos => "1",
                mapq => "30",
                cigar => "1=",
                rnext => "*",
                pnext => "0",
                tlen => "0",
                seq => "=",
                qual => "*",
                tags => {RG => {type => "Z", value => "id1"}},
            }, {
                qname => "seq",
                flag => "0",
                rname => "1",
                pos => "2",
                mapq => "30",
                cigar => "1=",
                rnext => "*",
                pnext => "0",
                tlen => "0",
                seq => "=",
                qual => "*",
                tags => {RG => {type => "Z", value => "id2"}},
            },
        ],
        "parses BAM reads"
    );

    $reads = HMF::Pipeline::Functions::Validate::bamReads($bam_path, 1, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}});
    is_deeply(
        $reads, [ {
                qname => "seq",
                flag => "0",
                rname => "1",
                pos => "1",
                mapq => "30",
                cigar => "1=",
                rnext => "*",
                pnext => "0",
                tlen => "0",
                seq => "=",
                qual => "*",
                tags => {RG => {type => "Z", value => "id1"}},
            },
        ],
        "limits number of BAM reads"
    );

    $bam_path = Util::convertToBam(catfile("t", "data", "morereads.sam"));
    $reads = HMF::Pipeline::Functions::Validate::bamReads($bam_path, 1, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}});
    is_deeply(
        $reads, [ {
                qname => "seq",
                flag => "0",
                rname => "1",
                pos => "1",
                mapq => "30",
                cigar => "1=",
                rnext => "*",
                pnext => "0",
                tlen => "0",
                seq => "=",
                qual => "*",
                tags => {RG => {type => "Z", value => "id1"}},
            },
        ],
        "limits number of BAM reads from longer file without failing due to SIGPIPE"
    );

    $bam_path = Util::convertToBam(catfile("t", "data", "reads.sam"), 1);
    my $contigs = HMF::Pipeline::Functions::Validate::indexContigs($bam_path, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}});
    is_deeply(
        $contigs, {
            1 => "2",
            2 => "3",
        },
        "parses BAI contigs"
    );

    my $contig_lengths = HMF::Pipeline::Functions::Validate::contigLengths(
        [
            "1	123\n", #
            "2	456\n", #
            "3	789",   #
        ]
    );
    is_deeply(
        $contig_lengths, {
            1 => "123",
            2 => "456",
            3 => "789",
        },
        "parses contig lengths"
    );

    warnings_are {
        $exception = exception {
            $result = HMF::Pipeline::Functions::Validate::verifyContigs({
                    1 => "42",
                    2 => "420",
                }, {
                    1 => "42",
                    2 => "420",
                }
                )
        }
    }
    [], "no warnings for matching contigs";
    is($exception, undef, "no exception for matching contigs");
    is($result, 1, "success for matching contigs");

    warnings_are {
        $exception = exception {
            HMF::Pipeline::Functions::Validate::verifyContigs({
                    1 => "42",
                    2 => "420",
                }, {
                    1 => "42",
                    2 => "420",
                    3 => "4200",
                }
                )
        }
    }
    [], "no warnings for extra ref contigs";
    is($exception, undef, "no exception for extra ref contigs");
    is($result, 1, "success for extra ref contigs");

    warning_is {
        $exception = exception {
            HMF::Pipeline::Functions::Validate::verifyContigs({
                    1 => "42",
                    2 => "420",
                }, {
                    1 => "42",
                }
                )
        }
    }
    "contig 2 missing", "warning for missing contig";
    like($exception, qr/^contigs do not match/, "detects overall contig mismatch");

    warning_is {
        $exception = exception {
            HMF::Pipeline::Functions::Validate::verifyContigs({
                    1 => "41",
                    2 => "420",
                }, {
                    1 => "42",
                    2 => "420",
                }
                )
        }
    }
    "contig 1 length 41 does not match 42", "warning for wrong contig length";
    like($exception, qr/^contigs do not match/, "detects overall contig mismatch");

    warnings_are {
        $exception = exception {
            HMF::Pipeline::Functions::Validate::verifyContigs({
                    1 => "41",
                    2 => "420",
                    3 => "4200",
                }, {
                    1 => "42",
                    2 => "420",
                }
                )
        }
    }
    [
        "contig 1 length 41 does not match 42", #
        "contig 3 missing",                     #
    ],
        "multiple warnings for multiple contig mismatches";
    like($exception, qr/^contigs do not match/, "detects overall contig mismatch");

    $exception = exception {
        $result = HMF::Pipeline::Functions::Validate::verifyReadGroups(
            [
                {tags => {ID => "id1"}},        #
                {tags => {ID => "id2"}},        #
            ], [
                {tags => {RG => {value => "id1"}}}, #
                {tags => {RG => {value => "id2"}}}, #
            ]
            )
    };
    is($exception, undef, "no exception for matching read groups");
    is($result, 1, "success for matching read groups");

    $exception = exception {
        $result = HMF::Pipeline::Functions::Validate::verifyReadGroups(
            [
                {tags => {ID => "id1"}},            #
                {tags => {ID => "id2"}},            #
                {tags => {ID => "id3"}},            #
            ], [
                {tags => {RG => {value => "id1"}}}, #
                {tags => {RG => {value => "id2"}}}, #
            ]
            )
    };
    is($exception, undef, "no exception for extra header read groups");
    is($result, 1, "success for extra header read groups");

    $exception = exception {
        $result = HMF::Pipeline::Functions::Validate::verifyReadGroups(
            [
                {tags => {ID => "id1"}},            #
            ], [
                {tags => {RG => {value => "id1"}}}, #
                {tags => {RG => {value => "id2"}}}, #
            ]
            )
    };
    like($exception, qr/^read group IDs from read tags not in BAM header:\n\tid2/, "detects missing header read groups");


    my $temp_file = File::Temp->new();
    $temp_file->DESTROY();

    $exception = exception { my $sample = HMF::Pipeline::Functions::Validate::verifyBam($temp_file->filename, {}) };
    like($exception, qr/^ERROR: $temp_file does not exist/, "detects missing BAM for verification");

    stderr_isnt {
        $exception = exception { $headers = HMF::Pipeline::Functions::Validate::bamHeaders($temp_file, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}}) }
    }
    "", "samtools shows errors for missing BAM for headers";
    like($exception, qr/^could not parse BAM headers from $temp_file/, "detects missing BAM for headers");

    stderr_isnt {
        $exception = exception { $reads = HMF::Pipeline::Functions::Validate::bamReads($temp_file, 1, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}}) }
    }
    "", "samtools shows errors for missing BAM for reads";
    like($exception, qr/^could not parse BAM reads from $temp_file/, "detects missing BAM for reads");

    stderr_isnt {
        $exception = exception { $reads = HMF::Pipeline::Functions::Validate::indexContigs($temp_file, {SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH}}) }
    }
    "", "samtools shows errors for missing BAM for index contigs";
    like($exception, qr/^could not read index stats from $temp_file/, "detects missing BAM for index contigs");

    stderr_isnt {
        $exception = exception { $reads = HMF::Pipeline::Functions::Validate::refGenomeContigs({REF_GENOME => $temp_file->filename}) }
    }
    "", "samtools shows errors for missing BAM for ref genome contigs";
    like($exception, qr/^could not read from ${temp_file}.fai/, "detects missing BAM for ref genome contigs");
}


done_testing();
