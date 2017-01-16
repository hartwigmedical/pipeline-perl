#!/usr/bin/env perl

use strict;
use warnings;

use File::Temp;
use File::Touch;
use Test::Fatal;
use Test::More;

use HMF::Pipeline::Config::Validate;


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_KEY => \&HMF::Pipeline::Config::Validate::key_not_present,
        }, {
            REQUIRED_KEY => 1,
        },
    ),
    [],
    "required key"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_KEY => \&HMF::Pipeline::Config::Validate::key_not_present,
        },
        {},
    ),
    ["No REQUIRED_KEY option found in config files"],
    "detects missing key"
);


my $temp_file = File::Temp->new();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_FILE => \&HMF::Pipeline::Config::Validate::missing_file,
        }, {
            REQUIRED_FILE => $temp_file->filename,
        },
    ),
    [],
    "required file"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_FILE => \&HMF::Pipeline::Config::Validate::missing_file,
        },
        {},
    ),
    ["No REQUIRED_FILE option found in config files"],
    "detects missing file key"
);

$temp_file->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_FILE => \&HMF::Pipeline::Config::Validate::missing_file,
        }, {
            REQUIRED_FILE => $temp_file->filename,
        },
    ),
    ["REQUIRED_FILE file $temp_file does not exist"],
    "detects missing file"
);


$temp_file = File::Temp->new();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILE => \&HMF::Pipeline::Config::Validate::missing_optional_file,
        }, {
            OPTIONAL_FILE => $temp_file->filename,
        },
    ),
    [],
    "optional file"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILE => \&HMF::Pipeline::Config::Validate::missing_optional_file,
        },
        {},
    ),
    [],
    "optional file not given"
);

$temp_file->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILE => \&HMF::Pipeline::Config::Validate::missing_optional_file,
        }, {
            OPTIONAL_FILE => $temp_file->filename,
        },
    ),
    ["OPTIONAL_FILE file $temp_file does not exist"],
    "detects missing optional file"
);


my $temp_file_a = File::Temp->new();
my $temp_file_b = File::Temp->new();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILES => \&HMF::Pipeline::Config::Validate::missing_optional_files,
        }, {
            OPTIONAL_FILES => join("\t", ($temp_file_a->filename, $temp_file_b->filename)),
        },
    ),
    [],
    "optional files"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILES => \&HMF::Pipeline::Config::Validate::missing_optional_files,
        },
        {},
    ),
    [],
    "optional files not given"
);

$temp_file_a->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILES => \&HMF::Pipeline::Config::Validate::missing_optional_files,
        }, {
            OPTIONAL_FILES => join("\t", ($temp_file_a->filename, $temp_file_b->filename)),
        },
    ),
    ["OPTIONAL_FILES file $temp_file_a does not exist"],
    "detects one missing optional files"
);

$temp_file_b->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILES => \&HMF::Pipeline::Config::Validate::missing_optional_files,
        }, {
            OPTIONAL_FILES => join("\t", ($temp_file_a->filename, $temp_file_b->filename)),
        },
    ), [
        "OPTIONAL_FILES file $temp_file_a does not exist and
OPTIONAL_FILES file $temp_file_b does not exist"
    ],
    "detects all missing optional files"
);


$temp_file_a = File::Temp->new();
$temp_file_b = File::Temp->new();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            FILES => \&HMF::Pipeline::Config::Validate::missing_files,
        }, {
            FILES => join("\t", ($temp_file_a->filename, $temp_file_b->filename)),
        },
    ),
    [],
    "files"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            FILES => \&HMF::Pipeline::Config::Validate::missing_files,
        },
        {},
    ),
    ["No FILES option found in config files"],
    "files not given"
);

$temp_file_a->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            FILES => \&HMF::Pipeline::Config::Validate::missing_files,
        }, {
            FILES => join("\t", ($temp_file_a->filename, $temp_file_b->filename)),
        },
    ),
    ["FILES file $temp_file_a does not exist"],
    "detects one missing optional files"
);

$temp_file_b->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            FILES => \&HMF::Pipeline::Config::Validate::missing_files,
        }, {
            FILES => join("\t", ($temp_file_a->filename, $temp_file_b->filename)),
        },
    ), [
        "FILES file $temp_file_a does not exist and
FILES file $temp_file_b does not exist"
    ],
    "detects all missing optional files"
);


is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILES => HMF::Pipeline::Config::Validate::invalid_choice([ "CHOICE_A", "CHOICE_B" ]),
        }, {
            OPTIONAL_FILES => "CHOICE_B",
        },
    ),
    [],
    "choice"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            OPTIONAL_FILES => HMF::Pipeline::Config::Validate::invalid_choice([ "CHOICE_A", "CHOICE_B" ]),
        },
        {},
    ),
    [],
    "choice not given"
);

my $temp_dir = File::Temp->newdir();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_DIRECTORY => \&HMF::Pipeline::Config::Validate::missing_directory,
        }, {
            REQUIRED_DIRECTORY => $temp_dir->dirname,
        },
    ),
    [],
    "required directory"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_DIRECTORY => \&HMF::Pipeline::Config::Validate::missing_directory,
        },
        {},
    ),
    ["No REQUIRED_DIRECTORY option found in config files"],
    "detects missing directory key"
);

$temp_dir->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_DIRECTORY => \&HMF::Pipeline::Config::Validate::missing_directory,
        }, {
            REQUIRED_DIRECTORY => $temp_dir->dirname,
        },
    ),
    ["REQUIRED_DIRECTORY directory $temp_dir does not exist"],
    "detects missing directory"
);

$temp_file = File::Temp->new();
$temp_file_a = "${temp_file}.fai";
$temp_file_b = "${temp_file}.bwt";
touch($temp_file_a, $temp_file_b);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => \&HMF::Pipeline::Config::Validate::missing_genome_files,
        }, {
            KEY => $temp_file,
        },
    ),
    [],
    "finds genome files"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => \&HMF::Pipeline::Config::Validate::missing_genome_files,
        },
        {},
    ),
    ["No KEY option found in config files"],
    "missing genome key"
);

unlink($temp_file_b);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => \&HMF::Pipeline::Config::Validate::missing_genome_files,
        }, {
            KEY => $temp_file,
        },
    ),
    ["KEY file $temp_file_b does not exist"],
    "missing genome bwt"
);

touch($temp_file_b);
unlink($temp_file_a);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => \&HMF::Pipeline::Config::Validate::missing_genome_files,
        }, {
            KEY => $temp_file,
        },
    ),
    ["KEY file $temp_file_a does not exist"],
    "missing genome fai"
);

touch($temp_file_a);
$temp_file->DESTROY();

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => \&HMF::Pipeline::Config::Validate::missing_genome_files,
        }, {
            KEY => $temp_file,
        },
    ),
    ["KEY file $temp_file does not exist"],
    "missing genome file"
);

unlink($temp_file_a);
unlink($temp_file_b);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => \&HMF::Pipeline::Config::Validate::missing_genome_files,
        }, {
            KEY => $temp_file,
        },
    ),
    ["KEY file $temp_file does not exist and KEY file $temp_file_a does not exist and KEY file $temp_file_b does not exist"],
    "multiple missing genome files"
);


is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            CHOICE => HMF::Pipeline::Config::Validate::invalid_choice([ "CHOICE_A", "CHOICE_B" ]),
        }, {
            CHOICE => "INVALID_CHOICE",
        },
    ),
    ["CHOICE must be one of CHOICE_A, CHOICE_B"],
    "detects invalid choice"
);


is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::key_not_present_and_not_present("OTHER_KEY"),
        }, {
            KEY => 1,
            OTHER_KEY => 2,
        },
    ),
    [],
    "both keys present"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::key_not_present_and_not_present("OTHER_KEY"),
        }, {
            OTHER_KEY => 2,
        },
    ),
    [],
    "detects primary key missing"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::key_not_present_and_not_present("OTHER_KEY"),
        }, {
            KEY => 1,
        },
    ),
    [],
    "detects secondary key missing"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::key_not_present_and_not_present("OTHER_KEY"),
        },
        {},
    ),
    ["No KEY or OTHER_KEY found in config files"],
    "detects all key possibilities missing"
);


is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::compare_to("OTHER_KEY", sub { $_[0] eq $_[1] }, "equal to"),
        }, {
            KEY => 1,
            OTHER_KEY => 1,
        },
    ),
    [],
    "comparison"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::compare_to("OTHER_KEY", sub { $_[0] eq $_[1] }, "equal to"),
        }, {
            OTHER_KEY => 1,
        },
    ),
    ["No KEY option found in config files"],
    "detects comparison primary key missing"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::compare_to("OTHER_KEY", sub { $_[0] eq $_[1] }, "equal to"),
        }, {
            KEY => 1,
        },
    ),
    ["No OTHER_KEY option found in config files"],
    "detects comparison secondary key missing"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::compare_to("OTHER_KEY", sub { $_[0] eq $_[1] }, "equal to"),
        }, {
            KEY => 1,
            OTHER_KEY => 2,
        },
    ),
    ["KEY (1) must be equal to OTHER_KEY (2)"],
    "detects failed comparison"
);


is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::if_enabled({
                    OTHER_KEY => \&HMF::Pipeline::Config::Validate::key_not_present,
                }
            ),
        }, {
            KEY => "yes",
            OTHER_KEY => 1,
        },
    ),
    [],
    "other checks if enabled"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::if_enabled({
                    OTHER_KEY => \&HMF::Pipeline::Config::Validate::key_not_present,
                }
            ),
        }, {
            KEY => "yes",
        },
    ),
    ["No OTHER_KEY option found in config files"],
    "other checks detects missing key if enabled"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::if_enabled({
                    OTHER_KEY => \&HMF::Pipeline::Config::Validate::key_not_present,
                }
            ),
        }, {
            KEY => "no",
        },
    ),
    [],
    "no other checks if not enabled"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            KEY => HMF::Pipeline::Config::Validate::if_enabled({
                    OTHER_KEY => \&HMF::Pipeline::Config::Validate::key_not_present,
                }
            ),
        },
        {},
    ),
    ["No KEY option found in config files"],
    "detects missing key with other checks"
);


is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_KEY_A => \&HMF::Pipeline::Config::Validate::key_not_present,
            REQUIRED_KEY_B => \&HMF::Pipeline::Config::Validate::key_not_present,
        }, {
            REQUIRED_KEY_A => 1,
            REQUIRED_KEY_B => 1,
        },
    ),
    [],
    "multiple checks"
);

is_deeply(
    HMF::Pipeline::Config::Validate::applyChecks({
            REQUIRED_KEY_A => \&HMF::Pipeline::Config::Validate::key_not_present,
            REQUIRED_KEY_B => \&HMF::Pipeline::Config::Validate::key_not_present,
        }, {
            REQUIRED_KEY_A => 1,
        },
    ),
    ["No REQUIRED_KEY_B option found in config files"],
    "detects missing key after required key"
);

is_deeply(
    [
        sort(@{
                HMF::Pipeline::Config::Validate::applyChecks({
                        REQUIRED_KEY_A => \&HMF::Pipeline::Config::Validate::key_not_present,
                        REQUIRED_KEY_B => \&HMF::Pipeline::Config::Validate::key_not_present,
                    },
                    {},
                )
        })
    ], [
        sort((
                "No REQUIRED_KEY_A option found in config files", #
                "No REQUIRED_KEY_B option found in config files", #
        ))
    ],
    "multiple failures"
);

is_deeply(
    HMF::Pipeline::Config::Validate::parseFastqName("CPCT12345678R_HJJLGCCXX_S1_L001_R1_001.fastq.gz"), {
        R1 => "CPCT12345678R_HJJLGCCXX_S1_L001_R1_001.fastq.gz",
        R2 => "CPCT12345678R_HJJLGCCXX_S1_L001_R2_001.fastq.gz",
        coreName => "CPCT12345678R_HJJLGCCXX_S1_L001_001",
        sampleName => "CPCT12345678R",
        flowcellID => "HJJLGCCXX",
        index => "S1",
        lane => "L001",
        tag => "R1",
        suffix => "001",
        ext => ".fastq.gz",
    },
    "parses R1 FASTQ name"
);

is_deeply(
    HMF::Pipeline::Config::Validate::parseFastqName("CPCT12345678R_HJJLGCCXX_S1_L001_R2_001.fastq.gz"), {
        R1 => "CPCT12345678R_HJJLGCCXX_S1_L001_R1_001.fastq.gz",
        R2 => "CPCT12345678R_HJJLGCCXX_S1_L001_R2_001.fastq.gz",
        coreName => "CPCT12345678R_HJJLGCCXX_S1_L001_001",
        sampleName => "CPCT12345678R",
        flowcellID => "HJJLGCCXX",
        index => "S1",
        lane => "L001",
        tag => "R2",
        suffix => "001",
        ext => ".fastq.gz",
    },
    "parses R2 FASTQ name"
);

my $exception = exception { my $fastq = HMF::Pipeline::Config::Validate::parseFastqName("INVALID.fastq.gz") };
like($exception, qr/^ERROR: FASTQ filename 'INVALID.fastq.gz' must match regex/, "detects bad FASTQ name");

$exception = exception { my $fastq = HMF::Pipeline::Config::Validate::parseFastqName("CPCT12345678T_HJJLGCCXX_S3_L001_001_2.fastq") };
like($exception, qr/^ERROR: FASTQ filename 'CPCT12345678T_HJJLGCCXX_S3_L001_001_2.fastq' must match regex/, "detects bad FASTQ name");


done_testing();
