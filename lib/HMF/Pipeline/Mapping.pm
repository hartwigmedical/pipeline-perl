package HMF::Pipeline::Mapping;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Config::Validate qw(parseFastqName verifyBai verifyFlagstat);
use HMF::Pipeline::Sge qw(qsubSimple qsubTemplate);
use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Job::Bam qw(sorted indexed flagstat readCountCheck);

use parent qw(Exporter);
our @EXPORT_OK = qw(run runBamPrep);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING MAPPING ###";

    my $samples = {};
    foreach my $input_fastq (keys %{$opt->{FASTQ}}) {
        my $fastq = parseFastqName($input_fastq);
        say "Skipping R2 sample $input_fastq" and next if $input_fastq eq $fastq->{R2};

        if (exists $opt->{FASTQ}->{$fastq->{R2}}) {
            say "Switching to paired end mode!";
        } else {
            say "Switching to fragment mode!";
            $opt->{SINGLE_END} = 1;
        }

        my ($job_id, $sorted_bam, $sorted_flagstat, $dirs) = runPerLane($fastq, $opt);
        push @{$samples->{$fastq->{sampleName}}->{job_ids}}, $job_id;
        push @{$samples->{$fastq->{sampleName}}->{bams}}, $sorted_bam;
        push @{$samples->{$fastq->{sampleName}}->{flagstats}}, $sorted_flagstat;
        $samples->{$fastq->{sampleName}}->{dirs} = $dirs;
    }

    say "";
    while (my ($sample, $info) = each %{$samples}) {
        $opt->{BAM_FILES}->{$sample} = "${sample}_dedup.bam";
        say "Creating $opt->{BAM_FILES}->{$sample}";

        my $output_bam = catfile($info->{dirs}->{mapping}, "${sample}_dedup.bam");
        my $output_flagstat = catfile($info->{dirs}->{mapping}, "${sample}_dedup.flagstat");

        my $done_file = checkReportedDoneFile("Mapping", $sample, $info->{dirs}, $opt) or next;
        my $merge_job_id = fromTemplate(
            "MergeMarkdup",
            $sample,
            0,
            qsubTemplate($opt, "MARKDUP"),
            $info->{job_ids},
            $info->{dirs},
            $opt,
            sample => $sample,
            input_bams => $info->{bams},
            output_bam => $output_bam,
        );

        my $flagstat_job_id = flagstat($sample, $output_bam, $output_flagstat, [$merge_job_id], $info->{dirs}, $opt);
        my $check_job_id = readCountCheck(
            #<<< no perltidy
            $sample,
            $info->{flagstats},
            $output_flagstat,
            "MergeMarkdupSuccess.tt",
            { input_bams => $info->{bams} },
            [ $flagstat_job_id ],
            $info->{dirs},
            $opt,
            #>>> no perltidy
        );
        my $job_id = markDone($done_file, [ $merge_job_id, $flagstat_job_id, $check_job_id ], $info->{dirs}, $opt);
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }
    return;
}

sub runPerLane {
    my ($fastq, $opt) = @_;

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $fastq->{sampleName}), mapping => "mapping");
    my $core_file = catfile($dirs->{mapping}, $fastq->{coreName});
    my $unsorted_bam = "${core_file}.bam";
    my $sorted_bam = "${core_file}_sorted.bam";
    my $unsorted_flagstat = "${core_file}.flagstat";
    my $sorted_flagstat = "${core_file}_sorted.flagstat";

    my $done_file = checkReportedDoneFile("PerLaneConvert", $fastq->{coreName}, $dirs, $opt) or return (undef, $sorted_bam, $sorted_flagstat, $dirs);

    say "Creating ${sorted_bam} with:";
    say "\t$fastq->{R1}";
    say "\t$fastq->{R2}" if $fastq->{R2};

    my $map_job_id = fromTemplate(
        "Map",
        $fastq->{coreName},
        0,
        qsubTemplate($opt, "MAPPING"),
        [],
        $dirs,
        $opt,
        fastq => $fastq,
        output_bam => $unsorted_bam,
    );

    my $flagstat_job_id = flagstat($fastq->{coreName}, $unsorted_bam, $unsorted_flagstat, [$map_job_id], $dirs, $opt);
    my $sort_job_id = sorted($fastq->{coreName}, $unsorted_bam, $sorted_bam, [$flagstat_job_id], $dirs, $opt);
    my $sorted_flagstat_job_id = flagstat($fastq->{coreName}, $sorted_bam, $sorted_flagstat, [$sort_job_id], $dirs, $opt);
    my $check_job_id = readCountCheck(
        #<<< no perltidy
        $fastq->{coreName},
        [ $unsorted_flagstat ],
        $sorted_flagstat,
        "MapSortSuccess.tt",
        { unsorted_bam => $unsorted_bam },
        [ $flagstat_job_id, $sorted_flagstat_job_id ],
        $dirs,
        $opt,
        #>>> no perltidy
    );

    my $job_id = markDone($done_file, [ $map_job_id, $flagstat_job_id, $sort_job_id, $sorted_flagstat_job_id, $check_job_id ], $dirs, $opt);
    return ($job_id, $sorted_bam, $sorted_flagstat, $dirs);
}

sub runBamPrep {
    my ($opt) = @_;

    $opt->{MAPPING} = "no";
    $opt->{PRESTATS} = "no";
    say "\n### SCHEDULING BAM PREP ###";

    while (my ($sample, $input_bams) = each %{$opt->{SAMPLES}}) {
        my $input_bam = ${$input_bams}[0];
        my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $sample), mapping => "mapping");

        (my $input_bai = $input_bam) =~ s/\.bam$/.bam.bai/;
        (my $input_flagstat = $input_bam) =~ s/\.bam$/.flagstat/;

        my $bam_file = "${sample}.bam";
        $opt->{BAM_FILES}->{$sample} = $bam_file;
        my $sample_bam = catfile($dirs->{mapping}, $bam_file);
        my $sample_bai = catfile($dirs->{mapping}, "${sample}.bam.bai");
        my $sample_flagstat = catfile($dirs->{mapping}, "${sample}.flagstat");

        my $bai_good = verifyBai($input_bai, $input_bam, $opt);
        my $flagstat_good = verifyFlagstat($input_flagstat, $input_bam);

        symlink($input_bam, $sample_bam);
        $bai_good and symlink($input_bai, $sample_bai);
        $flagstat_good and symlink($input_flagstat, $sample_flagstat);

        next if $bai_good and $flagstat_good;

        my $index_job_id = indexed($sample, $sample_bam, $sample_bai, [], $dirs, $opt);
        my $flagstat_job_id = flagstat($sample, $sample_bam, $sample_flagstat, [$index_job_id], $dirs, $opt);
        my $check_job_id = fromTemplate(
            "ContigCheck",
            $sample,
            0,
            qsubSimple($opt),
            [$index_job_id],
            $dirs,
            $opt,
            step => $sample,
            bam_path => $sample_bam,
        );

        if (defined $index_job_id) {
            push @{$opt->{RUNNING_JOBS}->{$sample}}, $index_job_id;
        }

        if (defined $flagstat_job_id) {
            push @{$opt->{RUNNING_JOBS}->{$sample}}, $flagstat_job_id;
        }

        if (defined $check_job_id) {
            push @{$opt->{RUNNING_JOBS}->{$sample}}, $check_job_id;
        }
    }
    return;
}

1;
