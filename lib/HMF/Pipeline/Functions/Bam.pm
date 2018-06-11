package HMF::Pipeline::Functions::Bam;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs);
use HMF::Pipeline::Functions::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate jobNative qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(
    sorted
    slice
    indexed
    flagstat
    readCountCheck
    diff
    prePostSliceAndDiff
    operationWithSliceChecks
    bamOperationWithSliceChecks
);

# SABR: "sort" would clash with Perl function
sub sorted {
    my ($step, $bam_path, $sorted_bam_path, $hold_jids, $dirs, $opt) = @_;

    my $sorted_bam_name = fileparse($sorted_bam_path);
    return fromTemplate(
        "BamSort",
        $sorted_bam_name,
        0,
        qsubTemplate($opt, "MAPPING"),
        $hold_jids,
        $dirs,
        $opt,
        step => $step,
        bam_path => $bam_path,
        sorted_bam_path => $sorted_bam_path,
        sorted_bam_name => $sorted_bam_name,
    );
}

sub slice {
    my ($step, $sample_bam, $sliced_bam, $bed_name, $hold_jids, $dirs, $opt) = @_;

    my $slice_name = fileparse($sliced_bam);
    return fromTemplate(
        "BamSlice",
        $slice_name,
        0,
        qsubTemplate($opt, "FLAGSTAT"),
        $hold_jids,
        $dirs,
        $opt,
        step => $step,
        input_bam => catfile($dirs->{mapping}, $sample_bam),
        bed_file => catfile($opt->{OUTPUT_DIR}, "settings", "slicing", $bed_name),
        sliced_bam => catfile($dirs->{mapping}, $sliced_bam),
        slice_name => $slice_name,
    );
}

# SABR: "index" would clash with Perl function
sub indexed {
    my ($step, $bam_path, $bai_path, $hold_jids, $dirs, $opt) = @_;

    my $bai_name = fileparse($bai_path);
    return fromTemplate(
        "BamIndex",
        $bai_name,
        0,
        qsubTemplate($opt, "MAPPING"),
        $hold_jids,
        $dirs,
        $opt,
        step => $step,
        bam_path => $bam_path,
        bai_path => $bai_path,
    );
}

sub flagstat {
    my ($step, $bam_path, $flagstat_path, $hold_jids, $dirs, $opt) = @_;

    my $flagstat_name = fileparse($flagstat_path);
    return fromTemplate(
        "BamFlagstat",
        $flagstat_name,
        0,
        qsubTemplate($opt, "FLAGSTAT"),
        $hold_jids,
        $dirs,
        $opt,
        step => $step,
        bam_path => $bam_path,
        flagstat_path => $flagstat_path,
        flagstat_name => $flagstat_name,
    );
}

sub readCountCheck {
    my ($step, $pre_flagstat_paths, $post_flagstat_path, $success_template, $extra_params, $hold_jids, $dirs, $opt) = @_;

    my $post_flagstat_name = fileparse($post_flagstat_path);
    return fromTemplate(
        "BamReadCountCheck",
        $post_flagstat_name,
        0,
        qsubTemplate($opt, "FLAGSTAT"),
        $hold_jids,
        $dirs,
        $opt,
        step => $step,
        pre_flagstat_paths => $pre_flagstat_paths,
        post_flagstat_path => $post_flagstat_path,
        post_flagstat_name => $post_flagstat_name,
        success_template => $success_template,
        %{$extra_params},
    );
}

sub diff {
    my ($step, $input_bam1, $input_bam2, $diff_name, $hold_jids, $dirs, $opt) = @_;

    return fromTemplate(
        "BamDiff",
        $diff_name,
        0,
        qsubTemplate($opt, "FLAGSTAT"),
        $hold_jids,
        $dirs,
        $opt,
        step => $step,
        diff_name => $diff_name,
        input_bam1 => catfile($dirs->{mapping}, $input_bam1),
        input_bam2 => catfile($dirs->{mapping}, $input_bam2),
        output_diff => catfile($dirs->{mapping}, $diff_name),
    );
}

sub prePostSliceAndDiff {
    my ($sample, $operation, $pre_bam, $post_bam, $hold_jids, $dirs, $opt) = @_;

    (my $pre_sliced_bam = $pre_bam) =~ s/\.bam$/.qc.pre${operation}.sliced.bam/;
    (my $post_sliced_bam = $pre_bam) =~ s/\.bam$/.qc.post${operation}.sliced.bam/;
    (my $post_sliced_flagstat = $pre_bam) =~ s/\.bam$/.qc.post${operation}.sliced.flagstat/;
    (my $pre_post_diff = $pre_bam) =~ s/\.bam$/.qc.prepost${operation}.diff/;

    my $post_sliced_bam_path = catfile($dirs->{mapping}, $post_sliced_bam);
    my $post_sliced_flagstat_path = catfile($dirs->{mapping}, $post_sliced_flagstat);

    my $pre_job_id = slice($sample, $pre_bam, $pre_sliced_bam, "QC_Slicing.bed", $hold_jids, $dirs, $opt);
    my $post_job_id = slice($sample, $post_bam, $post_sliced_bam, "QC_Slicing.bed", $hold_jids, $dirs, $opt);
    my $diff_job_id = diff($sample, $pre_sliced_bam, $post_sliced_bam, $pre_post_diff, [ $pre_job_id, $post_job_id ], $dirs, $opt);
    my $flagstat_job_id = flagstat($sample, $post_sliced_bam_path, $post_sliced_flagstat_path, [$post_job_id], $dirs, $opt);

    return [ $diff_job_id, $flagstat_job_id ];
}

sub operationWithSliceChecks {
    my ($job_template, $sample, $known_files, $post_tag, $slice_tag, $opt) = @_;

    my $sample_bam = $opt->{BAM_FILES}->{$sample};
    (my $post_bam = $sample_bam) =~ s/\.bam$/.${post_tag}.bam/;
    $opt->{BAM_FILES}->{$sample} = $post_bam;

    my (undef, $job_ids) = bamOperationWithSliceChecks($job_template, $sample, $sample_bam, $known_files, $post_tag, $slice_tag, $opt);

    push @{$opt->{RUNNING_JOBS}->{$sample}}, @{$job_ids};

    return;
}

sub bamOperationWithSliceChecks {
    my ($job_template, $sample, $sample_bam, $known_files, $post_tag, $slice_tag, $opt) = @_;

    (my $sample_flagstat = $sample_bam) =~ s/\.bam$/.flagstat/;
    (my $post_bai = $sample_bam) =~ s/\.bam$/.${post_tag}.bai/;
    (my $post_flagstat = $sample_bam) =~ s/\.bam$/.${post_tag}.flagstat/;
    (my $post_bam = $sample_bam) =~ s/\.bam$/.${post_tag}.bam/;

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $sample), mapping => "mapping");
    my $sample_bam_path = catfile($dirs->{mapping}, $sample_bam);
    (my $post_bam_path = $sample_bam_path) =~ s/\.bam$/.${post_tag}.bam/;
    say "\t${sample_bam_path}";

    my $done_file = checkReportedDoneFile($job_template, $sample, $dirs, $opt) or return ($post_bam_path, []);

    my $operation_job_id = fromTemplate(
        $job_template,
        $sample,
        0,
        qsubJava($opt, uc $job_template . "_MASTER"),
        $opt->{RUNNING_JOBS}->{$sample},
        $dirs,
        $opt,
        sample => $sample,
        sample_bam => $sample_bam,
        sample_bam_path => $sample_bam_path,
        job_native => jobNative($opt, uc $job_template),
        known_files => $known_files,
    ) or return ($post_bam_path, []);

    my $sample_flagstat_path = catfile($dirs->{mapping}, $sample_flagstat);
    my $post_flagstat_path = catfile($dirs->{mapping}, $post_flagstat);
    my $flagstat_job_id = flagstat($sample, catfile($dirs->{tmp}, $post_bam), $post_flagstat_path, [$operation_job_id], $dirs, $opt);
    my $check_job_id = readCountCheck(
        #<<< no perltidy
        $sample,
        [ $sample_flagstat_path ],
        $post_flagstat_path,
        "${job_template}Success.tt",
        { post_bam => $post_bam, post_bai => $post_bai, pre_bam => $sample_bam_path},
        [ $flagstat_job_id ],
        $dirs,
        $opt,
        #>>> no perltidy
    );

    my $job_id = markDone($done_file, [ $operation_job_id, $flagstat_job_id, $check_job_id ], $dirs, $opt);

    my $slicing_job_ids = prePostSliceAndDiff($sample, $slice_tag, $sample_bam, $post_bam, [$job_id], $dirs, $opt);
    push @{$opt->{RUNNING_JOBS}->{slicing}}, @{$slicing_job_ids};

    return ($post_bam_path, [$job_id]);
}

1;
