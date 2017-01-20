package HMF::Pipeline::Job;

use FindBin::libs;
use discipline;

use Carp;
use File::Basename;
use File::Spec::Functions;
use File::Temp qw(tmpnam);

use HMF::Pipeline::Template qw(writeFromTemplate);
use HMF::Pipeline::Sge qw(qsubSimple);

use parent qw(Exporter);
our @EXPORT_OK = qw(
    fromTemplate
    checkReportedDoneFile
    markDone
);


sub getId {
    my $name = tmpnam();
    my $id = fileparse($name);
    $id =~ s#(file|tmp\.[0-9]+\.)##;
    return $id;
}

sub fullname {
    my ($name, $step) = @_;
    return $step ? "${name}_${step}" : $name;
}

sub fromTemplate {
    my ($name, $step, $is_reported_job, $qsub, $hold_jids, $dirs, $opt, %params) = @_;

    my $job_name = fullname($name, $step);
    my $job_id = "${job_name}_" . getId();
    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    my $done_file = checkDoneFile($name, $step, $is_reported_job, $dirs, $opt);
    return unless $done_file;

    writeFromTemplate(
        "${name}.sh.tt",
        $bash_file,
        done_file => $done_file,
        dirs => $dirs,
        opt => $opt,
        %params,
    );

    submit($qsub, $job_name, $job_id, [ grep { defined } @{$hold_jids} ], $bash_file, $dirs);
    return $job_id;
}

sub submit {
    my ($qsub, $job_name, $job_id, $hold_jids, $bash_file, $dirs) = @_;

    my $stdout = catfile($dirs->{log}, "${job_name}.out");
    my $stderr = catfile($dirs->{log}, "${job_name}.err");
    my $hold_jid = hold_jid($hold_jids);

    my $option_file = File::Temp->new();
    open my $fh, ">", $option_file or confess "failed to write qsub parameter file $option_file";
    say $fh "-o $stdout -e $stderr -N $job_id $hold_jid";
    close $fh;
    system "$qsub -@ $option_file $bash_file";

    return;
}

sub hold_jid {
    my ($jids) = @_;
    return @{$jids} ? "-hold_jid " . join ",", @{$jids} : "";
}

# reported jobs are often the result of combining or post-processing
# steps which are: a) part of a multi-stage module b) checked before
# doing any work; c) run after the real work is done by unreported jobs
# (e.g. per-chromosome/lane -> combine). this function allows them to
# check a done file and get its name if it is not present. the name can
# then be used with markDone to quickly touch the (guaranteed same-name)
# done file after the combining job.
sub checkReportedDoneFile {
    my ($name, $step, $dirs, $opt) = @_;

    return checkDoneFile($name, $step, 1, $dirs, $opt);
}

# only used for multi-stage, reported jobs
# they want to report a specific name during finalize as described above
sub markDone {
    my ($done_file, $hold_job_ids, $dirs, $opt) = @_;

    return unless $done_file;
    my $done_file_name = fileparse($done_file);
    my $job_id = fromTemplate(
        "MarkDone",
        $done_file_name,
        0,
        qsubSimple($opt),
        $hold_job_ids,
        $dirs,
        $opt,
        done_file => $done_file,
        done_file_name => $done_file_name,
    );
    return $job_id;
}

sub checkDoneFile {
    my ($name, $step, $is_reported_job, $dirs, $opt) = @_;

    # need suffix to prevent hash clashes below, can use fullname() when removing backwards-compatibility
    my $suffix = "";
    $suffix = "_${step}" if $step;

    my $standard_done_name = "${name}${suffix}.done";
    my $standard_done_file = catfile($dirs->{log}, $standard_done_name);

    $step //= "";
    my ($sample_name) = split /_/, $step;
    $sample_name //= "";

    #<<< no perltidy
    my %old_done_files = (
        "Freec.done" => [
            catfile($dirs->{log}, "freec.done"),
        ],
        "QDNAseq.done" => [
            catfile($dirs->{log}, "qdnaseq.done"),
        ],
        "Strelka.done" => [
            catfile($dirs->{log}, "strelka.done"),
        ],
        "Varscan.done" => [
            catfile($dirs->{log}, "varscan.done"),
        ],
        "Freebayes.done" => [
            catfile($dirs->{log}, "freebayes.done"),
        ],
        "Mutect.done" => [
            catfile($dirs->{log}, "mutect.done"),
        ],
        "PerLaneConvert${suffix}.done" => [
            catfile($dirs->{mapping} // "", "${step}.done"),
        ],
        "Map${suffix}.done" => [
            catfile($dirs->{log}, "${step}_bwa.done"),
        ],
        "PreStats${suffix}.done" => [
            catfile($dirs->{log}, "PreStats_${sample_name}.done"),
        ],
        "GermlineCalling.done" => [
            catfile($dirs->{log}, "GermlineCaller.done"),
            catfile($dirs->{log}, "VariantCaller.done"),
        ],
        "GermlineFiltering.done" => [
            catfile($dirs->{log}, "GermlineFilter.done"),
            catfile($dirs->{log}, "VariantFilter.done"),
        ],
        "GermlineAnnotation.done" => [
            catfile($dirs->{log}, "VariantAnnotation.done"),
        ],
    );
    #>>> no perltidy

    my ($done_file) = grep { -f } ($standard_done_file, @{$old_done_files{$standard_done_name}});
    if ($done_file) {
        push @{$opt->{DONE_FILES}}, $done_file if $is_reported_job;
        say "WARNING: $done_file exists, skipping";
        return;
    }

    push @{$opt->{DONE_FILES}}, $standard_done_file if $is_reported_job;
    return $standard_done_file;
}

1;
