package HMF::Pipeline::Job;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;
use File::Temp qw(tmpnam);

use HMF::Pipeline::Template qw(writeFromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(
    getId
    fromTemplate
);


sub getId {
    my $name = tmpnam();
    my $id = fileparse($name);
    $id =~ s#(file|tmp\.[0-9]+\.)##;
    return $id;
}

sub fromTemplate {
    my ($name, $step, $qsub, $hold_jids, $dirs, $opt, %params) = @_;

    my $suffix = "";
    $suffix = "_${step}" if $step;

    my $job_id = "${name}${suffix}_" . getId();
    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    my $done_file = checkDone($dirs, $name, $suffix, $job_id);
    return unless $done_file;

    writeFromTemplate(
        "${name}.sh.tt",
        $bash_file,
        done_file => $done_file,
        dirs => $dirs,
        opt => $opt,
        %params,
    );

    my $stdout = catfile($dirs->{log}, "${name}${suffix}.out");
    my $stderr = catfile($dirs->{log}, "${name}${suffix}.err");
    my $hold_jid = "";
    $hold_jid = "-hold_jid " . join ",", @{$hold_jids} if @{$hold_jids};
    system "$qsub -o $stdout -e $stderr -N $job_id $hold_jid $bash_file";
    return $job_id;
}

sub checkDone {
    my ($dirs, $name, $suffix, $job_id) = @_;
    my $standard_done_name = "${name}${suffix}.done";
    my $standard_done_file = catfile($dirs->{log}, $standard_done_name);

    my $sample_name = (split /_/, $suffix)[1] // "";
    #<<< no perltidy
    my %old_done_files = (
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

    my $done_files = join ", ", grep { -f } ($standard_done_file, @{$old_done_files{$standard_done_name}});
    if ($done_files) {
        say "WARNING: $done_files exists, skipping $job_id";
        return;
    }
    return $standard_done_file;
}

1;
