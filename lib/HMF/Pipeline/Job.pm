package HMF::Pipeline::Job;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;
use POSIX qw(tmpnam);

use HMF::Pipeline::Template qw(writeFromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(
    getId
    fromTemplate
);


sub getId {
    my $id = fileparse(tmpnam());
    $id =~ s#(file|tmp\.[0-9]+\.)##;
    return $id;
}

sub fromTemplate {
    my ($name, $sample, $qsub, $hold_jids, $dirs, $opt, %params) = @_;

    my $suffix = "";
    $suffix = "_${sample}" if $sample;

    my $done_file = catfile($dirs->{log}, "${name}${suffix}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my $job_id = "${name}${suffix}_" . getId();
    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "${name}.sh.tt", $bash_file,
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

1;
