package Util;

use discipline;

use File::Spec::Functions;
use Test::Cmd;
use Test::More;

use parent qw(Exporter);
our @EXPORT_OK = qw(
    convertToBam
);


sub convertToBam {
    my ($sam_file, $with_indices) = @_;

    (my $bam_name = $sam_file) =~ s/\.sam$//;
    my $bam_file = "${bam_name}.bam";
    my $bai_file = "${bam_name}.bam.bai";
    my $flagstat_file = "${bam_name}.flagstat";

    unlink $bam_file, $bai_file, $flagstat_file;

    my $test_samtools = Test::Cmd->new(prog => catfile($ENV{SAMTOOLS_PATH}, "samtools"), workdir => "");
    $test_samtools->run(args => "view -bS -o ${bam_file} ${sam_file}");
    is($?, 0, "samtools successful (on ${sam_file})") or diag $test_samtools->stderr;

    if ($with_indices) {
        # required for stupid filesystem timestamp granularity
        sleep 1;

        $test_samtools->run(args => "index ${bam_file}");
        is($?, 0, "samtools index successful (on ${bam_file})") or diag $test_samtools->stderr;
        $test_samtools->run(args => "flagstat ${bam_file} > ${flagstat_file}");
        is($?, 0, "samtools flagstat successful (on ${bam_file})") or diag $test_samtools->stderr;
    }

    return $bam_file;
}

1;
