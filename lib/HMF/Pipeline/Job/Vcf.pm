package HMF::Pipeline::Job::Vcf;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Sge qw(qsubTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(
    sorted
    concat
);


sub sorted {
    my ($input_vcf, $output_vcf, $step, $qsub, $hold_job_ids, $dirs, $opt) = @_;

    my $vcf_name = fileparse($output_vcf);
    my $job_id = fromTemplate(
        "SortVcf",
        $vcf_name,
        qsubTemplate($opt, $qsub),
        $hold_job_ids,
        $dirs,
        $opt,
        input_vcf => $input_vcf,
        output_vcf => $output_vcf,
        step => $step,
    );
    return $job_id;
}

sub concat {
    my ($vcf_files, $output_vcf, $step, $qsub, $hold_job_ids, $dirs, $opt) = @_;

    my $vcf_name = fileparse($output_vcf);
    my $job_id = fromTemplate(
        "ConcatVcf",
        $vcf_name,
        qsubTemplate($opt, $qsub),
        $hold_job_ids,
        $dirs,
        $opt,
        vcf_files => $vcf_files,
        output_vcf => $output_vcf,
        step => $step,
    );
    return $job_id;
}

1;
