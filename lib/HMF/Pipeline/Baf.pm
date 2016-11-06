package HMF::Pipeline::Baf;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;
    my @baf_jobs;

    say "\n### SCHEDULING BAF ANALYSIS ###";

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $out_dir = catfile($opt->{OUTPUT_DIR}, $sample);
        my $dirs = createDirs($out_dir);
        my $sample_bam = catfile($dirs->{out}, "mapping", $opt->{BAM_FILES}->{$sample});

        my $job_id = "BAF_$sample\_" . getId();
        my $done_file = catfile($dirs->{log}, "BAF_${sample}.done");
        if (-f $done_file) {
            say "WARNING: $done_file exists, skipping $job_id";
        } else {
            my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
            my $output_vcf = "${sample}_BAF_SNPS.vcf";
            my $output_baf = "${sample}_BAF.txt";
            my $output_bafplot = "${sample}_BAF.pdf";

            my @running_jobs;
            push @running_jobs, @{$opt->{RUNNING_JOBS}->{$sample}} if @{$opt->{RUNNING_JOBS}->{$sample}};

            my $run_unified_genotyper = 0;
            $done_file = catfile($dirs->{log}, "BAF_UG_${sample}.done");
            if (-f $done_file) {
                say "WARNING: $done_file exists, skipping Unified Genotyper for $sample";
            } else {
                $run_unified_genotyper = 1;
            }

            my $create_baf_file = 0;
            $done_file = catfile($dirs->{log}, "BAF_FILE_${sample}.done");
            if (-f $done_file) {
                say "WARNING: $done_file exists, skipping BAF file for $sample";
            } else {
                $create_baf_file = 1;
            }
            my $create_baf_plots = 0;
            $done_file = catfile($dirs->{log}, "BAF_PLOT_${sample}.done");
            if (-f $done_file) {
                say "WARNING: $done_file exists, skipping BAF plot for $sample";
            } else {
                $create_baf_plots = 1;
            }

            writeFromTemplate(
                "BAF.sh.tt", $bash_file,
                sample => $sample,
                sample_bam => $sample_bam,
                output_vcf => $output_vcf,
                output_baf => $output_baf,
                output_bafplot => $output_bafplot,
                run_unified_genotyper => $run_unified_genotyper,
                create_baf_file => $create_baf_file,
                create_baf_plots => $create_baf_plots,
                dirs => $dirs,
                opt => $opt,
            );

            my $qsub = qsubJava($opt, "BAF");
            if (@running_jobs) {
                system "$qsub -o $dirs->{log}/BAF_${sample}.out -e $dirs->{log}/BAF_${sample}.err -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
            } else {
                system "$qsub -o $dirs->{log}/BAF_${sample}.out -e $dirs->{log}/BAF_${sample}.err -N $job_id $bash_file";
            }
            push @baf_jobs, $job_id;
        }
    }
    $opt->{RUNNING_JOBS}->{'baf'} = \@baf_jobs;
    return;
}

1;
