package illumina_baf;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge;
use illumina_template;


sub runBAF {
    my ($opt) = @_;
    my $run_name = basename($opt->{OUTPUT_DIR});
    my @baf_jobs;

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
		my $out_dir = catfile($opt->{OUTPUT_DIR}, $sample);
        my $baf_dirs = {
            out => $out_dir,
            log => catfile($out_dir, "logs"),
            tmp => catfile($out_dir, "tmp"),
            job => catfile($out_dir, "jobs"),
        };
		my $sample_bam = catfile($baf_dirs->{out}, "mapping", $opt->{BAM_FILES}->{$sample});

        my $done_file = catfile($baf_dirs->{log}, "BAF_${sample}.done");
		if (-f $done_file) {
			say "WARNING: $done_file exists, skipping BAF analysis for $sample";
		} else {
			my $job_id = "BAF_$sample\_".getJobId();
			my $bash_file = catfile($baf_dirs->{job}, "${job_id}.sh");
			my $output_vcf = "${sample}_BAF_SNPS.vcf";
			my $output_baf = "${sample}_BAF.txt";
			my $output_bafplot = "${sample}_BAF.pdf";

            my @running_jobs;
			push @running_jobs, @{$opt->{RUNNING_JOBS}->{$sample}} if @{$opt->{RUNNING_JOBS}->{$sample}};

			my $run_unified_genotyper = 0;
            my $done_file = catfile($baf_dirs->{log}, "BAF_UG_${sample}.done");
			if (-f $done_file) {
				say "WARNING: $done_file exists, skipping Unified Genotyper for $sample";
			} else {
				$run_unified_genotyper = 1;
			}

			my $create_baf_file = 0;
            $done_file = catfile($baf_dirs->{log}, "BAF_FILE_${sample}.done");
			if (-f $done_file) {
				say "WARNING: $done_file exists, skipping BAF file for $sample";
			} else {
				$create_baf_file = 1;
			}
			my $create_baf_plots = 0;
            $done_file = catfile($baf_dirs->{log}, "BAF_PLOT_${sample}.done");
			if (-f $done_file) {
				say "WARNING: $done_file exists, skipping BAF plot for $sample";
			} else {
				$create_baf_plots = 1;
			}

			from_template("BAF.sh.tt", $bash_file,
                          sample => $sample,
                          sample_bam => $sample_bam,
                          output_vcf => $output_vcf,
                          output_baf => $output_baf,
                          output_bafplot => $output_bafplot,
                          run_unified_genotyper => $run_unified_genotyper,
                          create_baf_file => $create_baf_file,
                          create_baf_plots => $create_baf_plots,
                          run_name => $run_name,
                          dirs => $baf_dirs,
                          opt => $opt);

			my $qsub = qsubJava($opt, "BAF");
			if (@running_jobs) {
				system "$qsub -o $baf_dirs->{log}/BAF_${sample}.out -e $baf_dirs->{log}/BAF_${sample}.err -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
			} else {
				system "$qsub -o $baf_dirs->{log}/BAF_${sample}.out -e $baf_dirs->{log}/BAF_${sample}.err -N $job_id $bash_file";
			}
			push @baf_jobs, $job_id;
		}
    }
    $opt->{RUNNING_JOBS}->{'baf'} = \@baf_jobs;
}

1;
