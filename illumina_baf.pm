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
    my $configuration = shift;
    my %opt = %{$configuration};
    my @baf_jobs;

    foreach my $sample (keys $opt{SAMPLES}) {
		my $sample_bam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
		my $log_dir = $opt{OUTPUT_DIR}."/".$sample."/logs/";
		my $tmp_dir = $opt{OUTPUT_DIR}."/".$sample."/tmp/";
		my $job_dir = $opt{OUTPUT_DIR}."/".$sample."/jobs/";
		my $output_dir = $opt{OUTPUT_DIR}."/".$sample."/";
		my $runName = basename($opt{OUTPUT_DIR});

		my @running_jobs;

		if (-e "$log_dir/BAF_$sample.done") {
			print "WARNING: $log_dir/BAF_$sample.done exists, skipping BAF analysis for $sample \n";
		} else {
			my $jobID = "BAF_$sample\_".getJobId();
			my $bashFile = $job_dir.$jobID.".sh";
			my $output_vcf = $sample."_BAF_SNPS.vcf";
			my $output_baf = $sample."_BAF.txt";
			my $output_bafplot = $sample."_BAF.pdf";

			if (@{$opt{RUNNING_JOBS}->{$sample}}) {
				push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample}});
			}

			my $run_unified_genotyper = 0;
			if (-e "$log_dir/BAF_UG_$sample.done") {
				print "WARNING: $log_dir/BAF_UG_$sample.done exists, skipping Unified Genotyper for $sample \n";
			} else {
				$run_unified_genotyper = 1;
			}

			my $create_baf_file = 0;
			if (-e "$log_dir/BAF_FILE_$sample.done") {
				print "WARNING: $log_dir/BAF_FILE_$sample.done exists, skipping BAF file for $sample \n";
			} else {
				$create_baf_file = 1;
			}
			my $create_baf_plots = 0;
			if (-e "$log_dir/BAF_PLOT_$sample.done") {
				print "WARNING: $log_dir/BAF_PLOT._$sample done exists, skipping BAF plot for $sample \n";
			} else {
				$create_baf_plots = 1;
			}

			from_template("BAF.sh.tt", $bashFile, tmp_dir => $tmp_dir, run_unified_genotyper => $run_unified_genotyper, log_dir => $log_dir, sample => $sample,
				sample_bam => $sample_bam, output_vcf => $output_vcf, output_dir => $output_dir, create_baf_file => $create_baf_file, output_baf => $output_baf,
				output_bafplot => $output_bafplot, create_baf_plots => $create_baf_plots, runName => $runName, opt => \%opt);

			my $qsub = &qsubJava(\%opt,"BAF");
			if (@running_jobs) {
				system "$qsub -o $log_dir/BAF_$sample.out -e $log_dir/BAF_$sample.err -N $jobID -hold_jid ".join(",",@running_jobs)." $bashFile";
			} else {
				system "$qsub -o $log_dir/BAF_$sample.out -e $log_dir/BAF_$sample.err -N $jobID $bashFile";
			}
			push(@baf_jobs, $jobID);
		}
    }
    return \@baf_jobs;
}

1;
