package illumina_poststats;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge;
use illumina_template;

sub runPostStats {
    my $configuration = shift;
    my %opt = %{$configuration};
    my @runningJobs;
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my $jobID = "PostStats_".getJobId();
    my $jobIDCheck = "PostStats_Check_".getJobId();

    if(! -e "$opt{OUTPUT_DIR}/logs/PostStats.done") {
		my $command = "perl $opt{BAMMETRICS_PATH}/bamMetrics.pl ";
		foreach my $sample (@{$opt{SAMPLES}}) {
			my $sampleBam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
			$command .= "-bam $sampleBam ";
			if (@{$opt{RUNNING_JOBS}->{$sample}}) {
				push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
			}
		}

		$command .= "-output_dir $opt{OUTPUT_DIR}/QCStats/ ";
		$command .= "-run_name $runName ";
		$command .= "-genome $opt{GENOME} ";
		$command .= "-queue $opt{POSTSTATS_QUEUE} ";
		$command .= "-queue_threads $opt{POSTSTATS_THREADS} ";
		$command .= "-queue_mem $opt{POSTSTATS_MEM} ";
		$command .= "-queue_time $opt{POSTSTATS_TIME} ";
		$command .= "-queue_project $opt{CLUSTER_PROJECT} ";
		$command .= "-picard_path $opt{PICARD_PATH} ";
		$command .= "-debug ";
		$command .= "-wgs ";
		$command .= "-coverage_cap 250 ";

		if ( $opt{SINGLE_END} ) {
			$command .= "-single_end ";
		}

		if ( $opt{CLUSTER_RESERVATION} eq "yes") {
			$command .= "-queue_reserve ";
		}

		my $bashFile = $opt{OUTPUT_DIR}."/jobs/".$jobID.".sh";
		my $logDir = $opt{OUTPUT_DIR}."/logs";

		from_template("PostStats.sh.tt", $bashFile, command => $command, runName => $runName, jobID => $jobID, jobIDCheck => $jobIDCheck,  opt => \%opt);

		my $qsub = &qsubTemplate(\%opt, "POSTSTATS");
		if (@runningJobs){
			system $qsub." -o ".$logDir."/PostStats_".$runName.".out -e ".$logDir."/PostStats_".$runName.".err -N ".$jobID." -hold_jid ".
			join(",",@runningJobs)." ".$bashFile;
		} else {
			system $qsub." -o ".$logDir."/PostStats_".$runName.".out -e ".$logDir."/PostStats_".$runName.".err -N ".$jobID." ".$bashFile;
		}

		my $bashFileCheck = $opt{OUTPUT_DIR}."/jobs/".$jobIDCheck.".sh";
		from_template("PostStatsCheck.sh.tt", $bashFileCheck, runName => $runName, opt => \%opt);

		system $qsub." -o ".$logDir."/PostStats_".$runName.".out -e ".$logDir."/PostStats_".$runName.".err -N ".$jobIDCheck.
			" -hold_jid bamMetrics_report_".$runName.",".$jobID." ".$bashFileCheck;
		return $jobIDCheck;
	} else {
		print "WARNING: $opt{OUTPUT_DIR}/logs/PostStats.done exists, skipping\n";
	}
}

1;
