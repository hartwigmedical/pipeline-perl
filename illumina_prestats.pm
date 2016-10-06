package illumina_prestats;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge;
use illumina_template;


sub runPreStats {
    my ($opt) = @_;
    my $run_name = basename($opt->{OUTPUT_DIR});

    say "Creating FASTQC report for the following fastq.gz files:";
    my $job_ids = {};
    foreach my $input_file (keys %{$opt->{FASTQ}}) {
		my $core_name = fileparse($input_file);
		$core_name =~ s/\.fastq.gz$//;
		my ($sample_name) = split "_", $core_name;
		say "\t$input_file";

        my $sample_dir= catfile($opt->{OUTPUT_DIR}, $sample_name);
        my $log_dir = catfile($sample_dir, "logs");
        my $done_file = catfile($log_dir, "PreStats_${sample_name}.done");
		if (!-f $done_file) {
			my $job_id = "PreStats_${core_name}_" . getJobId();
			my $bash_file = catfile($sample_dir, "jobs", "${job_id}.sh");
			from_template("PreStats.sh.tt", $bash_file,
                          sample_name => $sample_name,
                          core_name => $core_name,
                          input => $input_file,
                          run_name => $run_name,
                          opt => $opt);
			my $qsub = qsubTemplate($opt, "PRESTATS");
			system("$qsub -o $log_dir/PreStat_${core_name}.out -e $log_dir/PreStats_${core_name}.err -N $job_id $bash_file");
			push @{$job_ids->{$sample_name}}, $job_id;
		} else {
			say "\tWARNING: $done_file exists, skipping";
		}
    }
}

1;
