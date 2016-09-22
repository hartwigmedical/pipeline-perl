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
    my $configuration = shift;
    my %opt = %{$configuration};
    my $jobIds = {};

    my $mainJobID = "$opt{OUTPUT_DIR}/jobs/PreStatsMainJob_".getJobId().".sh";
    print "Creating FASTQC report for the following fastq.gz files:\n";

    my @qsubOut = ();
    my $runName = basename($opt{OUTPUT_DIR});
    foreach my $input (keys %{$opt{FASTQ}}) {
		my $coreName = undef;
		$coreName = fileparse($input);
		$coreName =~ s/\.fastq.gz$//;
		my ($sampleName) = split("_", $coreName);
		print "\t$input\n";

		if (!-e "$opt{OUTPUT_DIR}/$sampleName/logs/PreStats_$sampleName.done") {
			my $preStatsJobId = "PreStat_$coreName\_".getJobId();
			my $preStatsFile = "$opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh";
			push(@{$jobIds->{$sampleName}}, $preStatsJobId);
			from_template("PreStat.sh.tt", $preStatsFile, sampleName => $sampleName, coreName => $coreName, input => $input, opt => \%opt, runName => $runName);
			my $qsub = qsubTemplate(\%opt, "PRESTATS");
			push(@qsubOut, "$qsub -o $opt{OUTPUT_DIR}/$sampleName/logs/PreStat_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/PreStats_$coreName.err -N $preStatsJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh");
		} else {
			print "\t WARNING: FASTQC report for $input already exists, skipping.\n";
		}
    }

    from_template("PreStatsMainJob.sh.tt", $mainJobID, qsubOut => \@qsubOut, opt => \%opt);
    system("sh $mainJobID");
}

1;
