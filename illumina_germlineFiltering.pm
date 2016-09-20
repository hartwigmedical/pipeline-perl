package illumina_germlineFiltering;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge;
use illumina_template;

sub runFilterVariants {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "GermlineFilter_".getJobId();

    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-e "$opt{OUTPUT_DIR}/logs/GermlineFilter.done" || -e "$opt{OUTPUT_DIR}/logs/VariantFilter.done"){
		print "WARNING: $opt{OUTPUT_DIR}/logs/GermlineFilter.done exists, skipping \n";
		return $jobID;
    }

    my $command = "java -Xmx".$opt{FILTER_MASTER_MEM}."G -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp -jar $opt{QUEUE_PATH}/Queue.jar ";
    my $jobNative = &jobNative(\%opt,"FILTER");
    $command .= "-jobQueue $opt{FILTER_QUEUE} -jobNative \"$jobNative\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/logs/GermlineFilter.jobReport.txt ";

    $command .= "-S $opt{OUTPUT_DIR}/QScripts/$opt{FILTER_SCALA} -R $opt{GENOME} -V $opt{OUTPUT_DIR}/$runName\.raw_variants.vcf -O $runName -mem $opt{FILTER_MEM} -nsc $opt{FILTER_SCATTER} ";

	my @SNPFilterNames = split("\t",$opt{FILTER_SNPNAME});
	my @SNPFilterExprs = split("\t",$opt{FILTER_SNPEXPR});
	my @snpTypes = split(",",$opt{FILTER_SNPTYPES});

	foreach my $snpType (@snpTypes) {
		$command.= "-snpType $snpType ";
	}

	if (scalar(@SNPFilterNames) ne scalar(@SNPFilterExprs)) {
		die "FILTER_SNPNAME and FILTER_SNPEXPR do not have the same length.";
	}

	foreach my $i (0 .. scalar(@SNPFilterNames)-1 ){
		$command .= "-snpFilterName $SNPFilterNames[$i] -snpFilterExpression \"$SNPFilterExprs[$i]\" ";
	}

	if ($opt{FILTER_CLUSTERSIZE} and $opt{FILTER_CLUSTERWINDOWSIZE}){
		$command .= "-cluster $opt{FILTER_CLUSTERSIZE} -window $opt{FILTER_CLUSTERWINDOWSIZE} ";
	}

	my @INDELFilterNames = split("\t",$opt{FILTER_INDELNAME});
	my @INDELFilterExprs = split("\t",$opt{FILTER_INDELEXPR});
	my @indelTypes = split(",",$opt{FILTER_INDELTYPES});

	foreach my $indelType (@indelTypes) {
		$command.= "-indelType $indelType ";
	}

	if (scalar(@INDELFilterNames) ne scalar(@INDELFilterExprs)) {
		die "FILTER_INDELNAME and FILTER_INDELEXPR do not have the same length.";
	}

	foreach my $i (0 .. scalar(@INDELFilterNames)-1 ) {
		$command .= "-indelFilterName $INDELFilterNames[$i] -indelFilterExpression \"$INDELFilterExprs[$i]\" ";
	}

    $command .= "-run";

    my $bashFile = $opt{OUTPUT_DIR}."/jobs/".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";
    from_template("GermlineFiltering.sh.tt", $bashFile, runName => $runName, command => $command, opt => \%opt);

    foreach my $sample (@{$opt{SAMPLES}}) {
        if (exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}}) {
            push( @runningJobs, join( ",", @{$opt{RUNNING_JOBS}->{$sample}} ) );
        }
    }

    my $qsub = &qsubJava( \%opt, "FILTER_MASTER" );
    if (@runningJobs) {
        system "$qsub -o $logDir/GermlineFilter_$runName.out -e $logDir/GermlineFilter_$runName.err -N $jobID -hold_jid ".join(
                ",", @runningJobs )." $bashFile";
    } else {
        system "$qsub -o $logDir/GermlineFilter_$runName.out -e $logDir/GermlineFilter_$runName.err -N $jobID $bashFile";
    }

    return $jobID;
}

1;
