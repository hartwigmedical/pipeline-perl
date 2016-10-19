package illumina_germlineCalling;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_sge qw(jobNative qsubJava);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);


sub runVariantCalling {
    my ($opt) = @_;

    my @sampleBams;
    my @runningJobs;
    my $jobID = "GermlineCalling_" . getJobId();

    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-f "$opt->{OUTPUT_DIR}/logs/GermlineCaller.done" || -f "$opt->{OUTPUT_DIR}/logs/VariantCaller.done") {
	    say "WARNING: $opt->{OUTPUT_DIR}/logs/GermlineCaller.done exists, skipping";
	    return;
    }

    if ((!-d "$opt->{OUTPUT_DIR}/gvcf" && $opt->{CALLING_GVCF} eq 'yes')) {
	    mkdir("$opt->{OUTPUT_DIR}/gvcf") or die "Couldn't create directory $opt->{OUTPUT_DIR}/gvcf: $!";
    }

    my $jobNative = jobNative($opt, "CALLING");
    my $command = "java -Xmx".$opt->{CALLING_MASTER_MEM}."G -Djava.io.tmpdir=$opt->{OUTPUT_DIR}/tmp -jar $opt->{QUEUE_PATH}/Queue.jar ";
    $command .= "-jobQueue $opt->{CALLING_QUEUE} -jobNative \"$jobNative\" -jobRunner GridEngine -jobReport $opt->{OUTPUT_DIR}/logs/GermlineCaller.jobReport.txt -memLimit $opt->{CALLING_MEM} ";

    $command .= "-S $opt->{OUTPUT_DIR}/QScripts/$opt->{CALLING_SCALA} ";
    if ($opt->{CALLING_UGMODE}) {
	    $command .= " -glm $opt->{CALLING_UGMODE} ";
    }

    $command .= "-R $opt->{GENOME} -O $opt->{RUN_NAME} -mem $opt->{CALLING_MEM} -nct $opt->{CALLING_THREADS} -nsc $opt->{CALLING_SCATTER} -stand_call_conf $opt->{CALLING_STANDCALLCONF} -stand_emit_conf $opt->{CALLING_STANDEMITCONF} ";

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $sampleBam = "$opt->{OUTPUT_DIR}/$sample/mapping/$opt->{BAM_FILES}->{$sample}";

        $command .= "-I $sampleBam ";
        push(@sampleBams, $sampleBam);

        if (@{$opt->{RUNNING_JOBS}->{$sample}}) {
            push(@runningJobs, @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    if ($opt->{CALLING_DBSNP}) {
        $command .= "-D $opt->{CALLING_DBSNP} ";
    }

    if ($opt->{CALLING_TARGETS}) {
        $command .= "-L $opt->{CALLING_TARGETS} ";
        if ($opt->{CALLING_INTERVALPADDING}) {
            $command .= "-ip $opt->{CALLING_INTERVALPADDING} ";
        }
    }

    if ($opt->{CALLING_PLOIDY}) {
        $command .= "-ploidy $opt->{CALLING_PLOIDY} ";
    }

    $command .= "-run";

    my $bashFile = catfile($opt->{OUTPUT_DIR}, "jobs", "${jobID}.sh");
    my $logDir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $stdout = catfile($logDir, "GermlineCaller_$opt->{RUN_NAME}.out");
    my $stderr = catfile($logDir, "GermlineCaller_$opt->{RUN_NAME}.err");

    from_template("GermlineCalling.sh.tt", $bashFile,
                  command => $command,
                  sampleBams => \@sampleBams,
                  opt => $opt);

    my $qsub = qsubJava($opt, "CALLING_MASTER");
    if (@runningJobs) {
        system "$qsub -o $stdout -e $stderr -N $jobID -hold_jid " . join(",", @runningJobs) . " $bashFile";
    } else {
        system "$qsub -o $stdout -e $stderr -N $jobID $bashFile";
    }

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $jobID;
    }
    return;
}

1;
