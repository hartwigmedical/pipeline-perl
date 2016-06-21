#!/usr/bin/perl -w

package illumina_calling;

use strict;
use POSIX qw(tmpnam);
use illumina_sge;
use illumina_template;

sub runVariantCalling {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @sampleBams;
    my @runningJobs;
    my $jobID = "VC_".get_job_id();

    if (-e "$opt{OUTPUT_DIR}/logs/VariantCaller.done"){
	    print "WARNING: $opt{OUTPUT_DIR}/logs/VariantCaller.done exists, skipping \n";
	    return \%opt;
    }

    if((! -e "$opt{OUTPUT_DIR}/gvcf" && $opt{CALLING_GVCF} eq 'yes')){
	    mkdir("$opt{OUTPUT_DIR}/gvcf") or die "Couldn't create directory: $opt{OUTPUT_DIR}/gvcf\n";
    }

    my $jobNative = &jobNative(\%opt,"CALLING");
    my $command = "java -Xmx".$opt{CALLING_MASTER_MEM}."G -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp -jar $opt{QUEUE_PATH}/Queue.jar ";
    $command .= "-jobQueue $opt{CALLING_QUEUE} -jobNative \"$jobNative\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/logs/VariantCaller.jobReport.txt -memLimit $opt{CALLING_MEM} "; #Queue options

    $command .= "-S $opt{CALLING_SCALA} ";
    if ($opt{CALLING_UGMODE}) {
	    $command .= " -glm $opt{CALLING_UGMODE} ";
    }

    $command .= "-R $opt{GENOME} -O $runName -mem $opt{CALLING_MEM} -nct $opt{CALLING_THREADS} -nsc $opt{CALLING_SCATTER} -stand_call_conf $opt{CALLING_STANDCALLCONF} -stand_emit_conf $opt{CALLING_STANDEMITCONF} ";

    foreach my $sample (@{$opt{SAMPLES}}) {
        my $sampleBam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";

        $command .= "-I $sampleBam ";
        push( @sampleBams, $sampleBam);

        if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
            push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
        }
    }

    ### Optional settings
    if ( $opt{CALLING_DBSNP} ) {
        $command .= "-D $opt{CALLING_DBSNP} ";
    }

    if ( $opt{CALLING_TARGETS} ) {
        $command .= "-L $opt{CALLING_TARGETS} ";
        if ( $opt{CALLING_INTERVALPADDING} ) {
            $command .= "-ip $opt{CALLING_INTERVALPADDING} ";
        }
    }

    if ( $opt{CALLING_PLOIDY} ) {
        $command .= "-ploidy $opt{CALLING_PLOIDY} ";
    }

    $command .= "-run";

    my $bashFile = $opt{OUTPUT_DIR}."/jobs/VariantCalling_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";
    from_template("VariantCalling.sh.tt", $bashFile, runName => $runName, command => $command, sampleBams => \@sampleBams, opt => \%opt);

    my $qsub = &qsubJava(\%opt,"CALLING_MASTER");
    if (@runningJobs){
        system "$qsub -o $logDir/VariantCaller_$runName.out -e $logDir/VariantCaller_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
        system "$qsub -o $logDir/VariantCaller_$runName.out -e $logDir/VariantCaller_$runName.err -N $jobID $bashFile";
    }

    foreach my $sample (@{$opt{SAMPLES}}){
        push (@{$opt{RUNNING_JOBS}->{$sample}} , $jobID);
    }
    return \%opt;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
