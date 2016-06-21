#!/usr/bin/perl -w

package illumina_realign;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_sge;
use illumina_template;

sub runRealignment {

    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];

    print "Running single sample indel realignment for the following BAM-files:\n";
    
    my @knownIndelFiles;
    if($opt{REALIGNMENT_KNOWN}) {
		@knownIndelFiles = split('\t', $opt{REALIGNMENT_KNOWN});
    }

    foreach my $sample (@{$opt{SAMPLES}}){
	    my $bam = $opt{BAM_FILES}->{$sample};
	    (my $flagstat = $bam) =~ s/\.bam/.flagstat/;
	    (my $realignedBam = $bam) =~ s/\.bam/\.realigned\.bam/;
	    (my $realignedBai = $bam) =~ s/\.bam/\.realigned\.bai/;
	    (my $realignedBamBai = $bam) =~ s/\.bam/\.realigned\.bam\.bai/;
	    (my $realignedFlagstat = $bam) =~ s/\.bam/\.realigned\.flagstat/;

        (my $healthCheckPreRealignSlicedBam = $bam) =~ s/\.bam/\.qc.prerealign.sliced\.bam/;
        (my $healthCheckPreRealignSlicedBamBai = $bam) =~ s/\.bam/\.qc.prerealign.sliced\.bam\.bai/;
        (my $healthCheckPostRealignSlicedBam = $bam) =~ s/\.bam/\.qc.postrealign.sliced\.bam/;
        (my $healthCheckPostRealignSlicedBamBai = $bam) =~ s/\.bam/\.qc.postrealign.sliced\.bam\.bai/;
        (my $cpctSlicedBam = $bam) =~ s/\.bam/\.realigned.sliced\.bam/;
	    (my $cpctSlicedBamBai = $bam) =~ s/\.bam/\.realigned.sliced\.bam\.bai/;

        $opt{BAM_FILES}->{$sample} = $realignedBam;

	    print "\t$opt{OUTPUT_DIR}/$sample/mapping/$bam\n";

	    if (-e "$opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done"){
			print "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done exists, skipping\n";
			next;
	    }

	    my $logDir = $opt{OUTPUT_DIR}."/".$sample."/logs";
	    my $jobID = "Realign_".$sample."_".get_job_id();
	    my $bashFile = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobID.".sh";
	    my $jobNative = &jobNative(\%opt,"REALIGNMENT");

	    my $knownIndelFiles = "";
	    my @knownIndelFilesA = ();
	    if($opt{REALIGNMENT_KNOWN}) {
			foreach my $knownIndelFile (@knownIndelFiles) {
				if(! -e $knownIndelFile){ die"ERROR: $knownIndelFile does not exist\n" }
				else { push(@knownIndelFilesA, "-known $knownIndelFile"); }
			}
			$knownIndelFiles = join(" ", @knownIndelFilesA);
	    }

		from_template("Realign.sh.tt", $bashFile, sample => $sample, bam => $bam, logDir => $logDir, jobNative => $jobNative, knownIndelFiles => $knownIndelFiles,
            healthCheckPreRealignSlicedBam => $healthCheckPreRealignSlicedBam, healthCheckPreRealignSlicedBamBai => $healthCheckPreRealignSlicedBamBai, opt => \%opt, runName => $runName);

	    my $qsub = &qsubJava(\%opt,"REALIGNMENT_MASTER");
	    if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
			system $qsub." -o ".$logDir."/Realignment_".$sample.".out -e ".$logDir."/Realignment_".$sample.".err -N ".$jobID." -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." ".$bashFile;
	    } else {
			system $qsub." -o ".$logDir."/Realignment_".$sample.".out -e ".$logDir."/Realignment_".$sample.".err -N ".$jobID." ".$bashFile;
	    }

	    my $jobIDFS = "RealignFS_".$sample."_".get_job_id();
	    my $realignPostProcessScript = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobIDFS.".sh";

	    from_template("RealignPostProcess.sh.tt", $realignPostProcessScript, realignedBam => $realignedBam, realignedBai => $realignedBai, realignedBamBai => $realignedBamBai,
		    realignedFlagstat => $realignedFlagstat, flagstat => $flagstat, sample => $sample, logDir => $logDir, cpctSlicedBam => $cpctSlicedBam, cpctSlicedBamBai => $cpctSlicedBamBai,
            healthCheckPostRealignSlicedBam => $healthCheckPostRealignSlicedBam, healthCheckPostRealignSlicedBamBai => $healthCheckPostRealignSlicedBamBai, opt => \%opt, runName => $runName);

	    $qsub = &qsubTemplate(\%opt, "FLAGSTAT");
	    system $qsub." -o ".$logDir."/RealignmentPostProcess_".$sample.".out -e ".$logDir."/RealignmentPostProcess_".$sample.".err -N ".$jobIDFS." -hold_jid ".$jobID." ".$realignPostProcessScript;

	    push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobID);
	    push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobIDFS);
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
