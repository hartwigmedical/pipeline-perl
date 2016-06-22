#!/usr/bin/perl -w

package illumina_kinship;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin";
use illumina_sge;
use illumina_template;

sub runKinship {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "Kinship_".get_job_id();

    if (-e "$opt{OUTPUT_DIR}/logs/kinship.done") {
		print "WARNING: $opt{OUTPUT_DIR}/logs/kinship.done exists, skipping \n";
		return $jobID;
    }

    my $vcf = $runName.".filtered_variants.vcf";
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";

    open KINSHIP_SH, ">$bashFile" or die "cannot open file $bashFile \n";
    print KINSHIP_SH "#!/bin/bash\n\n";
    print KINSHIP_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
    print KINSHIP_SH "cd $opt{OUTPUT_DIR}/\n";
    print KINSHIP_SH "failed=false\n\n";
    print KINSHIP_SH "echo \"Start Kinship\t\" `date` \"\t$vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n\n";
    print KINSHIP_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
    print KINSHIP_SH "$opt{VCFTOOLS_PATH}/vcftools --temp $opt{OUTPUT_DIR}/tmp/ --vcf $opt{OUTPUT_DIR}/$vcf --plink\n";
    print KINSHIP_SH "$opt{PLINK_PATH}/plink --file out --make-bed --noweb\n";
    print KINSHIP_SH "$opt{KING_PATH}/king -b plink.bed --kinship\n";
    print KINSHIP_SH "cp king.kin0 $opt{OUTPUT_DIR}/$runName.kinship\n";
    print KINSHIP_SH "mv $opt{OUTPUT_DIR}/tmp/plink.log $opt{OUTPUT_DIR}/logs/\n";
    print KINSHIP_SH "mv $opt{OUTPUT_DIR}/tmp/out.log $opt{OUTPUT_DIR}/logs/\n";
    print KINSHIP_SH "if [ -s $opt{OUTPUT_DIR}/$runName.kinship ]; then\n";
    print KINSHIP_SH "\ttouch $opt{OUTPUT_DIR}/logs/Kinship.done\n";
    print KINSHIP_SH "else\n";
    print KINSHIP_SH "\tfailed=true\n";
    print KINSHIP_SH "fi\n\n";
    print KINSHIP_SH "if [ \"\$failed\" = false ]; then\n";
    print KINSHIP_SH "\ttouch $opt{OUTPUT_DIR}/logs/Kinship.done\n";
    print KINSHIP_SH "fi\n\n";
    print KINSHIP_SH "echo \"End Kinship\t\" `date` \"\t$vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";

    foreach my $sample (@{$opt{SAMPLES}}){
        if( exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}} ) {
            push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
        }
    }

    my $qsub = &qsubJava(\%opt, "KINSHIP");
    if (@runningJobs){
	    system "$qsub -o $logDir/Kinship_$runName.out -e $logDir/Kinship_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	    system "$qsub -o $logDir/Kinship_$runName.out -e $logDir/Kinship_$runName.err -N $jobID $bashFile";
    }

    return $jobID;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############
1;
