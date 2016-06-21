#!/usr/bin/perl -w

package illumina_vcfutils;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_sge;
use illumina_template;

sub runVcfUtils {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "VCFUTILS_".get_job_id();

    if (-e "$opt{OUTPUT_DIR}/logs/VCF_UTILS.done") {
		print "WARNING: $opt{OUTPUT_DIR}/logs/VCF_UTILS.done exists, skipping \n";
		return $jobID;
    }

    ### vcf file
    my $vcf;
    if($opt{FILTER_VARIANTS} eq "yes"){
	if($opt{FILTER_MODE} eq "BOTH"){ $vcf = $runName.".filtered_variants.vcf";}
	if($opt{FILTER_MODE} eq "SNP"){ $vcf = $runName.".filtered_snps.vcf";}
	if($opt{FILTER_MODE} eq "INDEL"){ $vcf = $runName.".filtered_indels.vcf";}
    } elsif ($opt{FILTER_VARIANTS} eq "no"){
	$vcf = $runName.".raw_variants.vcf";
    }

    ### Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";

    open VCFUTILS_SH, ">$bashFile" or die "cannot open file $bashFile \n";
    print VCFUTILS_SH "#!/bin/bash\n\n";
    print VCFUTILS_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
    print VCFUTILS_SH "cd $opt{OUTPUT_DIR}/\n";
    print VCFUTILS_SH "failed=false\n\n";
    print VCFUTILS_SH "echo \"Start VCF UTILS\t\" `date` \"\t$vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n\n";

    ### Run kinship analyses
    if ( $opt{VCFUTILS_KINSHIP} eq "yes" ) {
	if (-e "$opt{OUTPUT_DIR}/logs/Kinship.done"){
	    print "WARNING: $opt{OUTPUT_DIR}/logs/Kinship.done exists, skipping \n";
	} else {
	    print VCFUTILS_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
	    print VCFUTILS_SH "$opt{VCFTOOLS_PATH}/vcftools --temp $opt{OUTPUT_DIR}/tmp/ --vcf $opt{OUTPUT_DIR}/$vcf --plink\n";
	    print VCFUTILS_SH "$opt{PLINK_PATH}/plink --file out --make-bed --noweb\n";
	    print VCFUTILS_SH "$opt{KING_PATH}/king -b plink.bed --kinship\n";
	    print VCFUTILS_SH "cp king.kin0 $opt{OUTPUT_DIR}/$runName.kinship\n";
	    print VCFUTILS_SH "mv $opt{OUTPUT_DIR}/tmp/plink.log $opt{OUTPUT_DIR}/logs/\n";
	    print VCFUTILS_SH "mv $opt{OUTPUT_DIR}/tmp/out.log $opt{OUTPUT_DIR}/logs/\n";
	    print VCFUTILS_SH "if [ -s $opt{OUTPUT_DIR}/$runName.kinship ]; then\n";
	    print VCFUTILS_SH "\ttouch $opt{OUTPUT_DIR}/logs/Kinship.done\n";
	    print VCFUTILS_SH "else\n";
	    print VCFUTILS_SH "\tfailed=true\n";
	    print VCFUTILS_SH "fi\n\n";
	}
    }

    print VCFUTILS_SH "if [ \"\$failed\" = false ]; then\n";
    print VCFUTILS_SH "\ttouch $opt{OUTPUT_DIR}/logs/VCF_UTILS.done\n";
    print VCFUTILS_SH "fi\n\n";

    print VCFUTILS_SH "echo \"End VCF UTILS\t\" `date` \"\t$vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";

    ### Process runningjobs
    foreach my $sample (@{$opt{SAMPLES}}){
	if( exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}} ) {
	    push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	}
    }

    ### Run job
    my $qsub = &qsubJava(\%opt,"VCFUTILS");
    if (@runningJobs){
	system "$qsub -o $logDir/VCFUTILS_$runName.out -e $logDir/VCFUTILS_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "$qsub -o $logDir/VCFUTILS_$runName.out -e $logDir/VCFUTILS_$runName.err -N $jobID $bashFile";
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
