#!/usr/bin/perl -w

############################################################
### illumina_vcfutils.pm
### - Utility functions that can run after variant calling
###   - kinship analyses
###   - Phase by transmission
###
###Author: R.F.Ernst
###
############################################################

package illumina_vcfutils;

use strict;
use POSIX qw(tmpnam);

sub runVcfUtils {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "VCFUTILS_".get_job_id();

    if (-e "$opt{OUTPUT_DIR}/logs/VCF_UTILS.done"){
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
	    print VCFUTILS_SH "$opt{VCFTOOLS_PATH}/vcftools --vcf $opt{OUTPUT_DIR}/$vcf --plink\n";
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

    ### Phase by transmission
    if ( $opt{VCFUTILS_PHASE} eq "yes" ) {
	if (-e "$opt{OUTPUT_DIR}/logs/PhaseByTransmission.done"){
	    print "WARNING: $opt{OUTPUT_DIR}/logs/Phase.done exists, skipping \n";
	} else {
	    print VCFUTILS_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
	    print VCFUTILS_SH "java -Xmx8G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T PhaseByTransmission -R $opt{GENOME} -V $opt{OUTPUT_DIR}/$vcf -ped $opt{OUTPUT_DIR}/$runName.ped -o $runName.phased.vcf --MendelianViolationsFile $runName.MendelViol\n\n";
	    
	    ## Check output
	    print VCFUTILS_SH "if [ \"\$(tail -n 1 $opt{OUTPUT_DIR}/$vcf | cut -f 1,2)\" = \"\$(tail -n 1 $runName.phased.vcf | cut -f 1,2)\" -a -s $runName.MendelViol ]\n";
	    print VCFUTILS_SH "then\n";
	    print VCFUTILS_SH "\tmv $runName.phased.vcf $opt{OUTPUT_DIR}/\n";
	    print VCFUTILS_SH "\tmv $runName.MendelViol $opt{OUTPUT_DIR}/\n";
	    print VCFUTILS_SH "\ttouch $opt{OUTPUT_DIR}/logs/PhaseByTransmission.done\n";
	    print VCFUTILS_SH "else\n";
	    print VCFUTILS_SH "\tfailed=true\n";
	    print VCFUTILS_SH "fi\n\n";
	}
    }

    ### Gender check using plink
    if ( $opt{VCFUTILS_GENDERCHECK} eq "yes" ){
	if (-e "$opt{OUTPUT_DIR}/logs/Gender_check.done"){
	    print "WARNING: $opt{OUTPUT_DIR}/logs/Gender_check.done exists, skipping \n";
	} else {
	    print VCFUTILS_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
	    print VCFUTILS_SH "ln -sd $opt{OUTPUT_DIR}/$runName.ped $opt{OUTPUT_DIR}/$runName.fam\n";
	    print VCFUTILS_SH "java -Xmx8G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T VariantsToBinaryPed -R $opt{GENOME} -V $opt{OUTPUT_DIR}/$vcf -m $opt{OUTPUT_DIR}/$runName.fam -bed gender_check.bed -bim gender_check.bim -fam gender_check.fam -mgq 20\n";
	    print VCFUTILS_SH "$opt{PLINK_PATH}/plink -bfile gender_check --check-sex\n";
	    print VCFUTILS_SH "mv plink.sexcheck $opt{OUTPUT_DIR}/gender_check.out\n";
	    print VCFUTILS_SH "if [ -s $opt{OUTPUT_DIR}/gender_check.out ]; then\n";
	    print VCFUTILS_SH "\ttouch $opt{OUTPUT_DIR}/logs/Gender_check.done\n";
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
    if (@runningJobs){
	system "qsub -q $opt{VCFUTILS_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{VCFUTILS_THREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/VCFUTILS_$runName.out -e $logDir/VCFUTILS_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{VCFUTILS_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{VCFUTILS_THREADS} -P $opt{CLUSTER_PROJECT} -o $logDir/VCFUTILS_$runName.out -e $logDir/VCFUTILS_$runName.err -N $jobID $bashFile";
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