#!/usr/bin/perl -w

##################################################################################################################################################
### illumina_mapping.pm
### - Map illumina sequencing data using bwa-mem
### - Use sambamba to merge lanes to a sample bam and mark duplicates.
### - Generate flagstats after each step to check bam integrity.
###
### Authors: S.W.Boymans & R.F.Ernst
###
##################################################################################################################################################

package illumina_mapping;

use strict;
use POSIX qw(tmpnam);

sub runMapping {
    ###
    # Main mapping function, submits all separate mapping jobs.
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    
    my $FAI = "$opt{GENOME}\.fai";
    die "GENOME: $opt{GENOME} does not exists!\n" if !-e "$opt{GENOME}";
    die "GENOME BWT: $opt{GENOME}.bwt does not exists!\n" if !-e "$opt{GENOME}.bwt";
    die "GENOME FAI: $FAI does not exists!\n" if !-e $FAI;
    
    my $mainJobID = "$opt{OUTPUT_DIR}/jobs/MapMainJob_".get_job_id().".sh";

    open (my $QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print $QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";

    my $samples = {};
    my @jobs_to_wait;
    my $toMap = {};
    
    ### Try to search for matching pairs in the input FASTQ files
    foreach my $input (keys %{$opt{FASTQ}}){
	if($input =~ m/\_R1/){
	    my $pairName = $input;
	    $pairName =~ s/\_R1/\_R2/;
	    if(exists ($opt{FASTQ}->{$pairName})){
		$toMap->{$input."#".$pairName} = 1;
	    }else{
		$toMap->{$input} = 1;
	    }
	}elsif($input =~ m/\_R2/){
	    my $pairName = $input;
	    $pairName =~ s/\_R2/\_R1/;
	    if(exists ($opt{FASTQ}->{$pairName})){
		$toMap->{$pairName."#".$input} = 1;
	    }else{
		$toMap->{$input} = 1;
	    }
	}
    }

    foreach my $input (keys %{$toMap}){
	my @files = split("#",$input);
	my $R1 = undef;
        my $R2 = undef;
	my $coreName = undef;
    
	if(scalar(@files) == 2){
	    print "Switching to paired end mode!\n";
	    $R1 = $files[0];
	    $R2 = $files[1];
	    if($R1 !~ m/fastq.gz$/ or $R2 !~ m/fastq.gz$/){
		die "ERROR: Invalid input files:\n\t$R1\n\t$R2\n";
	    }
	}elsif(scalar(@files) == 1){
	    print "Switching to fragment mode!\n";
	    $R1 = $files[0];
	    $opt{SINGLE_END} = 1;
	    if($R1 !~ m/fastq.gz$/){
	        die "ERROR: Invalid input file:\n\t$R1\n";
	    }
	}else{
	    die "ERROR: Invalid input pair: $input\n";
	}

	$coreName = (split("/", $R1))[-1];
	my ($sampleName, $flowcellID, $index, $lane, $tag) =  split("_", $coreName);
	$coreName =~ s/\.fastq.gz//;
	$coreName =~ s/\_R1//;
	$coreName =~ s/\_R2//;
    
	my ($RG_PL, $RG_ID, $RG_LB, $RG_SM) = ('ILLUMINA', $coreName, $sampleName, $sampleName);
	
	if($opt{MAPPING_MARKDUP} eq "lane"){
	    print "Creating $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam with:\n";
	}elsif(($opt{MAPPING_MARKDUP} eq "no") || ($opt{MAPPING_MARKDUP} eq "sample")){
	    print "Creating $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam with:\n";
	}
	if($opt{MAPPING_MODE} eq 'batch'){
	    submitBatchJobs(\%opt,$QSUB,$samples, $sampleName, $coreName, $R1, $R2, $flowcellID);
	}elsif($opt{MAPPING_MODE} eq 'single'){
	    submitSingleJobs(\%opt,$QSUB, $samples, $sampleName, $coreName, $R1, $R2, $flowcellID);
	}
    }

    print "\n";
    ### Create a merging joblist for every sample
    foreach my $sample (keys %{$samples}){
	my @bamList = ();
	my @jobIds = ();
	my $pass = 1;

	foreach my $chunk (@{$samples->{$sample}}){
	    push(@bamList, $chunk->{'file'});
	    push(@jobIds, $chunk->{'jobId'});
	}

	if(($opt{MAPPING_MARKDUP} eq "lane") || ($opt{MAPPING_MARKDUP} eq "sample")){
	    $opt{BAM_FILES}->{$sample} = "$sample\_dedup.bam";
	    print "Creating $opt{BAM_FILES}->{$sample}\n";
	}elsif($opt{MAPPING_MARKDUP} eq "no"){
	    $opt{BAM_FILES}->{$sample} = "$sample.bam";
	    print "Creating $opt{BAM_FILES}->{$sample}\n";
	}

	### Skip mapping if dedup.done file already exists
	if (-e "$opt{OUTPUT_DIR}/$sample/logs/Mapping_$sample.done"){
	    print "\tWARNING: $opt{OUTPUT_DIR}/$sample/logs/Mapping_$sample.done exists, skipping\n";
	    next;
        }

	my $jobId = "Merge_$sample\_".get_job_id();
	push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobId);
	### Create final merge script
	open MERGE_SH,">$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	print MERGE_SH "\#!/bin/sh\n\n";
	print MERGE_SH "cd $opt{OUTPUT_DIR}/$sample\n";
	print MERGE_SH "echo \"Start merge \t\" `date` \"\t$sample\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n\n";
	print MERGE_SH "rm -f $opt{OUTPUT_DIR}/$sample/logs/$sample\_cleanup.err\n"; #rm old error file
	print MERGE_SH "BAMS=( ".join(" ",@bamList)." )\n";
	print MERGE_SH "PASS=0\n";
	print MERGE_SH "for i in \"\$\{BAMS\[\@\]\}\"\n";
	print MERGE_SH "do\n";
	print MERGE_SH "\tDONEFILE=\`echo \$i | sed -r 's/\(_sorted)*(_dedup)*\\\.bam/\\\.done/'\`\n";
	
	print MERGE_SH "\tif [ ! -f \$DONEFILE ]\n";
	print MERGE_SH "\tthen\n";
	print MERGE_SH "\t\techo \"ERROR: \$i is probably incomplete, no .done file found for it\" >> logs/merge.err\n";
	print MERGE_SH "\t\tPASS=1\n";
	print MERGE_SH "\tfi\n";
	print MERGE_SH "done\n\n";
	print MERGE_SH "if [ \$PASS -eq 1 ]\n";
	print MERGE_SH "then\n";
	print MERGE_SH "\techo \"ERROR: merging failed due to incomplete BAM-file(s)\" >> logs/merge.err\n";
	print MERGE_SH "else\n";
	
	if($opt{MAPPING_MARKDUP} eq "lane"){
	    if(scalar(@bamList) > 1){
		print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam ".join(" ",@bamList)."\n";
	    } else {
		print MERGE_SH "\tmv $bamList[0] mapping/$sample\_dedup.bam\n";
	    }
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam mapping/$sample\_dedup.bai\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam > mapping/$sample\_dedup.flagstat\n";
	
	} elsif($opt{MAPPING_MARKDUP} eq "sample") {
	    ### Use markdup to merge and markdup in one step, since sambamba v0.5.8
	    print MERGE_SH "\techo \"Start markdup\t\" `date` \"\t$sample.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba markdup --tmpdir=$opt{OUTPUT_DIR}/$sample/tmp/ --overflow-list-size=$opt{MAPPING_OVERFLOW_LIST_SIZE} -t $opt{MAPPING_THREADS} ".join(" ",@bamList)." mapping/$sample\_dedup.bam \n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam mapping/$sample\_dedup.bai\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam > mapping/$sample\_dedup.flagstat\n";
	    print MERGE_SH "\techo \"End markdup\t\" `date` \"\t$sample.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n";
	
	} elsif($opt{MAPPING_MARKDUP} eq "no") {
	    if(scalar(@bamList) > 1){
		print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{MAPPING_THREADS} mapping/$sample.bam ".join(" ",@bamList)."\n";
	    } else {
		print MERGE_SH "\tmv $bamList[0] mapping/$sample.bam\n";
	    }
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} mapping/$sample.bam mapping/$sample.bai\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} mapping/$sample.bam > mapping/$sample.flagstat\n";
	}
	
	print MERGE_SH "fi\n\n";
	
	print MERGE_SH "TOTALREADS=0\n";
	if($opt{MAPPING_MARKDUP} eq "lane"){
	    print MERGE_SH "for i in \$( find \$PWD/mapping -name '*sorted_dedup.flagstat')\n";
	}elsif(($opt{MAPPING_MARKDUP} eq "no") || ($opt{MAPPING_MARKDUP} eq "sample")){
	    print MERGE_SH "for i in \$( find \$PWD/mapping -name '*sorted.flagstat')\n";
	}
	print MERGE_SH "do\n";
	print MERGE_SH "\tVAL=\`grep -m 1 -P \"\\d+\" \$i | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
        print MERGE_SH "\tTOTALREADS=\$((\$TOTALREADS + \$VAL))\n";
	print MERGE_SH "done\n\n";
	
	if($opt{MAPPING_MARKDUP} eq "lane"){
	    print MERGE_SH "if [ -s mapping/$sample\_dedup.flagstat ]\n";
	    print MERGE_SH "then\n";
	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print MERGE_SH "\tif [ \$FS1 -eq \$TOTALREADS ]\n";
	    print MERGE_SH "\tthen\n";
    	    print MERGE_SH "\t\tfor i in \$( find \$PWD/mapping -name '*sorted_dedup.bam')\n";
    	    print MERGE_SH "\t\tdo\n";
    	    print MERGE_SH "\t\t\trm \$i\n";
    	    print MERGE_SH "\t\tdone\n";
	    print MERGE_SH "\telse\n";
    	    print MERGE_SH "\t\techo \"ERROR: read counts from *sorted_dedup.flagstat files and mapping/$sample\_dedup.flagstat do not match\" >> logs/$sample\_cleanup.err\n";
    	    print MERGE_SH "\tfi\n";
	    print MERGE_SH "else\n";
	    print MERGE_SH "\techo \"ERROR: mapping/$sample\_dedup.flagstat is empty.\" >> logs/$sample\_cleanup.err\n";
	    print MERGE_SH "fi\n\n";
	}elsif($opt{MAPPING_MARKDUP} eq "sample"){
	    print MERGE_SH "if [ -s mapping/$sample\_dedup.flagstat ]\n";
	    print MERGE_SH "then\n";
	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print MERGE_SH "\tif [ \$FS1 -eq \$TOTALREADS ]\n";
	    print MERGE_SH "\tthen\n";
	    print MERGE_SH "\t\trm mapping/$sample.bam\n";
	    print MERGE_SH "\t\trm mapping/$sample.bai\n";
    	    print MERGE_SH "\t\tfor i in \$( find \$PWD/mapping -name '*sorted.bam')\n";
    	    print MERGE_SH "\t\tdo\n";
    	    print MERGE_SH "\t\t\trm \$i\n";
    	    print MERGE_SH "\t\tdone\n";
	    print MERGE_SH "\telse\n";
    	    print MERGE_SH "\t\techo \"ERROR: read counts from *sorted_dedup.flagstat files and mapping/$sample\_dedup.flagstat do not match\" >> logs/$sample\_cleanup.err\n";
    	    print MERGE_SH "\tfi\n";
	    print MERGE_SH "else\n";
	    print MERGE_SH "\techo \"ERROR: mapping/$sample\_dedup.flagstat is empty.\" >> logs/$sample\_cleanup.err\n";
	    print MERGE_SH "fi\n\n";
	} elsif($opt{MAPPING_MARKDUP} eq "no"){
	    print MERGE_SH "if [ -s mapping/$sample.flagstat ]\n";
	    print MERGE_SH "then\n";
	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" mapping/$sample.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print MERGE_SH "\tif [ \$FS1 -eq \$TOTALREADS ]\n";
	    print MERGE_SH "\tthen\n";
    	    print MERGE_SH "\t\tfor i in \$( find \$PWD/mapping -name '*sorted.bam')\n";
    	    print MERGE_SH "\t\tdo\n";
    	    print MERGE_SH "\t\t\trm \$i\n";
    	    print MERGE_SH "\t\tdone\n";
	    print MERGE_SH "\telse\n";
    	    print MERGE_SH "\t\techo \"ERROR: read counts from *sorted.flagstat files and mapping/$sample.flagstat do not match\" >> logs/$sample\_cleanup.err\n";
    	    print MERGE_SH "\tfi\n";
	    print MERGE_SH "else\n";
	    print MERGE_SH "\techo \"ERROR: mapping/$sample.flagstat is empty.\" >> logs/$sample\_cleanup.err\n";
	    print MERGE_SH "fi\n\n";
	}
	print MERGE_SH "if [ ! -s logs/$sample\_cleanup.err ]\n";
	print MERGE_SH "then\n";
	print MERGE_SH "\ttouch logs/Mapping_$sample.done\n";
	print MERGE_SH "fi\n\n";
	
	print MERGE_SH "echo \"End merge \t\" `date` \"\t$sample\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n\n";
	close MERGE_SH;
    
	print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sample/logs/Merge_$sample.out -e $opt{OUTPUT_DIR}/$sample/logs/Merge_$sample.err -N $jobId -hold_jid ".join(",",@jobIds)." $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n\n";
    
    }

    close $QSUB;

    system("sh $mainJobID");
    return \%opt;
}

sub submitBatchJobs{
    ###
    # Submit batch jobs
    # One job per lane
    ###
    my ($opt,$QSUB ,$samples, $sampleName, $coreName, $R1, $R2, $flowcellID) = @_;
    my %opt = %$opt;
    my $jobId = "Map_$coreName\_".get_job_id();
    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM, $RG_PU) = ('ILLUMINA', $coreName, $sampleName, $sampleName, $flowcellID);
    
    if($opt{MAPPING_MARKDUP} eq "lane"){
	push(@{$samples->{$sampleName}}, {'jobId'=>$jobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam"});
    }elsif(($opt{MAPPING_MARKDUP} eq "no") || ($opt{MAPPING_MARKDUP} eq "sample")){
	push(@{$samples->{$sampleName}}, {'jobId'=>$jobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam"});
    }
    ### Skip mapping if coreName_sorted_dedup.done file already exists
    if (-e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done"){
        print "\tWARNING: $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done exists, skipping\n";
        return;
    }
    ### BWA mapping
    open BWA_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh\n";
    print BWA_SH "\#!/bin/sh\n\n";
    print BWA_SH "mkdir $opt{CLUSTER_TMP}/$jobId/\n";
    print BWA_SH "cd $opt{CLUSTER_TMP}/$jobId/ \n";
    print BWA_SH "echo \"Start mapping \t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n\n";
    
    if($R2){
        print "\t$R1\n\t$R2\n";
        print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} $opt{MAPPING_SETTINGS} -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU\" $opt{GENOME} $R1 $R2 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam.tmp /dev/stdin\n";
        print BWA_SH "mv $coreName.bam.tmp $coreName.bam\n";
    }else{
        print "\t$R1\n";
        print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} $opt{MAPPING_SETTINGS} -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU\" $opt{GENOME} $R1 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam.tmp /dev/stdin\n";
        print BWA_SH "mv $coreName.bam.tmp $coreName.bam\n";
    }
    print BWA_SH "echo \"End mapping\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n\n";
    print BWA_SH "echo \"Start flagstat\t\" `date` \"\t$coreName.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName.bam > $coreName.flagstat\n";
    print BWA_SH "echo \"End flagstat\t\" `date` \"\t$coreName.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    
    ### Sorting
    print BWA_SH "echo \"Start sort\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba sort --tmpdir=$opt{OUTPUT_DIR}/$sampleName/tmp/ -m $opt{MAPPING_MEM}GB -t $opt{MAPPING_THREADS} -o $coreName\_sorted.bam $coreName.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.err\n";
    print BWA_SH "echo \"End sort\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    print BWA_SH "echo \"Start flagstat\t \" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName\_sorted.bam > $coreName\_sorted.flagstat\n";
    print BWA_SH "echo \"End flagstat\t \" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    
    ### Index
    print BWA_SH "echo \"Start index\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted.bai 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.err\n";
    print BWA_SH "echo \"End index\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    
    ### MarkDup
    if($opt{MAPPING_MARKDUP} eq "lane"){
	print BWA_SH "echo \"Start markdup\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba markdup --overflow-list-size=$opt{MAPPING_OVERFLOW_LIST_SIZE} --tmpdir=$opt{OUTPUT_DIR}/$sampleName/tmp/ -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted_dedup.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.err\n";
	print BWA_SH "echo \"End markdup\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	print BWA_SH "echo \"Start flagstat\t\" `date` \"\t$coreName\_sorted_dedup.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName\_sorted_dedup.bam > $coreName\_sorted_dedup.flagstat\n";
	print BWA_SH "echo \"End flagstat\t\" `date` \"\t$coreName\_sorted_dedup.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    }
    ### Check and clean
    print BWA_SH "echo \"Start cleanup\t\" `date` \"\t $coreName \t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    print BWA_SH "rm -f $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n"; #rm old error file
    print BWA_SH "if [ -s $coreName.flagstat ] && [ -s $coreName\_sorted.flagstat ]\n";
    print BWA_SH "then\n";
    print BWA_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print BWA_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print BWA_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
    print BWA_SH "\tthen\n";
    print BWA_SH "\t\trm $coreName.bam\n";
    print BWA_SH "\telse\n";
    print BWA_SH "\t\techo \"ERROR: $coreName.flagstat and $coreName\_sorted.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print BWA_SH "\tfi\n";
    print BWA_SH "else\n";
    print BWA_SH "\techo \"ERROR: Either $coreName.flagstat or $coreName\_sorted.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print BWA_SH "fi\n\n";
    
    if($opt{MAPPING_MARKDUP} eq "lane"){
	print BWA_SH "if [ -s $coreName\_sorted.flagstat ] && [ -s $coreName\_sorted_dedup.flagstat ]\n";
	print BWA_SH "then\n";
	print BWA_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BWA_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BWA_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	print BWA_SH "\tthen\n";
	print BWA_SH "\t\trm $coreName\_sorted.bam\n";
	print BWA_SH "\t\trm $coreName\_sorted.bai\n";
	print BWA_SH "\telse\n";
	print BWA_SH "\t\techo \"ERROR: $coreName\_sorted.flagstat and $coreName\_sorted_dedup.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
	print BWA_SH "\tfi\n";
	print BWA_SH "else\n";
	print BWA_SH "\techo \"ERROR: Either $coreName\_sorted.flagstat or $coreName\_sorted_dedup.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
	print BWA_SH "fi\n\n";
    }
    
    print BWA_SH "if [ ! -s $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err ]\n";
    print BWA_SH "then\n";
    print BWA_SH "\ttouch $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done\n";
    print BWA_SH "fi\n\n";

    print BWA_SH "mv -f $opt{CLUSTER_TMP}/$jobId/* $opt{OUTPUT_DIR}/$sampleName/mapping/\n";
    if($opt{MAPPING_MARKDUP} eq "lane"){
	print BWA_SH "if [ -f $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam ];then\n";
    }elsif(($opt{MAPPING_MARKDUP} eq "no") || ($opt{MAPPING_MARKDUP} eq "sample")){
	print BWA_SH "if [ -f $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam ];then\n";
    }
    print BWA_SH "\trm -r $opt{CLUSTER_TMP}/$jobId/\n";
    print BWA_SH "fi\n";
    print BWA_SH "echo \"End cleanup\t\" `date` \"\t $coreName \t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";

    close BWA_SH;
	
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $jobId $opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh\n\n";

}

sub submitSingleJobs{
    ###
    # Submit single jobs
    # One job per lane per step
    ###
    my ($opt,$QSUB ,$samples, $sampleName, $coreName, $R1, $R2, $flowcellID) = @_;
    my %opt = %$opt;
    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM, $RG_PU) = ('ILLUMINA', $coreName, $sampleName, $sampleName, $flowcellID);
    
    my $mappingJobId = "Map_$coreName\_".get_job_id();
    my $mappingFSJobId = "MapFS_$coreName\_".get_job_id();
    my $sortJobId = "Sort_$coreName\_".get_job_id();
    my $sortFSJobId = "SortFS_$coreName\_".get_job_id();
    my $indexJobId = "Index_$coreName\_".get_job_id();
    my $markdupJobId = "MarkDup_$coreName\_".get_job_id();
    my $markdupFSJobId = "MarkDupFS_$coreName\_".get_job_id();
    my $cleanupJobId = "Clean_$coreName\_".get_job_id();
    
    if($opt{MAPPING_MARKDUP} eq "lane"){
	push(@{$samples->{$sampleName}}, {'jobId'=>$cleanupJobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam"});
    }elsif(($opt{MAPPING_MARKDUP} eq "no") || ($opt{MAPPING_MARKDUP} eq "sample")){
	push(@{$samples->{$sampleName}}, {'jobId'=>$cleanupJobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam"});
    }
    
    ###Skip mapping if coreName_sorted_dedup.done file already exists
    if (-e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done"){
        print "\tWARNING: $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done exists, skipping\n";
        return;
    }
    
    ### BWA JOB 
    if (! -e "$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa.done"){
	open BWA_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh\n";
	print BWA_SH "\#!/bin/sh\n\n";
	print BWA_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	print BWA_SH "echo \"Start mapping \t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n\n";

	if($R2){
	    print "\t$R1\n\t$R2\n";
	    print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} $opt{MAPPING_SETTINGS} -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU\" $opt{GENOME} $R1 $R2 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam.tmp /dev/stdin\n\n";
	}else{
	    print "\t$R1\n";
	    print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} $opt{MAPPING_SETTINGS} -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU\" $opt{GENOME} $R1 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam.tmp /dev/stdin\n\n";
	}
	### Check bwa completion
	print BWA_SH "if grep -Fq '[main] Real time:' $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log\n";
	print BWA_SH "then\n";
	print BWA_SH "\tmv $coreName.bam.tmp $coreName.bam\n";
	print BWA_SH "\ttouch $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa.done\n";
	print BWA_SH "fi\n";
	print BWA_SH "echo \"End mapping\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n\n";
	close BWA_SH;

	print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $mappingJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa.done exists, skipping bwa\n";
    }
    ### BWA Flagstat
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat")){
	open FS1_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh\n";
	print FS1_SH "\#!/bin/sh\n\n";
	print FS1_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	print FS1_SH "echo \"Start flagstat\t\" `date` \"\t$coreName.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	print FS1_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{FLAGSTAT_THREADS} $coreName.bam > $coreName.flagstat\n";
	print FS1_SH "echo \"End Flagstat \t\" `date` \"\t$coreName.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n\n";
	close FS1_SH;

	print $QSUB "qsub -q $opt{FLAGSTAT_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{FLAGSTAT_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $mappingFSJobId -hold_jid $mappingJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat exist and is not empty, skipping bwa flagstat\n";
    }

    ### Sort bam
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam")) {
	open SORT_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh\n";
	print SORT_SH "\#!/bin/sh\n\n";
	print SORT_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	print SORT_SH "echo \"Start sort\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	print SORT_SH "$opt{SAMBAMBA_PATH}/sambamba sort --tmpdir=$opt{OUTPUT_DIR}/$sampleName/tmp/ -m $opt{MAPPING_MEM}GB -t $opt{MAPPING_THREADS} -o $coreName\_sorted.bam $coreName.bam\n";
	print SORT_SH "echo \"End sort\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	close SORT_SH;

	print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $sortJobId -hold_jid $mappingJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh\n";
    } else {
        print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam exist and is not empty, skipping sort\n";
    }
    ### Sorted bam flagstat
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat")){
	open FS2_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh\n";
	print FS2_SH "\#!/bin/sh\n\n";
	print FS2_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	print FS2_SH "echo \"Start flagstat\t \" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	print FS2_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{FLAGSTAT_THREADS} $coreName\_sorted.bam > $coreName\_sorted.flagstat\n";
	print FS2_SH "echo \"End flagstat\t \" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n\n";
	close FS2_SH;

	print $QSUB "qsub -q $opt{FLAGSTAT_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{FLAGSTAT_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $sortFSJobId -hold_jid $sortJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat exist and is not empty, skipping sorted bam flagstat\n";
    }
    ### Sorted bam index
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai")) {
	open INDEX_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh\n";
	print INDEX_SH "\#!/bin/sh\n\n";
	print INDEX_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	print INDEX_SH "echo \"Start index\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	print INDEX_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted.bai \n";
	print INDEX_SH "echo \"End index\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	close INDEX_SH;

	print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $indexJobId -hold_jid $sortJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai exist and is not empty, skipping sorted bam index\n";
    }

    ### Mark duplicates
    if($opt{MAPPING_MARKDUP} eq "lane"){
	if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam")) {
	    open MARKDUP_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh\n";
	    print MARKDUP_SH "\#!/bin/sh\n\n";
	    print MARKDUP_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	    print MARKDUP_SH "echo \"Start markdup\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    print MARKDUP_SH "$opt{SAMBAMBA_PATH}/sambamba markdup --overflow-list-size=$opt{MAPPING_OVERFLOW_LIST_SIZE} --tmpdir=$opt{OUTPUT_DIR}/$sampleName/tmp/ -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted_dedup.bam \n";
	    print MARKDUP_SH "echo \"End markdup\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    close MARKDUP_SH;
	
	    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $markdupJobId -hold_jid $indexJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh\n";
	} else {
	    print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam exist and is not empty, skipping sort\n";
	}
    ### Mark duplicates flagstat
	if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.flagstat") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.flagstat")){
	    open FS3_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$markdupFSJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupFSJobId.sh\n";
	    print FS3_SH "\#!/bin/sh\n\n";
	    print FS3_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	    print FS3_SH "echo \"Start flagstat\t\" `date` \"\t$coreName\_sorted_dedup.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    print FS3_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{FLAGSTAT_THREADS} $coreName\_sorted_dedup.bam > $coreName\_sorted_dedup.flagstat\n";
	    print FS3_SH "echo \"End flagstat\t\" `date` \"\t$coreName\_sorted_dedup.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    close FS3_SH;
	    print $QSUB "qsub -q $opt{FLAGSTAT_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{FLAGSTAT_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $markdupFSJobId -hold_jid $markdupJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupFSJobId.sh\n";
	} else {
	    print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.flagstat exist and is not empty, skipping dedup flagstat\n";
	}
    }

    ### Cleanup job, rm intermediate bams after checking flagstat
    open CLEAN_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$cleanupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$cleanupJobId.sh\n";
    print CLEAN_SH "\#!/bin/sh\n\n";
    print CLEAN_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print CLEAN_SH "echo \"Start cleanup\t\" `date` \"\t $coreName \t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    print CLEAN_SH "rm -f $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n"; #rm old error file

    print CLEAN_SH "if [ -s $coreName.flagstat ] && [ -s $coreName\_sorted.flagstat ]\n";
    print CLEAN_SH "then\n";
    print CLEAN_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print CLEAN_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print CLEAN_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
    print CLEAN_SH "\tthen\n";
    print CLEAN_SH "\t\trm $coreName.bam\n";
    print CLEAN_SH "\telse\n";
    print CLEAN_SH "\t\techo \"ERROR: $coreName.flagstat and $coreName\_sorted.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print CLEAN_SH "\tfi\n";
    print CLEAN_SH "else\n";
    print CLEAN_SH "\techo \"ERROR: Either $coreName.flagstat or $coreName\_sorted.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print CLEAN_SH "fi\n\n";
    
    if($opt{MAPPING_MARKDUP} eq "lane"){
	print CLEAN_SH "if [ -s $coreName\_sorted.flagstat ] && [ -s $coreName\_sorted_dedup.flagstat ]\n";
	print CLEAN_SH "then\n";
	print CLEAN_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print CLEAN_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print CLEAN_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	print CLEAN_SH "\tthen\n";
	print CLEAN_SH "\t\trm $coreName\_sorted.bam\n";
	print CLEAN_SH "\t\trm $coreName\_sorted.bam.bai\n";
	print CLEAN_SH "\telse\n";
	print CLEAN_SH "\t\techo \"ERROR: $coreName\_sorted.flagstat and $coreName\_sorted_dedup.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
	print CLEAN_SH "\tfi\n";
	print CLEAN_SH "else\n";
	print CLEAN_SH "\techo \"ERROR: Either $coreName\_sorted.flagstat or $coreName\_sorted_dedup.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
	print CLEAN_SH "fi\n\n";
    }

    print CLEAN_SH "if [ ! -s $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err ]\n";
    print CLEAN_SH "then\n";
    print CLEAN_SH "\ttouch $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done\n";
    print CLEAN_SH "fi\n\n";
    print CLEAN_SH "echo \"End cleanup\t\" `date` \"\t $coreName \t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
    close CLEAN_SH;

    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/Mapping_$coreName.err -N $cleanupJobId -hold_jid $mappingFSJobId,$sortFSJobId,$markdupFSJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$cleanupJobId.sh\n\n";

}

sub runBamPrep {
    ###
    # Bam prep function, runs when starting the pipeline with bam files
    # Symlink or create index and flagstat files.
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $jobIds = {};
    foreach my $input (keys %{$opt{BAM}}){
	# Input Bam
	my $bamFile = (split("/", $input))[-1];
	my $input_bai = $input;
	my $input_flagstat = $input;
	$input_bai =~ s/\.bam/\.bai/g;
	$input_flagstat =~ s/\.bam/\.flagstat/g;
	# Ouput bam used in IAP
	my $sample = $bamFile;
	$sample =~ s/\.bam//g;
	$opt{BAM_FILES}->{$sample} = "$sample.bam";
	my $sample_bai = "$opt{OUTPUT_DIR}/$sample/mapping/$sample.bai";
	my $sample_flagstat = "$opt{OUTPUT_DIR}/$sample/mapping/$sample.flagstat";
	
	### Check input and symlink
	if (-e $input){
	    symlink($input,"$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}");
	} else {
	    die "ERROR: $input does not exist.";
	}
	if (-e $input_bai){
	    symlink($input_bai,"$sample_bai");
	}
	if (-e $input_flagstat){
	    symlink($input_flagstat,"$sample_flagstat");
	}
	
	### Index and Flagstat
	if (-e "$sample_bai" && -e "$sample_flagstat" ) {
	    next;
	} else {
	    my $jobId = "PrepBam_$sample\_".get_job_id();
	    open BAM_SH,">$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	    print BAM_SH "cd $opt{OUTPUT_DIR}/$sample/\n";
	    if (! -e "$sample_bai" ) { print BAM_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} mapping/$sample.bam $sample_bai\n"; }
	    if (! -e "$sample_flagstat" ) { print BAM_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} mapping/$sample.bam > $sample_flagstat\n"; }
	    close BAM_SH;
	
	    system "qsub -q $opt{MAPPING_QUEUE} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{MAPPING_THREADS} -R $opt{CLUSTER_RESERVATION} -o $opt{OUTPUT_DIR}/$sample/logs/PrepBam_$sample.out -e $opt{OUTPUT_DIR}/$sample/logs/PrepBam_$sample.err -N $jobId $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh";
	
	    push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobId);
	}
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