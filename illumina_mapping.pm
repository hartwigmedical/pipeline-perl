#!/usr/bin/perl -w

##################################################################################################################################################
### illumina_mapping.pm
### - Map illumina sequencing data using bwa-mem
### - Use sambamba to merge lanes to a sample bam and mark duplicates.
### - Generate flagstats after each step to check bam integrity.
###
### Authors: S.W.Boymans, R.F.Ernst, H.H.D.kerstens
###
##################################################################################################################################################

package illumina_mapping;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_sge;
use illumina_template;

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

	if($opt{MARKDUP_LEVEL} eq "lane"){
	    print "Creating $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam with:\n";
	}elsif(($opt{MARKDUP_LEVEL} eq "no") || ($opt{MARKDUP_LEVEL} eq "sample")){
	    print "Creating $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam with:\n";
	}

	## Submit mapping jobs per fastq (pair)
	submitMappingJobs(\%opt,$QSUB, $samples, $sampleName, $coreName, $R1, $R2, $flowcellID);
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

	if(($opt{MARKDUP_LEVEL} eq "lane") || ($opt{MARKDUP_LEVEL} eq "sample")){
	    $opt{BAM_FILES}->{$sample} = "$sample\_dedup.bam";
	    print "Creating $opt{BAM_FILES}->{$sample}\n";
	}elsif($opt{MARKDUP_LEVEL} eq "no"){
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
	my $bams = join(" ", @bamList);

	### Create final merge script
	from_template("Merge.sh.tt", "$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh", sample => $sample, bamList => \@bamList, bams => $bams, opt => \%opt);

#	open MERGE_SH,">$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
#	print MERGE_SH "\#!/bin/sh\n\n";
#	print MERGE_SH "cd $opt{OUTPUT_DIR}/$sample\n";
#	print MERGE_SH "echo \"Start merge \t\" `date` \"\t$sample\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n\n";
#	print MERGE_SH "rm -f $opt{OUTPUT_DIR}/$sample/logs/$sample\_cleanup.err\n"; #rm old error file
#	print MERGE_SH "BAMS=( ".join(" ",@bamList)." )\n";
#	print MERGE_SH "PASS=0\n";
#	print MERGE_SH "for i in \"\$\{BAMS\[\@\]\}\"\n";
#	print MERGE_SH "do\n";
#
#	print MERGE_SH "\tDONEFILE=\`echo \$i | sed -r 's/\(_sorted)*(_dedup)*\\\.bam/\\\.done/'\`\n";
#
#
#	print MERGE_SH "\tif [ ! -f \$DONEFILE ]\n";
#	print MERGE_SH "\tthen\n";
#	print MERGE_SH "\t\techo \"ERROR: \$i is probably incomplete, no .done file found for it\" >> logs/merge.err\n";
#	print MERGE_SH "\t\tPASS=1\n";
#	print MERGE_SH "\tfi\n";
#	print MERGE_SH "done\n\n";
#	print MERGE_SH "if [ \$PASS -eq 1 ]\n";
#	print MERGE_SH "then\n";
#	print MERGE_SH "\techo \"ERROR: merging failed due to incomplete BAM-file(s)\" >> logs/merge.err\n";
#	print MERGE_SH "else\n";
#
#	if($opt{MARKDUP_LEVEL} eq "lane"){
#	    if(scalar(@bamList) > 1){
#		print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam ".join(" ",@bamList)."\n";
#	    } else {
#		print MERGE_SH "\tmv $bamList[0] mapping/$sample\_dedup.bam\n";
#	    }
#	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam mapping/$sample\_dedup.bai\n";
#	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam > mapping/$sample\_dedup.flagstat\n";
#
#	} elsif($opt{MARKDUP_LEVEL} eq "sample") {
#	    ### Use markdup to merge and markdup in one step, since sambamba v0.5.8
#	    print MERGE_SH "\techo \"Start markdup\t\" `date` \"\t$sample.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n";
#	    print MERGE_SH "ulimit -n 4096\n";
#	    ### Centos7 hpc: Use $TMPDIR as tmpdir variable.
#	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba markdup --tmpdir=$opt{OUTPUT_DIR}/$sample/tmp/ --overflow-list-size=$opt{MARKDUP_OVERFLOW_LIST_SIZE} -t $opt{MARKDUP_THREADS} ".join(" ",@bamList)." mapping/$sample\_dedup.bam\n";
#	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MARKDUP_THREADS} mapping/$sample\_dedup.bam mapping/$sample\_dedup.bai\n";
#	    ### compute resource efficient alternative
#	    #print MERGE_SH "module load sambamcram/biobambam\n"
#	    #print MERGE_SH "bammerge level=0 tmpfile=\$TMPDIR/",$sample,"merge I=",join(" I=",@bamList)," |bammarkduplicates2 tmpfile=\$TMPDIR/",$sample;
#	    #print MERGE_SH "mark markthreads=",$opt{MARKDUP_THREADS}," index=1 M=mapping/",$sample,"_metrics.txt O=mapping/",$sample,"_dedup.bam\n";
#	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MARKDUP_THREADS} mapping/$sample\_dedup.bam > mapping/$sample\_dedup.flagstat\n";
#	    print MERGE_SH "\techo \"End markdup\t\" `date` \"\t$sample.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n";
#
#	} elsif($opt{MARKDUP_LEVEL} eq "no") {
#	    if(scalar(@bamList) > 1){
#		print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{MARKDUP_THREADS} mapping/$sample.bam ".join(" ",@bamList)."\n";
#	    } else {
#		print MERGE_SH "\tmv $bamList[0] mapping/$sample.bam\n";
#	    }
#	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MARKDUP_THREADS} mapping/$sample.bam mapping/$sample.bai\n";
#	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MARKDUP_THREADS} mapping/$sample.bam > mapping/$sample.flagstat\n";
#	}
#
#	print MERGE_SH "fi\n\n";
#
#	print MERGE_SH "TOTALREADS=0\n";
#	if($opt{MARKDUP_LEVEL} eq "lane"){
#	    print MERGE_SH "for i in \$( find \$PWD/mapping -name '*sorted_dedup.flagstat')\n";
#	}elsif(($opt{MARKDUP_LEVEL} eq "no") || ($opt{MARKDUP_LEVEL} eq "sample")){
#	    print MERGE_SH "for i in \$( find \$PWD/mapping -name '*sorted.flagstat')\n";
#	}
#	print MERGE_SH "do\n";
#	print MERGE_SH "\tVAL=\`grep -m 1 -P \"\\d+\" \$i | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
#        print MERGE_SH "\tTOTALREADS=\$((\$TOTALREADS + \$VAL))\n";
#	print MERGE_SH "done\n\n";
#
#	if($opt{MARKDUP_LEVEL} eq "lane"){
#	    print MERGE_SH "if [ -s mapping/$sample\_dedup.flagstat ]\n";
#	    print MERGE_SH "then\n";
#	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
#	    print MERGE_SH "\tif [ \$FS1 -eq \$TOTALREADS ]\n";
#	    print MERGE_SH "\tthen\n";
#    	    print MERGE_SH "\t\tfor i in \$( find \$PWD/mapping -name '*sorted_dedup.bam')\n";
#    	    print MERGE_SH "\t\tdo\n";
#    	    print MERGE_SH "\t\t\trm \$i\n";
#    	    print MERGE_SH "\t\tdone\n";
#	    print MERGE_SH "\telse\n";
#    	    print MERGE_SH "\t\techo \"ERROR: read counts from *sorted_dedup.flagstat files and mapping/$sample\_dedup.flagstat do not match\" >> logs/$sample\_cleanup.err\n";
#    	    print MERGE_SH "\tfi\n";
#	    print MERGE_SH "else\n";
#	    print MERGE_SH "\techo \"ERROR: mapping/$sample\_dedup.flagstat is empty.\" >> logs/$sample\_cleanup.err\n";
#	    print MERGE_SH "fi\n\n";
#	}elsif($opt{MARKDUP_LEVEL} eq "sample"){
#	    print MERGE_SH "if [ -s mapping/$sample\_dedup.flagstat ]\n";
#	    print MERGE_SH "then\n";
#	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
#	    print MERGE_SH "\tif [ \$FS1 -eq \$TOTALREADS ]\n";
#	    print MERGE_SH "\tthen\n";
#	    print MERGE_SH "\t\trm mapping/$sample.bam 2> /dev/null\n";
#	    print MERGE_SH "\t\trm mapping/$sample.bai 2> /dev/null\n";
#    	    print MERGE_SH "\t\tfor i in \$( find \$PWD/mapping -name '*sorted.bam')\n";
#    	    print MERGE_SH "\t\tdo\n";
#    	    print MERGE_SH "\t\t\trm \$i\n";
#    	    print MERGE_SH "\t\tdone\n";
#	    print MERGE_SH "\telse\n";
#    	    print MERGE_SH "\t\techo \"ERROR: read counts from *sorted_dedup.flagstat files and mapping/$sample\_dedup.flagstat do not match\" >> logs/$sample\_cleanup.err\n";
#    	    print MERGE_SH "\tfi\n";
#	    print MERGE_SH "else\n";
#	    print MERGE_SH "\techo \"ERROR: mapping/$sample\_dedup.flagstat is empty.\" >> logs/$sample\_cleanup.err\n";
#	    print MERGE_SH "fi\n\n";
#	} elsif($opt{MARKDUP_LEVEL} eq "no"){
#	    print MERGE_SH "if [ -s mapping/$sample.flagstat ]\n";
#	    print MERGE_SH "then\n";
#	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" mapping/$sample.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
#	    print MERGE_SH "\tif [ \$FS1 -eq \$TOTALREADS ]\n";
#	    print MERGE_SH "\tthen\n";
#    	    print MERGE_SH "\t\tfor i in \$( find \$PWD/mapping -name '*sorted.bam')\n";
#    	    print MERGE_SH "\t\tdo\n";
#    	    print MERGE_SH "\t\t\trm \$i\n";
#    	    print MERGE_SH "\t\tdone\n";
#	    print MERGE_SH "\telse\n";
#    	    print MERGE_SH "\t\techo \"ERROR: read counts from *sorted.flagstat files and mapping/$sample.flagstat do not match\" >> logs/$sample\_cleanup.err\n";
#    	    print MERGE_SH "\tfi\n";
#	    print MERGE_SH "else\n";
#	    print MERGE_SH "\techo \"ERROR: mapping/$sample.flagstat is empty.\" >> logs/$sample\_cleanup.err\n";
#	    print MERGE_SH "fi\n\n";
#	}
#	print MERGE_SH "if [ ! -s logs/$sample\_cleanup.err ]\n";
#	print MERGE_SH "then\n";
#	print MERGE_SH "\ttouch logs/Mapping_$sample.done\n";
#	print MERGE_SH "fi\n\n";
#
#	print MERGE_SH "echo \"End merge \t\" `date` \"\t$sample\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sample/logs/$sample.log\n\n";
#	close MERGE_SH;
	my $qsub = &qsubTemplate(\%opt,"MARKDUP");
	print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sample,"/logs/Merge_",$sample,".out -e ",$opt{OUTPUT_DIR},"/",$sample,"/logs/Merge_",$sample,".err -N ",$jobId,
	" -hold_jid ",join(",",@jobIds)," ",$opt{OUTPUT_DIR},"/",$sample,"/jobs/",$jobId,".sh\n\n";

    }

    close $QSUB;

    system("sh $mainJobID");
    return \%opt;
}

sub submitMappingJobs{
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

    if($opt{MARKDUP_LEVEL} eq "lane"){
	push(@{$samples->{$sampleName}}, {'jobId'=>$cleanupJobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam"});
    }elsif(($opt{MARKDUP_LEVEL} eq "no") || ($opt{MARKDUP_LEVEL} eq "sample")){
	push(@{$samples->{$sampleName}}, {'jobId'=>$cleanupJobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam"});
    }

    ###Skip mapping if coreName_sorted_dedup.done file already exists
    if (-e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done"){
        print "\tWARNING: $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.done exists, skipping\n";
        return;
    }

    ### BWA JOB
    if (! -e "$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa.done"){
	print $R2 ? "\t$R1\n\t$R2\n" : "\t$R1\n";
	from_template("Map.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh", coreName => $coreName, sampleName => $sampleName, R1 => $R1, R2 => $R2, 
		RG_ID => $RG_ID, RG_SM => $RG_SM, RG_PL => $RG_PL, RG_LB => $RG_LB, RG_PU => $RG_PU, opt => \%opt);

	my $qsub = &qsubTemplate(\%opt,"MAPPING");
	print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
	$mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$mappingJobId,".sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa.done exists, skipping bwa\n";
    }
    ### BWA Flagstat
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat")){
	from_template("MapFS.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh", sampleName => $sampleName, coreName => $coreName, opt => \%opt);

	my $qsub = &qsubTemplate(\%opt,"FLAGSTAT");
	print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",
	$coreName,".err -N ",$mappingFSJobId," -hold_jid ",$mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$mappingFSJobId,".sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName.flagstat exist and is not empty, skipping bwa flagstat\n";
    }

    ### Sort bam
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam")) {
	from_template("Sort.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh", coreName => $coreName, sampleName => $sampleName, opt => \%opt);
	die "$opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh";
	#open SORT_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh\n";
	#print SORT_SH "\#!/bin/sh\n\n";
	#print SORT_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	#print SORT_SH "echo \"Start sort\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	#print SORT_SH "$opt{SAMBAMBA_PATH}/sambamba sort --tmpdir=$opt{OUTPUT_DIR}/$sampleName/tmp/ -m $opt{MAPPING_MEM}GB -t $opt{MAPPING_THREADS} -o $coreName\_sorted.bam $coreName.bam\n";
	#print SORT_SH "echo \"End sort\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	#close SORT_SH;

	my $qsub = &qsubTemplate(\%opt,"MAPPING");
	print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",$sortJobId,
	" -hold_jid ",$mappingJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$sortJobId,".sh\n";
    } else {
        print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bam exist and is not empty, skipping sort\n";
    }
    ### Sorted bam flagstat
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat")){
	from_template("SortFS.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh", sampleName => $sampleName, coreName => $coreName, opt => \%opt);

	#open FS2_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh\n";
	#print FS2_SH "\#!/bin/sh\n\n";
	#print FS2_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	#print FS2_SH "echo \"Start flagstat\t \" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	#print FS2_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{FLAGSTAT_THREADS} $coreName\_sorted.bam > $coreName\_sorted.flagstat\n";
	#print FS2_SH "echo \"End flagstat\t \" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n\n";
	#close FS2_SH;

	my $qsub = &qsubTemplate(\%opt,"FLAGSTAT");
	print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
	$sortFSJobId," -hold_jid ",$sortJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$sortFSJobId,".sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.flagstat exist and is not empty, skipping sorted bam flagstat\n";
    }
    ### Sorted bam index
    if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai")) {
	from_template("Index.sh.tt", "$opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh", sampleName => $sampleName, coreName => $coreName, opt => \%opt);
	my $qsub = &qsubTemplate(\%opt,"MAPPING");
	print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
	$indexJobId," -hold_jid ",$sortJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$indexJobId,".sh\n";
    } else {
	print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted.bai exist and is not empty, skipping sorted bam index\n";
    }

    ### Mark duplicates
    if($opt{MARKDUP_LEVEL} eq "lane"){
	if ((! -e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam") || (-z "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam")) {
	    open MARKDUP_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh\n";
	    print MARKDUP_SH "\#!/bin/sh\n\n";
	    print MARKDUP_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
	    print MARKDUP_SH "echo \"Start markdup\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    print MARKDUP_SH "$opt{SAMBAMBA_PATH}/sambamba markdup --overflow-list-size=$opt{MARKDUP_OVERFLOW_LIST_SIZE} --tmpdir=$opt{OUTPUT_DIR}/$sampleName/tmp/ -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted_dedup.bam \n";
	    print MARKDUP_SH "echo \"End markdup\t\" `date` \"\t$coreName\_sorted.bam\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    close MARKDUP_SH;
	    my $qsub = &qsubTemplate(\%opt,"MAPPING");
	    print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",
	    $markdupJobId," -hold_jid ",$indexJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$markdupJobId,".sh\n";
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

	    my $qsub = &qsubTemplate(\%opt,"FLAGSTAT");
	    print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR}/$sampleName,"/logs/Mapping_",$coreName,".err -N ",
	    $markdupFSJobId," -hold_jid ",$markdupJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$markdupFSJobId,".sh\n";
	} else {
	    print "\t$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.flagstat exist and is not empty, skipping dedup flagstat\n";
	}
    }

    ### Cleanup job, rm intermediate bams after checking flagstat
    my $cleanSh = "$opt{OUTPUT_DIR}/$sampleName/jobs/$cleanupJobId.sh";
    from_template("Clean.sh.tt", $cleanSh, sampleName => $sampleName, coreName => $coreName, opt => \%opt);

    my $qsub = &qsubTemplate(\%opt,"MAPPING");
    print $QSUB $qsub," -o ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".out -e ",$opt{OUTPUT_DIR},"/",$sampleName,"/logs/Mapping_",$coreName,".err -N ",$cleanupJobId,
    " -hold_jid ",$mappingFSJobId,",",$sortFSJobId,",",$markdupFSJobId," ",$opt{OUTPUT_DIR},"/",$sampleName,"/jobs/",$cleanupJobId,".sh\n\n";

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

	    my $qsub = &qsubTemplate(\%opt,"MAPPING");
	    system $qsub." -o ".$opt{OUTPUT_DIR}."/".$sample."/logs/PrepBam_".$sample.".out -e ".$opt{OUTPUT_DIR}."/".$sample."/logs/PrepBam_".$sample.".err -N ".$jobId
		." ".$opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobId.".sh";

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
