#!/usr/bin/perl -w

##################################################################
### illumina_structuralVariants.pm
### - Run structural variant caller Delly
###
### Author: R.F.Ernst & M. van Roosmalen
##################################################################

package illumina_structuralVariants;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);

my @runningJobs;
my @sampleBams;
my ($delly_out_dir, $delly_log_dir, $delly_job_dir, $delly_tmp_dir);
my %opt;

sub runDelly {
    ###
    # Run structural variant caller Delly
    ###
    my $configuration = shift;
    %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my $jobID = "SV_".get_job_id();

    # Skip sv calling if .done file exists
    if (-e "$opt{OUTPUT_DIR}/logs/StructuralVariants.done"){
	print "WARNING: $opt{OUTPUT_DIR}/logs/StructuralVariants.done exists, skipping \n";
	return \%opt;
    }

    # Create output, log, job and tmp directories.
    $delly_out_dir = "$opt{OUTPUT_DIR}/DELLY/";
    $delly_log_dir = "$delly_out_dir/logs/";
    $delly_job_dir = "$delly_out_dir/jobs/";
    $delly_tmp_dir = "$delly_out_dir/tmp/";
    if (! -e $delly_log_dir) {
	make_path($delly_log_dir) or die "Couldn't create directory: $delly_log_dir\n $! \n";
    }
    if (! -e $delly_job_dir) {
	make_path($delly_job_dir) or die "Couldn't create directory: $delly_job_dir\n $! \n";
    }
    if (! -e $delly_tmp_dir) {
	make_path($delly_tmp_dir) or die "Couldn't create directory: $delly_tmp_dir\n $! \n";
    }

    # Get sample bam files and store running jobs
    foreach my $sample (@{$opt{SAMPLES}}){
	my $sampleBam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
	push @sampleBams, $sampleBam;
	if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
    	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
        }
    }

    # Get SV types and split settings configured in ini/config
    my @svTypes = split/\t/, $opt{DELLY_SVTYPE};
    my @svSplit = split/\t/, $opt{DELLY_SPLIT};

    my @jobIDs_concat;
    my %chrs;

    for(my $i = 0; $i <= $#svTypes; $i++) {
	my $type = $svTypes[$i];

	# Skip SV calling for type if .done exist
	if (-e "$delly_log_dir/DELLY_$type.done"){
	    print "WARNING: $delly_log_dir/DELLY_$type.done exists, skipping \n";
	    next;
	}

	# Split per chromosome
	if ($svSplit[$i] eq "yes") {
	    get_chrs_from_dict(\%chrs) unless scalar(keys %chrs);
	    my ( $jobIDs_chunks, $logFiles );
	    ( $jobIDs_chunks, $logFiles ) = create_interchromosomal_chunks(\@sampleBams, \%chrs, $type) if $type eq "TRA";
	    ( $jobIDs_chunks, $logFiles ) = create_intrachromosomal_chunks(\@sampleBams, \%chrs, $type) if $type =~ /DEL|DUP|INV/;

	    # Translocation jobs
	    if ($type eq "TRA") {
		my $jobID = "CONVERT_".get_job_id();
		my $convert_file = "$delly_job_dir/$type\_".$jobID.".sh";
		open CONVERT, ">$convert_file";
		print CONVERT "#!/bin/bash\n\n";
		print CONVERT "FILES=(".join(" ", @$logFiles).")\n";
		print CONVERT "FINISHED=()\n";
		print CONVERT "FAILED=()\n";
		print CONVERT "for FILE in \${FILES[\@]}\n";
		print CONVERT "do\n";
		print CONVERT "\tTAIL=`tail -n 1 \$FILE`\n";
		print CONVERT "\tif [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print CONVERT "\t\tDONE_FILE=\${FILE//DELLY_[a-zA-Z0-9]*.out/DELLY.done}\n";
		print CONVERT "\t\ttouch \$DONE_FILE\n";
		print CONVERT "\t\tVCF_FILE=\$FILE\n";
		print CONVERT "\t\tVCF_FILE=\${VCF_FILE//\/logs\//\/tmp\/}\n";
		print CONVERT "\t\tVCF_FILE=\${VCF_FILE//_DELLY_[a-zA-Z0-9]*.out/.vcf}\n";
		print CONVERT "\t\tif [ -f \"\$VCF_FILE\" ]\n";
		print CONVERT "\t\tthen\n";
		print CONVERT "\t\t\techo \$VCF_FILE >> $delly_tmp_dir/$type\_vcf_files.txt\n";
		print CONVERT "\t\tfi\n";
		print CONVERT "\t\tFINISHED+=(\$FILE)\n";
	        print CONVERT "\telse\n";
	        print CONVERT "\t\tFAILED+=(\$FILE)\n";
	        print CONVERT "\tfi\n";
	        print CONVERT "done\n\n";
	        print CONVERT "if [[ \${#FAILED[@]} > 0 ]] ; then\n";
	        print CONVERT "\t>&2 echo \"error\"\n";
	        print CONVERT "else\n";
	        print CONVERT "\t$opt{VCFTOOLS_PATH}/vcf-concat -f $delly_tmp_dir/$type\_vcf_files.txt | $opt{VCFTOOLS_PATH}/vcf-sort -c > $delly_tmp_dir/$runName\_$type.vcf\n";
	        print CONVERT "\t$opt{IAP_PATH}/scripts/delly_TRA_convert.pl $delly_tmp_dir/$runName\_$type.vcf\n";
	        print CONVERT "fi\n";
		close CONVERT;
		
		system "qsub -q $opt{DELLY_MERGE_QUEUE} -m a -M $opt{MAIL} -P $opt{CLUSTER_PROJECT} -o $delly_log_dir/$type\_CONVERT.out -e $delly_log_dir/$type\_CONVERT.err -N $jobID -hold_jid ".join(",",@$jobIDs_chunks). " $convert_file";
		push @$jobIDs_chunks, $jobID;
		
		my $jobID2 = "VCF_CONCAT_".get_job_id();
	        my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID2.".sh";
		open VCF_CONCAT, ">$vcf_concat_file";
		print VCF_CONCAT "#!/bin/bash\n\n";
		print VCF_CONCAT "FILE=$delly_log_dir/$type\_CONVERT.out\n";
		print VCF_CONCAT "TAIL=`tail -n 1 \$FILE`\n";
		print VCF_CONCAT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print VCF_CONCAT "\t$opt{VCFTOOLS_PATH}/vcf-sort -c $delly_tmp_dir/$runName\_$type\_CONVERT.vcf > $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf\n";
		print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf $delly_out_dir/$runName\_$type.vcf\n";
	        print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
	        print VCF_CONCAT "fi\n\n";
		close VCF_CONCAT;
		
		system "qsub -q $opt{DELLY_MERGE_QUEUE} -m a -M $opt{MAIL} -P $opt{CLUSTER_PROJECT} -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID2 -hold_jid ".join(",",@$jobIDs_chunks). " $vcf_concat_file";

		push @jobIDs_concat, $jobID2;
	    # Other sv types
	    } else {
		my $jobID = "VCF_CONCAT_".get_job_id();
	        my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID.".sh";
		open VCF_CONCAT, ">$vcf_concat_file";
		print VCF_CONCAT "#!/bin/bash\n\n";
		print VCF_CONCAT "FILES=(".join(" ", @$logFiles).")\n";
		print VCF_CONCAT "FINISHED=()\n";
		print VCF_CONCAT "FAILED=()\n";
		print VCF_CONCAT "for FILE in \${FILES[\@]}\n";
		print VCF_CONCAT "do\n";
		print VCF_CONCAT "\tTAIL=`tail -n 1 \$FILE`\n";
		print VCF_CONCAT "\tif [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print VCF_CONCAT "\t\tDONE_FILE=\${FILE//DELLY_[a-zA-Z0-9]*.out/DELLY.done}\n";
		print VCF_CONCAT "\t\ttouch \$DONE_FILE\n";
		print VCF_CONCAT "\t\tVCF_FILE=\$FILE\n";
		print VCF_CONCAT "\t\tVCF_FILE=\${VCF_FILE//\/logs\//\/tmp\/}\n";
		print VCF_CONCAT "\t\tVCF_FILE=\${VCF_FILE//_DELLY_[a-zA-Z0-9]*.out/.vcf}\n";
		print VCF_CONCAT "\t\tif [ -f \"\$VCF_FILE\" ]\n";
		print VCF_CONCAT "\t\tthen\n";
		print VCF_CONCAT "\t\t\techo \$VCF_FILE >> $delly_tmp_dir/$type\_vcf_files.txt\n";
		print VCF_CONCAT "\t\tfi\n";
		print VCF_CONCAT "\t\tFINISHED+=(\$FILE)\n";
	        print VCF_CONCAT "\telse\n";
	        print VCF_CONCAT "\t\tFAILED+=(\$FILE)\n";
	        print VCF_CONCAT "\tfi\n";
	        print VCF_CONCAT "done\n\n";
	        print VCF_CONCAT "if [[ \${#FAILED[@]} > 0 ]] ; then\n";
	        print VCF_CONCAT "\t>&2 echo \"error\"\n";
	        print VCF_CONCAT "else\n";
	        print VCF_CONCAT "\t$opt{VCFTOOLS_PATH}/vcf-concat -f $delly_tmp_dir/$type\_vcf_files.txt | $opt{VCFTOOLS_PATH}/vcf-sort -c > $delly_tmp_dir/$runName\_$type.vcf\n";
	        print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type.vcf $delly_out_dir/$runName\_$type.vcf\n";
	        print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
	        print VCF_CONCAT "fi\n\n";
	        close VCF_CONCAT;
	        system "qsub -q $opt{DELLY_MERGE_QUEUE} -m a -M $opt{MAIL} -P $opt{CLUSTER_PROJECT} -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID -hold_jid ".join(",",@$jobIDs_chunks). " $vcf_concat_file";
	        push @jobIDs_concat, $jobID;
	    }
	# Non split jobs
	} else {
	    my $jobID = "DELLY_".get_job_id();
	    my $dellyFile = "$delly_job_dir/$type\_".$jobID.".sh";
	    my $logFile = $dellyFile;
	    $logFile =~ s/jobs/logs/;
	    $logFile =~ s/.sh$/.out/;

	    submit_delly($dellyFile, $jobID, $type, "", "$delly_tmp_dir/$runName\_$type.vcf");
	    # Translocation jobs
	    if ($type eq "TRA") {
		    my $jobID2 = "CONVERT_".get_job_id();
		    my $convert_file = "$delly_job_dir/$type\_".$jobID2.".sh";
		    open CONVERT, ">$convert_file";
		    print CONVERT "#!/bin/bash\n\n";
		    print CONVERT "TAIL=`tail -n 1 $logFile`\n";
    		    print CONVERT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
    		    print CONVERT "\t$opt{IAP_PATH}/scripts/delly_TRA_convert.pl $delly_tmp_dir/$runName\_$type.vcf\n";
    		    print CONVERT "fi\n";
    		    close CONVERT;

		    system "qsub -q $opt{DELLY_MERGE_QUEUE} -m a -M $opt{MAIL} -P $opt{CLUSTER_PROJECT} -o $delly_log_dir/$type\_CONVERT.out -e $delly_log_dir/$type\_CONVERT.err -N $jobID2 -hold_jid $jobID $convert_file";

		    my $jobID3 = "VCF_CONCAT_".get_job_id();
		    my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID3.".sh";
		    open VCF_CONCAT, ">$vcf_concat_file";
		    print VCF_CONCAT "#!/bin/bash\n\n";
    		    print VCF_CONCAT "FILE=$delly_log_dir/$type\_CONVERT.out\n";
		    print VCF_CONCAT "TAIL=`tail -n 1 \$FILE`\n";
		    print VCF_CONCAT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
		    print VCF_CONCAT "\t$opt{VCFTOOLS_PATH}/vcf-sort -c $delly_tmp_dir/$runName\_$type\_CONVERT.vcf > $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf\n";
		    print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf $delly_out_dir/$runName\_$type.vcf\n";
	    	    print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
	    	    print VCF_CONCAT "fi\n\n";
		    close VCF_CONCAT;

		    system "qsub -q $opt{DELLY_MERGE_QUEUE} -m a -M $opt{MAIL} -P $opt{CLUSTER_PROJECT} -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID3 -hold_jid $jobID2 $vcf_concat_file";
		
		    push @jobIDs_concat, $jobID3;
	    # Other sv types
	    } else {
		my $jobID2 = "VCF_CONCAT_".get_job_id();
		my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID2.".sh";
		open VCF_CONCAT, ">$vcf_concat_file";
		print VCF_CONCAT "#!/bin/bash\n\n";
		print VCF_CONCAT "TAIL=`tail -n 1 $logFile`\n";
		print VCF_CONCAT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type.vcf $delly_out_dir/$runName\_$type.vcf\n";
		print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
		print VCF_CONCAT "else\n";
		print VCF_CONCAT "\t>&2 echo \"error\"\n";
		print VCF_CONCAT "fi\n\n";
		close VCF_CONCAT;

		system "qsub -q $opt{DELLY_MERGE_QUEUE} -m a -M $opt{MAIL} -P $opt{CLUSTER_PROJECT} -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID2 -hold_jid $jobID $vcf_concat_file";
		
		push @jobIDs_concat, $jobID2;
	    }
	}
    }
    return(\@jobIDs_concat);
}

### Submit delly jobs
sub submit_delly {
    my ($bashFile, $jobID, $type, $excludeFile, $vcfFile) = @_;
    my ($logFile, $errorFile) = ($bashFile, $bashFile);
    $logFile =~ s/jobs/logs/;
    $logFile =~ s/.sh$/.out/;
    $errorFile =~ s/jobs/logs/;
    $errorFile =~ s/.sh$/.err/;

    open DELLY_SH , ">$bashFile" or die "cannot open file $bashFile\n $! \n";
    print DELLY_SH "#!/bin/bash\n\n";
    print DELLY_SH "export OMP_NUM_THREADS=".$opt{DELLY_THREADS}."\n";
    print DELLY_SH "$opt{DELLY_PATH}/delly";
    print DELLY_SH " -t " . $type;
    print DELLY_SH " -g " . $opt{GENOME};
    print DELLY_SH " -x " . $excludeFile if $excludeFile;
    print DELLY_SH " -q " . $opt{DELLY_MAPQUAL};
    print DELLY_SH " -s " . $opt{DELLY_MAD};
    print DELLY_SH " -m " . $opt{DELLY_FLANK};
    print DELLY_SH " -u " . $opt{DELLY_GENO_QUAL};
    print DELLY_SH " -v " . $opt{DELLY_VCF_GENO} if $opt{DELLY_VCF_GENO};
    print DELLY_SH " -o " . $vcfFile;
    print DELLY_SH " ".join(" ", @sampleBams);
    close DELLY_SH;

    if (@runningJobs) {
	system "qsub -q $opt{DELLY_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{DELLY_THREADS} -P $opt{CLUSTER_PROJECT} -o $logFile -e $errorFile -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
        system "qsub -q $opt{DELLY_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{DELLY_THREADS} -P $opt{CLUSTER_PROJECT} -o $logFile -e $errorFile -N $jobID $bashFile";
    }
    return($logFile);
}

### Create inter chromosomal chunks
sub create_interchromosomal_chunks {
    my ($bams, $chrs, $type) = @_;
    my @jobIDs;
    my @logFiles;
    open VCF_FILES, ">$delly_tmp_dir/$type\_vcf_files.txt";
    foreach my $chr1 (keys %{$chrs}) {
	foreach my $chr2 (keys %{$chrs}) {
	    next unless $chr2 gt $chr1;
	    my $vcfFile = "$delly_tmp_dir/$type\_$chr1\_$chr2.vcf";
	    print VCF_FILES $vcfFile . "\n" if -e $vcfFile;
	    if (-e "$delly_log_dir/$type\_$chr1\_$chr2\_DELLY.done") {
		print "WARNING: $delly_log_dir/$type\_$chr1\_$chr2\_DELLY.done exists, skipping \n";
		next;
	    }
	    my $excludeFile = "$delly_tmp_dir/$type\_$chr1\_$chr2\_exclude.txt";
	    open EXC, ">$excludeFile";
	    foreach my $chrom (keys %{$chrs}) {
		print EXC join("\t", $chrom, -1, $chrs->{$chrom}) . "\n" unless $chrom =~ /^$chr1$/ or $chrom =~ /^$chr2$/;
	    }
	    close EXC;
	    my $jobID = "DELLY_".get_job_id();
	    my $dellyFile = "$delly_job_dir/$type\_$chr1\_$chr2\_".$jobID.".sh";
	    push @jobIDs, $jobID;
	    my ( $logFile ) = submit_delly($dellyFile, $jobID, $type, $excludeFile, $vcfFile);
	    push @logFiles, $logFile;
	}
    }
    close VCF_FILES;
    return(\@jobIDs, \@logFiles);
}

### Create intra chromosomal chunks
sub create_intrachromosomal_chunks {
    my ($bams, $chrs, $type) = @_;
    my @jobIDs;
    my @logFiles;
    open VCF_FILES, ">$delly_tmp_dir/$type\_vcf_files.txt";
    foreach my $chr (keys %{$chrs}) {
	my $vcfFile = "$delly_tmp_dir/$type\_$chr.vcf";
	print VCF_FILES $vcfFile . "\n" if -e $vcfFile;
	if ( -e "$delly_log_dir/$type\_$chr\_DELLY.done" ) {
	    print "WARNING: $delly_log_dir/$type\_$chr\_DELLY.done exists, skipping \n";
	    next;
	}
	my $excludeFile = "$delly_tmp_dir/$type\_$chr\_exclude.txt";
	open EXC, ">$excludeFile";
	foreach my $chrom (keys %{$chrs}) {
	    print EXC join("\t", $chrom, -1, $chrs->{$chrom}) . "\n" unless $chrom =~ /^$chr$/;
	}
	close EXC;
	
	my $jobID = "DELLY_".get_job_id();
	my $dellyFile = "$delly_job_dir/$type\_$chr\_".$jobID.".sh";
	push @jobIDs, $jobID;
	my ( $logFile ) = submit_delly($dellyFile, $jobID, $type, $excludeFile, $vcfFile);
	push @logFiles, $logFile;
    }
    close VCF_FILES;
    return(\@jobIDs, \@logFiles);
}

### Get chromosomes from genome.dict
sub get_chrs_from_dict {
    my ($chrs) = @_;
    my $dictFile = $opt{GENOME};
    $dictFile =~ s/.fasta$/.dict/;
    open DICT, $dictFile;
    while(<DICT>) {
	chomp;
	my ($chr, $length) = ($1, $2) if $_ =~ /SN:(\w+)\s*LN:(\d+)/;
	$chrs->{$chr} = $length if $chr;
    }
    close DICT;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;