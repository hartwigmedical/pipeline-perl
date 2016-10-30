package illumina_germlineCalling;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_sge qw(jobNative qsubJava);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);
use illumina_metadata;


sub runVariantCalling {
    my ($opt) = @_;

    my @sample_bams;
    my @running_jobs;
    my $job_id = "GermlineCalling_" . getJobId();

    # maintain backward-compatibility with old naming for now, useful for re-running somatics without re-running germline
    if (-f "$opt->{OUTPUT_DIR}/logs/GermlineCaller.done" || -f "$opt->{OUTPUT_DIR}/logs/VariantCaller.done") {
	    say "WARNING: $opt->{OUTPUT_DIR}/logs/GermlineCaller.done exists, skipping";
	    return;
    }

    my $gvcf_dir = catfile($opt->{OUTPUT_DIR}, "gvcf");
    if ((!-d $gvcf_dir && $opt->{CALLING_GVCF} eq "yes")) {
	    mkdir($gvcf_dir) or die "Couldn't create directory $gvcf_dir: $!";
    }

    my $job_native = jobNative($opt, "CALLING");
    my $command = "java -Xmx".$opt->{CALLING_MASTER_MEM}."G -Djava.io.tmpdir=$opt->{OUTPUT_DIR}/tmp -jar $opt->{QUEUE_PATH}/Queue.jar ";
    $command .= "-jobQueue $opt->{CALLING_QUEUE} -jobNative \"$job_native\" -jobRunner GridEngine -jobReport $opt->{OUTPUT_DIR}/logs/GermlineCaller.jobReport.txt -memLimit $opt->{CALLING_MEM} ";

    $command .= "-S $opt->{OUTPUT_DIR}/QScripts/$opt->{CALLING_SCALA} ";
    if ($opt->{CALLING_UGMODE}) {
	    $command .= " -glm $opt->{CALLING_UGMODE} ";
    }

    $command .= "-R $opt->{GENOME} -O $opt->{RUN_NAME} -mem $opt->{CALLING_MEM} -nct $opt->{CALLING_THREADS} -nsc $opt->{CALLING_SCATTER} -stand_call_conf $opt->{CALLING_STANDCALLCONF} -stand_emit_conf $opt->{CALLING_STANDEMITCONF} ";

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $sample_bam = "$opt->{OUTPUT_DIR}/$sample/mapping/$opt->{BAM_FILES}->{$sample}";

        $command .= "-I $sample_bam ";
        push(@sample_bams, $sample_bam);

        if (@{$opt->{RUNNING_JOBS}->{$sample}}) {
            push(@running_jobs, @{$opt->{RUNNING_JOBS}->{$sample}});
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

    my $bash_file = catfile($opt->{OUTPUT_DIR}, "jobs", "${job_id}.sh");
    my $log_dir = catfile($opt->{OUTPUT_DIR}, "logs");
    my $stdout = catfile($log_dir, "GermlineCaller_$opt->{RUN_NAME}.out");
    my $stderr = catfile($log_dir, "GermlineCaller_$opt->{RUN_NAME}.err");

    from_template("GermlineCalling.sh.tt", $bash_file,
                  command => $command,
                  gvcf_dir => $gvcf_dir,
                  sample_bams => \@sample_bams,
                  opt => $opt);

    my $qsub = qsubJava($opt, "CALLING_MASTER");
    if (@running_jobs) {
        system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
    }

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }

    linkArtefacts($gvcf_dir, $opt);
    return;
}

# naming dependent on GermlineCaller.scala, could fix to be explicit
sub linkArtefacts {
    my ($gvcf_dir, $opt) = @_;

    if ($opt->{CALLING_GVCF} eq "yes") {
        foreach my $sample (keys %{$opt->{SAMPLES}}) {
            my $bam_file = $opt->{BAM_FILES}->{$sample};
            (my $gvcf_file = $bam_file) =~ s/\.bam$/.g.vcf.gz/;
            my $gvcf_path = catfile($gvcf_dir, $gvcf_file);
            my $portal_name = illumina_metadata::portalName($sample, $opt);
            illumina_metadata::linkArtefact($gvcf_path, "${sample}.g.vcf.gz", "${portal_name} gVCF", $opt);
        }
    }
    my $germline_vcf_path = catfile($opt->{OUTPUT_DIR}, "$opt->{RUN_NAME}.raw_variants.vcf");
    illumina_metadata::linkArtefact($germline_vcf_path, "germline.vcf", "Germline VCF", $opt);
    return;
}

1;
