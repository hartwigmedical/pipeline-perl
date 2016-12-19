package HMF::Pipeline::Mapping;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Config::Validate qw(parseFastqName verifyBai verifyFlagstat);
use HMF::Pipeline::Metadata qw(linkExtraArtefact);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run runBamPrep);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING MAPPING ###";

    die "GENOME: $opt->{GENOME} does not exist!" if !-f $opt->{GENOME};
    die "GENOME BWT: $opt->{GENOME}.bwt does not exist!" if !-f "$opt->{GENOME}.bwt";
    die "GENOME FAI: $opt->{GENOME}.fai does not exist!" if !-f "$opt->{GENOME}.fai";

    my $samples = {};
    foreach my $input_fastq (keys %{$opt->{FASTQ}}) {
        my $fastq = parseFastqName($input_fastq);
        say "Skipping R2 sample $input_fastq" and next if $input_fastq eq $fastq->{R2};
        if (exists $opt->{FASTQ}->{$fastq->{R2}}) {
            say "Switching to paired end mode!";
        } else {
            say "Switching to fragment mode!";
            $opt->{SINGLE_END} = 1;
        }

        my $out_dir = catfile($opt->{OUTPUT_DIR}, $fastq->{sampleName});
        $samples->{$fastq->{sampleName}}->{dirs} = createDirs($out_dir, mapping => "mapping");
        say "Creating $samples->{$fastq->{sampleName}}{dirs}{mapping}/$fastq->{coreName}_sorted.bam with:";
        createIndividualMappingJobs($opt, $fastq, $samples);
    }

    say "";
    while (my ($sample, $info) = each %{$samples}) {
        $opt->{BAM_FILES}->{$sample} = "${sample}_dedup.bam";
        say "Creating $opt->{BAM_FILES}->{$sample}";

        my $done_file = catfile($info->{dirs}->{log}, "Mapping_${sample}.done");
        if (-f $done_file) {
            say "WARNING: $done_file exists, skipping";
            next;
        }

        my $bams = join " ", map { $_->{file} } @{$info->{jobs}};
        my $output_bam = catfile($info->{dirs}->{mapping}, "${sample}_dedup.bam");
        my $output_flagstat = catfile($info->{dirs}->{mapping}, "${sample}_dedup.flagstat");
        my $hold_jids = join ",", map { $_->{job_id} } @{$info->{jobs}};
        my $job_id = "MergeMarkdup_${sample}_" . getId();
        my $bash_file = catfile($info->{dirs}->{job}, "${job_id}.sh");

        writeFromTemplate(
            "MergeMarkdup.sh.tt", $bash_file,
            sample => $sample,
            bams => $bams,
            output_bam => $output_bam,
            output_flagstat => $output_flagstat,
            dirs => $info->{dirs},
            opt => $opt
        );

        my $qsub = qsubTemplate($opt, "MARKDUP");
        my $stdout = catfile($info->{dirs}->{log}, "MergeMarkdup_${sample}.out");
        my $stderr = catfile($info->{dirs}->{log}, "MergeMarkdup_${sample}.err");
        system("$qsub -o $stdout -e $stderr -N $job_id -hold_jid $hold_jids $bash_file");

        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }
    return;
}

sub createIndividualMappingJobs {
    my ($opt, $fastq, $samples) = @_;

    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM, $RG_PU) = ('ILLUMINA', $fastq->{coreName}, $fastq->{sampleName}, $fastq->{sampleName}, $fastq->{flowcellID});

    my %jid = (
        mapping => "Map_$fastq->{coreName}_" . getId(),
        mapping_flagstat => "MapFS_$fastq->{coreName}_" . getId(),
        sort => "Sort_$fastq->{coreName}_" . getId(),
        sort_flagstat => "SortFS_$fastq->{coreName}_" . getId(),
        check_clean => "CheckAndClean_$fastq->{coreName}_" . getId(),
    );

    my $dirs = $samples->{$fastq->{sampleName}}->{dirs};
    push @{$samples->{$fastq->{sampleName}}->{jobs}}, {
        job_id => $jid{check_clean},
        file => catfile($dirs->{mapping}, "$fastq->{coreName}_sorted.bam"),
        };

    my $done_file = catfile($dirs->{mapping}, "$fastq->{coreName}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my $stdout = catfile($dirs->{log}, "Mapping_$fastq->{coreName}.out");
    my $stderr = catfile($dirs->{log}, "Mapping_$fastq->{coreName}.err");
    my $core_file = catfile($dirs->{mapping}, $fastq->{coreName});

    $done_file = catfile($dirs->{log}, "$fastq->{coreName}_bwa.done");
    if (!-f $done_file) {
        say $fastq->{R2} ? "\t$fastq->{R1}\n\t$fastq->{R2}" : "\t$fastq->{R1}";

        my $bash_file = catfile($dirs->{job}, "$jid{mapping}.sh");
        writeFromTemplate(
            "PerLaneMap.sh.tt", $bash_file,
            coreName => $fastq->{coreName},
            sampleName => $fastq->{sampleName},
            R1 => $fastq->{R1},
            R2 => $fastq->{R2},
            RG_ID => $RG_ID,
            RG_SM => $RG_SM,
            RG_PL => $RG_PL,
            RG_LB => $RG_LB,
            RG_PU => $RG_PU,
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "MAPPING");
        system("$qsub -o $stdout -e $stderr -N $jid{mapping} $bash_file");
    } else {
        say "WARNING: $done_file exists, skipping BWA";
    }

    if (!-s "${core_file}.flagstat") {
        my $bash_file = catfile($dirs->{job}, "$jid{mapping_flagstat}.sh");
        writeFromTemplate(
            "PerLaneMapFS.sh.tt", $bash_file,
            sampleName => $fastq->{sampleName},
            coreName => $fastq->{coreName},
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "FLAGSTAT");
        system("$qsub -o $stdout -e $stderr -N $jid{mapping_flagstat} -hold_jid $jid{mapping} $bash_file");
    } else {
        say "\t${core_file}.flagstat exists and is not empty, skipping BWA flagstat";
    }

    if (!-s "${core_file}_sorted.bam") {
        my $bash_file = catfile($dirs->{job}, "$jid{sort}.sh");
        writeFromTemplate(
            "PerLaneSort.sh.tt", $bash_file,
            coreName => $fastq->{coreName},
            sampleName => $fastq->{sampleName},
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "MAPPING");
        system("$qsub -o $stdout -e $stderr -N $jid{sort} -hold_jid $jid{mapping} $bash_file");
    } else {
        say "\t${core_file}_sorted.bam exists and is not empty, skipping sort";
    }

    if (!-s "${core_file}_sorted.flagstat") {
        my $bash_file = catfile($dirs->{job}, "$jid{sort_flagstat}.sh");
        writeFromTemplate(
            "PerLaneSortFS.sh.tt", $bash_file,
            sampleName => $fastq->{sampleName},
            coreName => $fastq->{coreName},
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "FLAGSTAT");
        system("$qsub -o $stdout -e $stderr -N $jid{sort_flagstat} -hold_jid $jid{sort} $bash_file");
    } else {
        say "\t${core_file}_sorted.flagstat exists and is not empty, skipping sorted BAM flagstat";
    }

    my $bash_file = catfile($dirs->{job}, "$jid{check_clean}.sh");
    writeFromTemplate(
        "PerLaneCheckAndClean.sh.tt", $bash_file,
        sampleName => $fastq->{sampleName},
        coreName => $fastq->{coreName},
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "MAPPING");
    system("$qsub -o $stdout -e $stderr -N $jid{check_clean} -hold_jid $jid{mapping_flagstat},$jid{sort_flagstat} $bash_file");
    return;
}

sub runBamPrep {
    my ($opt) = @_;

    $opt->{MAPPING} = "no";
    $opt->{PRESTATS} = "no";
    say "\n### SCHEDULING BAM PREP ###";

    while (my ($sample, $input_bams) = each %{$opt->{SAMPLES}}) {
        my $input_bam = ${$input_bams}[0];
        my $out_dir = catfile($opt->{OUTPUT_DIR}, $sample);
        my $dirs = createDirs($out_dir, mapping => "mapping");

        (my $input_bai = $input_bam) =~ s/\.bam$/.bam.bai/;
        (my $input_flagstat = $input_bam) =~ s/\.bam$/.flagstat/;

        my $bam_file = "${sample}.bam";
        $opt->{BAM_FILES}->{$sample} = $bam_file;
        my $sample_bam = catfile($dirs->{mapping}, $bam_file);
        my $sample_bai = catfile($dirs->{mapping}, "${sample}.bam.bai");
        my $sample_flagstat = catfile($dirs->{mapping}, "${sample}.flagstat");

        my $bai_good = verifyBai($input_bai, $input_bam, $opt);
        my $flagstat_good = verifyFlagstat($input_flagstat, $input_bam);

        symlink($input_bam, $sample_bam);
        $bai_good and symlink($input_bai, $sample_bai);
        $flagstat_good and symlink($input_flagstat, $sample_flagstat);

        next if $bai_good and $flagstat_good;

        my $job_id = "PrepBam_${sample}_" . getId();
        my $bash_file = catfile($dirs->{job}, "$job_id.sh");

        writeFromTemplate(
            "PrepBam.sh.tt", $bash_file,
            sample => $sample,
            sample_bam => $sample_bam,
            sample_bai => $sample_bai,
            sample_flagstat => $sample_flagstat,
            dirs => $dirs,
            opt => $opt,
        );

        my $qsub = qsubTemplate($opt, "MAPPING");
        system "$qsub -o $dirs->{log}/PrepBam_${sample}.out -e $dirs->{log}/PrepBam_${sample}.err -N $job_id $bash_file";
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
    }
    return;
}

1;
