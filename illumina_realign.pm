package illumina_realign;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge qw(jobNative qsubJava qsubTemplate);
use illumina_jobs;


sub runRealignment {
    my ($opt) = @_;

    say "Running single sample indel realignment for the following BAM-files:";

    my $known_files;
    $known_files = join " ", map { "-known $_" } split '\t', $opt->{REALIGNMENT_KNOWN} if $opt->{REALIGNMENT_KNOWN};

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $sample_bam = $opt->{BAM_FILES}->{$sample};
        (my $sample_flagstat = $sample_bam) =~ s/\.bam$/.flagstat/;
        (my $realigned_bam = $sample_bam) =~ s/\.bam$/.realigned.bam/;
        (my $realigned_bai = $sample_bam) =~ s/\.bam$/.realigned.bai/;
        (my $realigned_flagstat = $sample_bam) =~ s/\.bam$/.realigned.flagstat/;
        (my $cpct_sliced_bam = $sample_bam) =~ s/\.bam$/.realigned.sliced.bam/;

        $opt->{BAM_FILES}->{$sample} = $realigned_bam;

        my $out_dir = catfile($opt->{OUTPUT_DIR}, $sample);
        my $dirs = {
            out => $out_dir,
            log => catfile($out_dir, "logs"),
            tmp => catfile($out_dir, "tmp"),
            job => catfile($out_dir, "jobs"),
            mapping => catfile($out_dir, "mapping"),
        };

        my $sample_bam_path = catfile($dirs->{mapping}, $sample_bam);
        say "\t${sample_bam_path}";

        my $realign_job_id = illumina_jobs::fromTemplate(
            "Realignment",
            $sample,
            qsubJava($opt, "REALIGNMENT_MASTER"),
            $opt->{RUNNING_JOBS}->{$sample},
            $dirs,
            $opt,
            sample => $sample,
            sample_bam => $sample_bam,
            sample_bam_path => $sample_bam_path,
            job_native => jobNative($opt, "REALIGNMENT"),
            known_files => $known_files);

        next unless $realign_job_id;

        my $flagstat_job_id = illumina_jobs::flagstatBam(
            $sample,
            catfile($dirs->{tmp}, $realigned_bam),
            catfile($dirs->{mapping}, $realigned_flagstat),
            [$realign_job_id],
            $dirs,
            $opt);

        my $check_job_id = illumina_jobs::fromTemplate(
            "ReadCountCheck",
            $sample,
            qsubTemplate($opt, "FLAGSTAT"),
            [$flagstat_job_id],
            $dirs,
            $opt,
            sample => $sample,
            pre_flagstat_path => catfile($dirs->{mapping}, $sample_flagstat),
            post_flagstat_path => catfile($dirs->{mapping}, $realigned_flagstat),
            post_bam => $realigned_bam,
            post_bai => $realigned_bai,
            success_template => "RealignmentSuccess.tt");

        push @{$opt->{RUNNING_JOBS}->{$sample}}, $check_job_id;

        my $qc_job_ids = illumina_jobs::prePostSliceAndDiffBams(
            $sample,
            "realign",
            $sample_bam,
            $realigned_bam,
            [$check_job_id],
            $dirs,
            $opt);

        my $cpct_job_id = illumina_jobs::sliceBam(
            $sample,
            $realigned_bam,
            $cpct_sliced_bam,
            "CPCT_Slicing.bed",
            [$check_job_id],
            $dirs,
            $opt);

        push @{$opt->{RUNNING_JOBS}->{slicing}}, @{$qc_job_ids}, $cpct_job_id;
    }
    return;
}

1;
