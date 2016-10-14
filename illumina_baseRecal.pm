package illumina_baseRecal;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge qw(jobNative qsubJava qsubTemplate);
use illumina_jobs;


sub runBaseRecalibration {
    my ($opt) = @_;

    say "Running base recalibration for the following BAM-files:";

    my $known_files;
    $known_files = join " ", map { "-knownSites $_" } split '\t', $opt->{BASERECALIBRATION_KNOWN} if $opt->{BASERECALIBRATION_KNOWN};

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $sample_bam = $opt->{BAM_FILES}->{$sample};
        (my $sample_flagstat = $sample_bam) =~ s/\.bam$/.flagstat/;
        (my $recal_bam = $sample_bam) =~ s/\.bam$/.recalibrated.bam/;
        (my $recal_bai = $sample_bam) =~ s/\.bam$/.recalibrated.bai/;
        (my $recal_flagstat = $sample_bam) =~ s/\.bam$/.recalibrated.flagstat/;
        (my $cpct_sliced_bam = $sample_bam) =~ s/\.bam$/.recalibrated.sliced.bam/;

        $opt->{BAM_FILES}->{$sample} = $recal_bam;

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

        my $recal_job_id = illumina_jobs::fromTemplate(
            "BaseRecalibration",
            $sample,
            qsubJava($opt, "BASERECALIBRATION_MASTER"),
            $opt->{RUNNING_JOBS}->{$sample},
            $dirs,
            $opt,
            sample => $sample,
            sample_bam => $sample_bam,
            sample_bam_path => $sample_bam_path,
            job_native => jobNative($opt, "BASERECALIBRATION"),
            known_files => $known_files);

        next unless $recal_job_id;

        my $flagstat_job_id = illumina_jobs::flagstatBam(
            $sample,
            catfile($dirs->{tmp}, $recal_bam),
            catfile($dirs->{mapping}, $recal_flagstat),
            [$recal_job_id],
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
            post_flagstat_path => catfile($dirs->{mapping}, $recal_flagstat),
            post_bam => $recal_bam,
            post_bai => $recal_bai,
            success_template => "BaseRecalibrationSuccess.tt");

        push @{$opt->{RUNNING_JOBS}->{$sample}}, $check_job_id;

        my $qc_job_ids = illumina_jobs::prePostSliceAndDiffBams(
            $sample,
            "recal",
            $sample_bam,
            $recal_bam,
            [$check_job_id],
            $dirs,
            $opt);

        my $cpct_job_id = illumina_jobs::sliceBam(
            $sample,
            $recal_bam,
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
