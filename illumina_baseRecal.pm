package illumina_baseRecal;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge qw(getJobId jobNative qsubJava qsubTemplate);
use illumina_template qw(from_template);


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
        (my $recal_bambai = $sample_bam) =~ s/\.bam$/.recalibrated.bam.bai/;
        (my $recal_flagstat = $sample_bam) =~ s/\.bam$/.recalibrated.flagstat/;

        my $out_dir = catfile($opt->{OUTPUT_DIR}, $sample);
        my $dirs = {
            out => $out_dir,
            log => catfile($out_dir, "logs"),
            mapping => catfile($out_dir, "mapping"),
            tmp => catfile($out_dir, "tmp"),
            job => catfile($out_dir, "jobs"),
        };

        my $sample_bam_path = catfile($dirs->{mapping}, $sample_bam);
        say "\t$sample_bam_path";

        $opt->{BAM_FILES}->{$sample} = $recal_bam;

        my $done_file = catfile($dirs->{log}, "BaseRecalibration_${sample}.done");
        if (-f $done_file) {
            say "\tWARNING: $done_file exists, skipping";
            next;
        }

        my $job_native = jobNative($opt, "BASERECALIBRATION");
        my $job_id = "BQSR_${sample}_" . getJobId();
        my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
        my $stdout = catfile($dirs->{log}, "BaseRecalibration_${sample}.out");
        my $stderr = catfile($dirs->{log}, "BaseRecalibration_${sample}.err");
        my $qsub = qsubJava($opt, "BASERECALIBRATION_MASTER");

        from_template("BaseRecalibration.sh.tt", $bash_file,
                      sample => $sample,
                      sample_bam => $sample_bam,
                      sample_bam_path => $sample_bam_path,
                      job_native => $job_native,
                      known_files => $known_files,
                      dirs => $dirs,
                      opt => $opt
                  );

        if (@{$opt->{RUNNING_JOBS}->{$sample}}) {
            system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @{$opt->{RUNNING_JOBS}->{$sample}}) . " $bash_file";
        } else {
            system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
        }

        my $job_id_fs = "BQSRFS_${sample}_" . getJobId();
        $bash_file = catfile($dirs->{job}, "${job_id_fs}.sh");

        from_template("BaseRecalibrationFS.sh.tt", $bash_file,
                      sample => $sample,
                      sample_bam => $sample_bam,
                      sample_flagstat => $sample_flagstat,
                      recal_bam => $recal_bam,
                      recal_bai => $recal_bai,
                      recal_bambai => $recal_bambai,
                      recal_flagstat => $recal_flagstat,
                      dirs => $dirs,
                      opt => $opt
                  );

        $qsub = qsubTemplate($opt, "FLAGSTAT");
        system "$qsub -o $stdout -e $stderr -N $job_id_fs -hold_jid $job_id $bash_file";

        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id;
        push @{$opt->{RUNNING_JOBS}->{$sample}}, $job_id_fs;
    }
    return;
}

1;
