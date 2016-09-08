package illumina_sge;

use strict;
use warnings;
use POSIX qw(tmpnam);

BEGIN {
    use Exporter;

    our @ISA = ('Exporter');

    our @EXPORT = qw(
                        &qsubTemplate
                        &qsubJava
                        &jobNative
                        &getJobId
                   );
}

sub generic(@) {
    my ($opt, $function) = @_;
    my $qsub = qq|qsub -P $opt->{CLUSTER_PROJECT} -pe threaded $opt->{$function . "_THREADS"} -q $opt->{$function . "_QUEUE"} -l h_rt=$opt->{$function . "_TIME"}|;
    $qsub .= qq|,tmpspace=$opt->{$function . "_TMP"}G| if $opt->{$function . "_TMP"};
    return $qsub;
}

sub qsubTemplate(@) {
    my ($opt, $function) = @_;
    my $qsub = generic($opt, $function) . " -m a -M $opt->{MAIL} -R $opt->{CLUSTER_RESERVATION}";
    return $qsub;
}

sub qsubJava(@) {
    my ($opt, $function) = @_;
    # KODU: h_vmem setting leads to issues with Queue.jar wants to spawn jobs
    # my $h_vmem = (4 + $opt->{$function."_MEM"})."G";
    my $qsub = generic($opt, $function) . " -m a -M $opt->{MAIL} -R $opt->{CLUSTER_RESERVATION}";
    return $qsub
}

sub jobNative(@) {
    my ($opt, $function) = @_;
    # KODU: h_vmem setting leads to issues with Queue.jar wants to spawn jobs
    # my $h_vmem = (4 + $opt->{$function."_MEM"})."G";
    my $qsub = generic($opt, $function);
    $qsub =~ s/^qsub//;
    return $qsub
}

sub getJobId {
    my $id = tmpnam();
    $id =~ s#/tmp/file##;
    return $id;
}

1;
