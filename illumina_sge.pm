#
# Copyright (c) Hindrik H.D. Kerstens and
# UMC Utrecht. All rights reserved
#


package illumina_sge;
use strict;
use warnings;

BEGIN {

    use Exporter;

    our @ISA = ('Exporter');

    our @EXPORT = qw (
        &qsubTemplate
        &qsubJava
        &jobNative
   );
}

#EXCEPTIONS
#my %excepts=();
#$excepts{"MEM_THREADS"}{"REALIGNMENT"}=1;
#$excepts{"MEM_THREADS"}{"CALLING"}=1;

sub generic(@){
    my ($opt,$function)=@_;
    my $tmp="";
    $tmp = ",tmpspace=".$$opt{$function."_TMP"}."G" if ($$opt{$function."_TMP"});
    my $qsub = "qsub -P ".$$opt{CLUSTER_PROJECT}." -pe threaded ".$$opt{$function."_THREADS"}." -q ".$$opt{$function."_QUEUE"}." -l h_rt=".$$opt{$function."_TIME"}.$tmp;
    return ($qsub);
}

sub qsubTemplate(@){
    my ($opt,$function)=@_;
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION}." -l h_vmem=".$$opt{$function."_MEM"}."G";
    return ($qsub)
}

sub qsubJava(@){
    my ($opt,$function)=@_;
    my $h_vmem = (4 + $$opt{$function."_MEM"})."G";
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION}." -l h_vmem=".$h_vmem;
    return ($qsub)
}

sub jobNative(@){
    my ($opt,$function)=@_;
    my $h_vmem = (4 + $$opt{$function."_MEM"})."G";
    my $qsub = &generic($opt,$function).",h_vmem=".$h_vmem;
    $qsub=~s/^qsub//;
    return ($qsub)
}

1;
