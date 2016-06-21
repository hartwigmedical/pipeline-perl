#!/usr/bin/perl -w

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

sub generic(@){
    my ($opt,$function)=@_;
    my $tmp="";
    $tmp = ",tmpspace=".$$opt{$function."_TMP"}."G" if ($$opt{$function."_TMP"});
    my $qsub = "qsub -P ".$$opt{CLUSTER_PROJECT}." -pe threaded ".$$opt{$function."_THREADS"}." -q ".$$opt{$function."_QUEUE"}." -l h_rt=".$$opt{$function."_TIME"}.$tmp;
    return ($qsub);
}

sub qsubTemplate(@){
    my ($opt,$function)=@_;
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION};
    return ($qsub)
}

sub qsubJava(@){
    my ($opt,$function)=@_;
    # KODU: h_vmem setting leads to issues with Queue.jar wants to spawn jobs
#    my $h_vmem = (4 + $$opt{$function."_MEM"})."G";
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION};
    return ($qsub)
}

sub jobNative(@){
    my ($opt,$function)=@_;
    # KODU: h_vmem setting leads to issues with Queue.jar wants to spawn jobs
#    my $h_vmem = (4 + $$opt{$function."_MEM"})."G";
    my $qsub = &generic($opt,$function);
    $qsub=~s/^qsub//;
    return ($qsub)
}

1;
