#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);


my $email = "i.nijman\@umcutrecht.nl";
my $snpEff_path = "/hpc/cog_bioinf/data/ies/src/snpEff";
my $SAP42_path = "/hpc/cog_bioinf/common_scripts/SAP42-HPC";

my $vcf = $ARGV[0];
die "vcf file invalid or non-existent...\n" unless -e $vcf;

my @jids;

mkdir "SAP42_annotation_parts" or die "cannot create output dir\n";
chdir "SAP42_annotation_parts" or die "cannot change into dir\n";;

my $pwd = `pwd`;
chomp($pwd);
print $pwd,"\n";

print "splitting vcf into chunks...\n";

#system "cp ../$vcf ./";
system "java -Xmx5g -jar $snpEff_path/SnpSift.jar split -l 500 ../$vcf";
print "done; submitting annotation chunks...\n";

foreach my $chunk (<*.vcf>) {
    next if $chunk eq $vcf;
    my $command = "$SAP42_path/annotator.pl -in $chunk -out $chunk\n";

    my $jid = cluster($command,$pwd);
    push(@jids,$jid);

}


#make hold job for merging snv files and cleanup.


#===========================================================================================================================
sub cluster {
    my $jid2 = tmpnam();
    $jid2=~s/\/tmp\/file//;
    my $comm = shift;
    my $pwd = shift;
    
    open OUT, ">$pwd/SAP42ann_$jid2.sh" or die "cannot open file $pwd/HC_$jid2.sh\n\n";
    print OUT "#!/bin/bash\n\n";
    print OUT "cd $pwd\n";
    print OUT "$comm\n";
    
    system "qsub -q veryshort -M $email -m a -l cog_bioinf_mysql=1 -o $pwd -e $pwd -N SAP42ann_$jid2 $pwd/SAP42ann_$jid2.sh"; #require two slots for memory reasons
    return $jid2;


}