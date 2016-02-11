#!/usr/bin/perl -w
use strict;


#run kinship alalysis on vcf input file.
#requires plink, king and vcftools;

#to be incorporated: read in proposed ped file with fam structure and parse results file. warn if kinship values disagree with the ped file!




#Relationship        Kinship coefficient   Coeffcient of relatedness
#Self or Monozygotic twins   0.5000                1.000
#Parent-child        0.2500                0.500
#Full siblings       0.2500                0.500
#Half siblings       0.1250                0.250
#First cousins       0.0625                0.1250
#Unrelated           0.0000                0.0000

my $vcf = $ARGV[0];
die "vcf file invalid or non-existent...\n" unless -e $vcf;

my $cfn = $vcf;
$cfn =~ s/.vcf.gz//;
$cfn =~ s/.vcf//;

my $vcftools="/hpc/cog_bioinf/common_scripts/vcftools/bin/vcftools";
my $plink="/hpc/local/CentOS6/cog_bioinf/plink-1.07-x86_64/plink";
my $king="/hpc/local/CentOS6/cog_bioinf/plink-1.07-x86_64/king";


my $pwd = `pwd`;
chomp($pwd);

mkdir "kinship_analyses" or die "cannot make output dir...\n";
chdir "kinship_analyses";

if ($vcf =~ /vcf.gz$/) {
    system "$vcftools --gzvcf $pwd/$vcf --plink";
}else{
    system "$vcftools --vcf $pwd/$vcf --plink";
}

system "$plink --file out --make-bed --noweb";
system "$king -b plink.bed --kinship";

system "cat king.kin0";
system "cp king.kin0 ../$cfn.kinship";
chdir "../";
system "rm -rf kinship_analyses";
