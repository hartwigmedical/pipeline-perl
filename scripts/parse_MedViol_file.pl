#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#parse (or run) output from PhaseByTransmission MedelianViolation file
#CHROM	POS	AC	FAMILY	TP	MOTHER_GT	MOTHER_DP	MOTHER_AD	MOTHER_PL	FATHER_GT	FATHER_DP	FATHER_AD	FATHER_PL	CHILD_GT	CHILD_DP	CHILD_AD	CHILD_PL

my $pwd = `pwd`;
chomp($pwd);

my (%vars);
my ($help, $debug,$mv_file,$output_name, $vcf_file);

my %SET = ( # contains filter settings
    'removeKnowns'     => 0,
);

my %LIM = ( # contains all filter limits
    'min_cov'       => 10,
    'min_av_cov'    => 5,
    'min_aff_pos'   => 1,
    'max_aff_neg'   => 0,
);



GetOptions (
	'h'                => \$help,
	'debug'            => \$debug,
	'in|i=s'           => \$mv_file,
	'vcf=s'			=> \$vcf_file,
	'out|o=s'          => \$output_name,
	
	'min_cov=i'        => \$LIM{ 'min_cov' },
	'min_av_cov=i'     => \$LIM{ 'min_av_cov' },
	'min_aff=f'        => \$LIM{ 'min_aff_pos' },
	'max_uaf=f'        => \$LIM{ 'max_aff_neg' },
	
	'novel'            => \$SET{ 'removeKnowns' },
) 
or die "ERROR: Illegal arguments or parameters: @ARGV\n";
die usage() if $help;
die "[ERROR] Provide Mendelian Violationfile?\n" unless $mv_file;
die "[ERROR] No output name given?\n" unless ( $output_name );

open IN, $mv_file or die "cannot open infile:$!\n";
open OUT, ">$output_name.table" or die "cannot open outfile table\n";
open OUTVCF, ">$output_name.vcf" or die "cannot open outfile vcf\n";

my @headers;

while (my $line=<IN>) {
    chomp($line);
    if ($line =~ /^CHROM/) {
	@headers=split("\t",$line);
    }else{
	my @vals = split("\t",$line);
	my $c=0;
	foreach my $col (@headers) {
	    $vars{$vals[0].'_'.$vals[1]}{$col}=$vals[$c] unless $col =~ /CHROM|POS/i;;
	    $c++;
	}
	$vars{$vals[0].'_'.$vals[1]}{full_line} = $line;
    }
}

print scalar(@headers), " headers stored\n";
print scalar keys %vars," positions stored\n";


#OUTPUT
print OUT join("\t",@headers),"\n";

system "grep -P \"#+\" $vcf_file > $output_name.vcf";


foreach my $loc (sort keys %vars) {

    #no data filter
    next if ( ($vars{$loc}{MOTHER_GT} =~ /\./) or ($vars{$loc}{FATHER_GT} =~ /\./) or ($vars{$loc}{CHILD_GT} =~ /\./) );


    #genetic filter
    my ($c_allele1,$c_allele2) = split(/\||\//, $vars{$loc}{CHILD_GT});
    my ($f_allele1,$f_allele2) = split(/\||\//, $vars{$loc}{FATHER_GT});
    my ($m_allele1,$m_allele2) = split(/\||\//, $vars{$loc}{MOTHER_GT});
    
    my %alleles;
    $alleles{$_}++ foreach ( $f_allele1,$f_allele2,$m_allele1,$m_allele2 );
    #print $vars{$loc}{MOTHER_GT},", ",$vars{$loc}{FATHER_GT},"\t",$vars{$loc}{CHILD_GT},"\n";;
    
    #pass only de-novo
    #next if ( ( ($allele1 =~ /$vars{$loc}{MOTHER_GT}/) and ($allele2 =~ /$vars{$loc}{FATHER_GT}/) ) and ( ($allele2 =~ /$vars{$loc}{MOTHER_GT}/) and ($allele1 =~ /$vars{$loc}{FATHER_GT}/) ) );
    next if ( (exists $alleles{$c_allele1}) and (exists $alleles{$c_allele2}) );
    
    
    
    #data filters
    next if $vars{$loc}{CHILD_DP} <= $LIM{'min_av_cov'};
    next if $vars{$loc}{FATHER_DP} <= $LIM{'min_cov'};
    next if $vars{$loc}{MOTHER_DP} <= $LIM{'min_cov'};
    
    print OUT $vars{$loc}{full_line},"\n";
    
    
    #filter lines from org VCf file
    my ($chr,$pos) = split("_",$loc);
    system "grep -P \"^$chr\t$pos\" $vcf_file >> $output_name.vcf";

}

print "done\n";