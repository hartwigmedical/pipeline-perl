#!/usr/bin/perl -w

###################################
### illumina_template.pm
### - Use templates to generate files
###
### Author: awolfs@schubergphilis.com
###################################

package illumina_template;
require Exporter;
@ISA = qw(Exporter);
@EXPORT= qw(from_template);

use strict;
use Text::Template;
use Template;
use Env qw($TEMPLATES);

my $template_dir = $TEMPLATES;

sub from_template {
	my $tname = shift || return undef;
	my $outname = shift || return undef;
	my %data = @_;
	my $t = Template->new(INCLUDE_PATH => $template_dir);

	my $tout;
	open($tout, ">", "$outname") or die "Unable to open $tout for writing";
	$t->process($tname, \%data, \*$tout) or die $t->error();
	close($tout);
};

sub from_template3 {
	my $tname = shift || return undef;
	my $outname = shift || return undef;
	my %data = @_;
	my $t = new Text::Template(TYPE => 'FILE', SOURCE=> "$template_dir/$tname");

	my $tout;
	open($tout, ">", "$outname") or die "Unable to open $tout for writing";
	$t->fill_in(OUTPUT=>\*$tout, HASH => \%data) or die "Unable to process template $template_dir/$tname:$Text::Template::ERROR";  
	close($tout);
};

1;
