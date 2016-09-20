package illumina_template;
require Exporter;
@ISA = qw(Exporter);
@EXPORT= qw(from_template);

use 5.16.0;
use strict;
use warnings;

use Carp;
use File::Spec::Functions;
use Template;
use Env qw($TEMPLATES);

use FindBin;

my $template_dir = $TEMPLATES ? $TEMPLATES : catfile("$FindBin::Bin", "templates");

sub from_template {
    my $tname = shift || return undef;
    my $outname = shift || return undef;
    my %data = @_;
    my $t = Template->new(INCLUDE_PATH => $template_dir);

    my $tout;
    open $tout, ">$outname" or confess "Unable to open $tout for writing";
    $t->process($tname, \%data, \*$tout) or confess $t->error();
    close $tout;
};

1;
