package illumina_template;

use 5.16.0;
use strict;
use warnings;

use Carp;
use File::Spec::Functions;
use Template;
use Env qw($TEMPLATES);

use FindBin;


my $template_dir = $TEMPLATES ? $TEMPLATES : catfile("$FindBin::Bin", "templates");

BEGIN {
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT= qw(from_template);
}


sub from_template {
    my $name = shift || return;
    my $output_file = shift || return;
    my %data = @_;

    my $t = Template->new(INCLUDE_PATH => $template_dir, STRICT => 1);
    open my $tout, ">", $output_file or confess "Unable to open $output_file for writing";
    $t->process($name, \%data, \*$tout) or confess $t->error();
    close $tout;
};

1;
