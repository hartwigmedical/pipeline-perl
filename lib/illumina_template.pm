package illumina_template;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use Carp;
use File::Spec::Functions;
use Template;
use Env qw($TEMPLATES);


use parent qw(Exporter);
our @EXPORT_OK = qw(from_template);


my $template_dir = $TEMPLATES ? $TEMPLATES : catfile("$FindBin::Bin", "templates");


sub from_template {
    my ($name, $output_file, %data) = @_;

    my $t = Template->new(INCLUDE_PATH => $template_dir, STRICT => 1);
    open my $tout, ">", $output_file or confess "Unable to open $output_file for writing";
    $t->process($name, \%data, \*$tout) or confess $t->error();
    close $tout;
    return;
};

1;
