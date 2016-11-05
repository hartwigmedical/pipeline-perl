package illumina_template;

use FindBin::libs;
use discipline;

use Carp;
use File::Spec::Functions;
use FindBin;
use Template;
use Env qw($TEMPLATES);

use illumina_config qw(pipelinePath);

use parent qw(Exporter);
our @EXPORT_OK = qw(from_template);


sub from_template {
    my ($name, $output_file, %data) = @_;

    my $t = Template->new(INCLUDE_PATH => templateDir(), STRICT => 1);
    open my $tout, ">", $output_file or confess "Unable to open $output_file for writing";
    $t->process($name, \%data, \*$tout) or confess $t->error();
    close $tout;
    return;
}

sub templateDir {
    my $source_template_dir = catfile(pipelinePath(), "templates");
    return $TEMPLATES ? $TEMPLATES : $source_template_dir;
}

1;
