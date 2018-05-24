package HMF::Pipeline::Functions::Template;

use FindBin::libs;
use discipline;

use Carp;
use Env qw($TEMPLATES);
use File::Spec::Functions;
use FindBin;
use Template;

use HMF::Pipeline::Functions::Config;

use parent qw(Exporter);
our @EXPORT_OK = qw(writeFromTemplate);

sub writeFromTemplate {
    my ($name, $output_file, %data) = @_;

    my $t = Template->new(INCLUDE_PATH => templateDir(), STRICT => 1);
    open my $tout, ">", $output_file or confess "unable to open $output_file for writing";
    $t->process($name, \%data, \*$tout) or confess $t->error();
    close $tout;
    return;
}

sub templateDir {
    ## no critic (Subroutines::ProhibitCallsToUnexportedSubs)
    my $source_template_dir = catfile(HMF::Pipeline::Functions::Config::pipelinePath(), "templates");
    ## use critic
    return $TEMPLATES ? $TEMPLATES : $source_template_dir;
}

1;