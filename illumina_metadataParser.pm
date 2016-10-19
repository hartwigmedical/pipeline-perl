package illumina_metadataParser;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Spec::Functions;
use JSON;
use Carp;

use parent qw(Exporter);
our @EXPORT_OK = qw(metadataParse);


sub metadataParse {
    my $directory = shift || return;
    my $metadata_file = catfile($directory, "metadata");
    my $json_conf;
    {
        open my $json_fh, "<:encoding(UTF-8)", $metadata_file or confess "Can't open $metadata_file";
        local $/ = undef;
        $json_conf = <$json_fh>;
        close $json_fh;
    }
    my $config = decode_json($json_conf);
    return $config;
}

1;
