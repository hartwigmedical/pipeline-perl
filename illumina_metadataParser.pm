package illumina_metadataParser;

use 5.16.0;
use strict;
use warnings;

use File::Spec::Functions;
use JSON;
use Carp;


BEGIN {
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT = qw(metadataParse);
}


sub metadataParse {
    my $directory = shift || return;
    my $metadata_file = catfile($directory, "metadata");
    my $json_conf = do {
        open my $json_fh, "<:encoding(UTF-8)", $metadata_file or confess "Can't open $metadata_file";
        local $/;
        <$json_fh>;
    };
    my $config = decode_json($json_conf);
    return $config;
}

1;
