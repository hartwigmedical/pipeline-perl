package illumina_metadata;

use FindBin::libs;
use discipline;

use File::Spec::Functions qw(:ALL);
use File::Basename;

use JSON;
use Carp;

use parent qw(Exporter);
our @EXPORT_OK = qw(parse metaSampleName linkArtefact writeLinks);


sub readJson {
    my ($path) = @_;

    my $json_conf;
    {
        open my $json_fh, "<:encoding(UTF-8)", $path or confess "Can't open $path: $!";
        local $/ = undef;
        $json_conf = <$json_fh>;
        close $json_fh;
    }
    my $config = decode_json($json_conf);
    return $config;
}

sub writeJson {
    my ($path, $data) = @_;

    my $json = to_json($data, {utf8 => 1, pretty => 1});
    open my $json_fh, ">:encoding(UTF-8)", $path or confess "Can't open $path: $!";
    print $json_fh $json;
    close $json_fh;
    return;
}

sub parse {
    my ($opt, $parse_opts) = @_;
    $parse_opts //= {};

    my $metadata_path = catfile($opt->{OUTPUT_DIR}, "metadata");
    if (-z $metadata_path and exists $parse_opts->{required} and not $parse_opts->{required}) {
        return {};
    } else {
        return readJson($metadata_path);
    }
}

sub metaSampleName {
    my ($sample, $opt) = @_;

    my %name_map = reverse %{parse($opt, {required => 0})};
    $name_map{$sample} //= "sample";
    return $name_map{$sample};
}

sub linkArtefact {
    my ($source_path, $canonical_name, $opt) = @_;

    $opt->{LINKS} = {} if not exists $opt->{LINKS};
    $opt->{LINKS}->{$canonical_name} = $source_path;
    return;
}

sub writeLinks {
    my ($opt) = @_;
    my $outfile = catfile($opt->{OUTPUT_DIR}, "logs", "links.json");
    writeJson($outfile, $opt->{LINKS});
    return;
}

1;
