package illumina_metadata;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Spec::Functions qw(:ALL);
use File::Basename;
use File::chdir;
use JSON;
use Carp;

use parent qw(Exporter);
our @EXPORT_OK = qw(parse linkArtefact);


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

    my $json = to_json($data, { utf8 => 1, pretty => 1 });
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

sub portalName {
    my ($sample, $opt) = @_;

    my %name_map = reverse %{parse($opt, { required => 0 })};
    $name_map{$sample} //= "sample";
    my $portal_name = ucfirst ($name_map{$sample} =~ tr/_/ /r);
    return $portal_name;
}

sub linkArtefact {
    my ($source_path, $target_path, $portal_name, $opt) = @_;

    $opt->{PORTAL_LINKS} = {} if not exists $opt->{PORTAL_LINKS};
    not -l $source_path or die "$source_path is a symlink and should not be provided to the portal";
    $opt->{PORTAL_LINKS}->{$portal_name} = $source_path;

    my ($name, $parent_dir) = fileparse($source_path);
    if ($source_path ne $target_path) {
        local $CWD = $parent_dir;
        not -l $target_path or unlink $target_path or confess "Couldn't replace previous symlink $target_path: $!";
        symlink $name, $target_path or confess "Couldn't create symlink $target_path: $!";
    }
    return fileparse($target_path);
}

sub writePortalLinks {
    my ($opt) = @_;
    my $outfile = catfile($opt->{OUTPUT_DIR}, "logs", "portal.json");
    writeJson($outfile, $opt->{PORTAL_LINKS});
    return;
}

1;
