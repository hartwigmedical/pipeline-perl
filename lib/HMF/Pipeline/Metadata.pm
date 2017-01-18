package HMF::Pipeline::Metadata;

use FindBin::libs;
use discipline;

use Carp;
use File::Spec::Functions qw(:ALL);
use JSON;

use HMF::Pipeline::Config qw(sampleBamAndJobs);

use parent qw(Exporter);
our @EXPORT_OK = qw(
    parse
    linkArtefact
    linkExtraArtefact
    linkVcfArtefacts
    linkBamArtefacts
    metaSampleName
    sampleControlNames
    readLinks
    writeLinks
);


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

    my $json = to_json($data, {utf8 => 1, pretty => 1, canonical => 1});
    open my $json_fh, ">:encoding(UTF-8)", $path or confess "Can't open $path: $!";
    print $json_fh $json;
    close $json_fh;
    return;
}

sub parse {
    my ($opt, $parse_opts) = @_;
    $parse_opts //= {};

    my $metadata_path = catfile($opt->{OUTPUT_DIR}, "metadata");
    if (not -s $metadata_path and exists $parse_opts->{required} and not $parse_opts->{required}) {
        return {};
    } else {
        return readJson($metadata_path);
    }
}

sub metaSampleName {
    my ($sample, $opt) = @_;

    my $metadata = parse($opt, {required => 0});
    defined $metadata->{$_} or delete $metadata->{$_} for keys %{$metadata};
    my %name_map = reverse %{$metadata};
    $name_map{$sample} //= "sample";
    return $name_map{$sample};
}

sub sampleControlNames {
    my ($opt) = @_;

    my $metadata = parse($opt);
    my $ref_sample = $metadata->{ref_sample} or die "metadata missing ref_sample";
    my $tumor_sample = $metadata->{tumor_sample} or die "metadata missing tumor_sample";
    my $joint_name = "${ref_sample}_${tumor_sample}";
    return ($ref_sample, $tumor_sample, $joint_name);
}

sub linkArtefact {
    my ($source_path, $canonical_name, $opt) = @_;

    $opt->{LINKS} = {} if not exists $opt->{LINKS};
    $opt->{LINKS}->{$canonical_name} = $source_path;
    return;
}

sub linkExtraArtefact {
    my ($source_path, $opt) = @_;

    $opt->{EXTRAS} = [] if not exists $opt->{EXTRAS};
    push @{$opt->{EXTRAS}}, abs2rel($source_path, $opt->{OUTPUT_DIR});
    return;
}

sub linkVcfArtefacts {
    my ($source_path, $vcf_type, $opt) = @_;

    linkArtefact($source_path, "${vcf_type}_vcf", $opt);
    linkArtefact("${source_path}.idx", "${vcf_type}_vcf_index", $opt);
    return;
}

sub linkBamArtefacts {
    my ($opt) = @_;

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my ($bam_path) = HMF::Pipeline::Config::sampleBamAndJobs($sample, $opt);
        my $sample_name = metaSampleName($sample, $opt);
        linkArtefact($bam_path, "${sample_name}_bam", $opt);
        linkArtefact("${bam_path}.bai", "${sample_name}_bai", $opt);
    }
    return;
}

sub linksPath {
    my ($opt) = @_;

    return catfile($opt->{OUTPUT_DIR}, "logs", "links.json");
}

sub readLinks {
    my ($opt) = @_;

    my $links_path = linksPath($opt);
    return -s $links_path ? readJson($links_path) : {};
}

sub writeLinks {
    my ($opt) = @_;

    writeJson(linksPath($opt), $opt->{LINKS});
    return;
}

1;
