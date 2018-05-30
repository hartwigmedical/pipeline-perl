package HMF::Pipeline::Functions::Metadata;

use FindBin::libs;
use discipline;

use Carp;
use File::Basename;
use File::Spec::Functions qw(:ALL);
use JSON;

use HMF::Pipeline::Functions::Config;

use parent qw(Exporter);
our @EXPORT_OK = qw(
    parse
    linkArtefact
    linkVcfArtefacts
    linkBamArtefacts
    metaSampleName
    refSampleName
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
    open my $json_fh, ">:encoding(UTF-8)", $path or confess "Can't open $path: $!";
    if ($data) {
        my $json = to_json($data, {utf8 => 1, pretty => 1, canonical => 1});
        print $json_fh $json;
    } else {
        print $json_fh "{\n}";
    }
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

sub refSampleName {
    my ($opt) = @_;
    my $metadata = parse($opt);
    my $ref_sample = $metadata->{ref_sample} or die "metadata missing ref_sample";
    return ($ref_sample);
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

    $opt->{LINKS} = readLinks($opt) if not exists $opt->{LINKS};
    $source_path = abs2rel($source_path, $opt->{OUTPUT_DIR}) if file_name_is_absolute($source_path);
    $opt->{LINKS}->{$canonical_name} = $source_path;
    return;
}

sub linkVcfArtefacts {
    my ($source_path, $vcf_type, $opt) = @_;

    linkArtefact($source_path, "${vcf_type}_vcf", $opt);
    linkArtefact("${source_path}.tbi", "${vcf_type}_vcf_index", $opt);
    return;
}

sub linkBamArtefacts {
    my ($opt) = @_;

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my ($bam_path) = HMF::Pipeline::Functions::Config::sampleBamAndJobs($sample, $opt);
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
    my $links = -s $links_path ? readJson($links_path) : {};
    while (my ($name, $link) = each %{$links}) {
        $links->{$name} = stripPath($link, $opt->{RUN_NAME}) if file_name_is_absolute($link);
    }
    return $links;
}

sub writeLinks {
    my ($opt) = @_;

    writeJson(linksPath($opt), $opt->{LINKS});
    return;
}

sub stripPath {
    my ($path, $last_irrelevant_segment) = @_;

    my $relevant_part = basename($path);
    $path = dirname($path);
    while (basename($path) ne $last_irrelevant_segment and not file_name_is_absolute($relevant_part)) {
        $relevant_part = catfile(basename($path), $relevant_part);
        $path = dirname($path);
    }
    return $relevant_part;
}

1;
