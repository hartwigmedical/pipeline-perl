#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import sys
import argparse
import vcf

FORMAT_FIELDS = ["GT", "AD", "DP"]
FORMAT = ":".join(FORMAT_FIELDS)
CALL_DATA_FORMAT = vcf.model.make_calldata_tuple(FORMAT_FIELDS)


def melt_somatic_vcf(in_vcf_path, remove_filtered, tumor_sample):
    tumor_call_parsers = {
        "TUMOR": StrelkaParser(),
    }

    with open(in_vcf_path) as in_fh:
        in_vcf = vcf.Reader(in_fh, strict_whitespace=True)
        out_vcf = reheader_vcf(in_vcf, tumor_sample)
        first_record = next(in_vcf)
        validate_sample_names(first_record, tumor_call_parsers.keys(), in_vcf_path)
        melt_record(first_record, tumor_call_parsers, tumor_sample, remove_filtered, out_vcf)
        for record in in_vcf:
            melt_record(record, tumor_call_parsers, tumor_sample, remove_filtered, out_vcf)


def melt_record(record, tumor_call_parsers, tumor_sample, remove_filtered, out_vcf):
    # can be record.is_filtered() when PyVCF releases again
    if remove_filtered and record.FILTER is not None and len(record.FILTER) != 0:
        return
    support = [0] * len(record.alleles)
    ads = [[] for _ in xrange(len(record.alleles))]
    dps = []
    for call in record.samples:
        if call.sample not in tumor_call_parsers or call.gt_alleles is None or not any(call.gt_alleles):
            continue
        dps.append(call.data.DP)
        for allele_num in set(int(an) for an in call.gt_alleles if an is not None):
            support[allele_num] += 1
        call_parser = tumor_call_parsers[call.sample]
        for allele_num, allele in enumerate(record.alleles):
            try:
                ad = call_parser.ad(allele_num, allele, call)
                if ad is not None:
                    ads[allele_num].append(ad)
            except Exception:
                print("exception for allele num {} of {}, {}, {}".format(allele_num, record, record.INFO, call), file=sys.stderr)
                raise

    if sum(support[1:]) == 0:
        return

    melted_record = make_melted_record(record, tumor_sample, support, ads, dps)
    out_vcf.write_record(melted_record)


def make_melted_record(record, tumor_sample, support, ads, dps):
    gts = [str(allele_num) for allele_num, count in enumerate(support) if count != 0]
    melted_gt = "/".join(gts * 2 if len(gts) == 1 else gts)
    melted_ad = ",".join(str(int(round(sum(ds) / len(ds))) if len(ds) > 0 else 0) for ds in ads)
    melted_dp = int(round(sum(dps) / len(dps)))

    record.FORMAT = FORMAT
    data = CALL_DATA_FORMAT(melted_gt, melted_ad, melted_dp)

    melted_call = vcf.model._Call(record, tumor_sample, data)
    record.samples = [melted_call]
    return record


def validate_sample_names(record, required_sample_names, in_vcf_path):
    sample_names = [sample.sample for sample in record.samples]
    missing = [sample for sample in required_sample_names if sample not in sample_names]
    if missing:
        sys.exit("Error: {0} somatic samples missing from {1}. Samples found: {2}".format(
            ", ".join(missing),
            in_vcf_path,
            record.samples,
        ))


def reheader_vcf(in_vcf, tumor_sample):
    if "source" not in in_vcf.metadata:
        in_vcf.metadata["source"] = []
    in_vcf.metadata["source"].append(sys.argv[0])
    # writer outputs headers from constructor
    orig_samples = in_vcf.samples
    in_vcf.samples = [tumor_sample]
    out_vcf = vcf.Writer(sys.stdout, template=in_vcf)
    in_vcf.samples = orig_samples
    return out_vcf


class StrelkaParser(object):
    ##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
    ##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
    ##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
    ##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
    ##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2">
    ##FORMAT=<ID=TOR,Number=2,Type=Integer,Description="Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2">
    ##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">
    def ad(self, allele_num, allele, call):
        if not hasattr(call.data, "TAR"):
            snp_field_name = "{}U".format(allele)
            return getattr(call.data, snp_field_name)[0] if hasattr(call.data, snp_field_name) else None
        elif allele_num == 0:
            return call.data.TAR[0]
        elif len(allele) > 1 or len(call.site.REF) > 1:
            return call.data.TIR[0]
        else:
            return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A tool to combine a subset of fields from multi-sample VCFs into a single-sample VCF.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--remove_filtered", action="store_true", help="Skip variants with a filter flag other than PASS.")
    parser.add_argument("-v", "--vcf_file", required=True, help="path/to/file.vcf")
    parser.add_argument("-t", "--tumor_sample", required=True, help="Tumor sample name.")
    args = parser.parse_args()

    melt_somatic_vcf(args.vcf_file, args.remove_filtered, args.tumor_sample)
