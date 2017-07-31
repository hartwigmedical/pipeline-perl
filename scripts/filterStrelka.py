#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from collections import defaultdict

import argparse
import sys
import vcf

TUMOR_SAMPLE = "TUMOR"
THRESHOLD = 1.3
LC_AF_THRESHOLD = 0.1
LC_QS_THRESHOLD = 20


def filter_strelka(input_fh, bed_content):
    input_vcf = vcf.Reader(input_fh, strict_whitespace=True)
    output_vcf = vcf.Writer(sys.stdout, template=input_vcf)
    passed_records = (record for record in input_vcf if check_record(record, bed_content))
    for record in passed_records:
        output_vcf.write_record(record)


def check_record(record, bed_content):
    if len(record.ALT) > 1:
        print("WARN: more than one alt, keeping: {} {} {}".format(record, record.INFO, record.samples), file=sys.stderr)
        return True
    elif record.ALT[0] is None:
        print("WARN: alt is None, skipping: {} {} {}".format(record, record.INFO, record.samples), file=sys.stderr)
        return False
    try:
        qs, tier = quality_score(record)
        for call in record.samples:
            if call.sample != TUMOR_SAMPLE:
                continue
            return allelic_frequency(call, tier - 1) * qs > THRESHOLD or (is_in_lc_region(bed_content, record.CHROM, record.POS) and lc_filter(call, tier - 1, qs))
    except:
        print("exception for {}, {}, {}".format(record, record.INFO, call), file=sys.stderr)
        raise


def quality_score(record):
    allele = record.ALT[0]
    if len(record.REF) == 1 and len(allele) == 1:
        return record.INFO["QSS_NT"], record.INFO["TQSS_NT"]
    else:
        return record.INFO["QSI_NT"], record.INFO["TQSI_NT"]


def allelic_frequency(call, tier_index):
    allele = call.site.ALT[0]
    if len(call.site.REF) == 1 and len(allele) == 1:
        ref_field_name = "{}U".format(call.site.REF)
        alt_field_name = "{}U".format(allele)
        sumAltFieldRefField = getattr(call.data, alt_field_name)[tier_index] + getattr(call.data, ref_field_name)[tier_index]
        # prevent division by 0
        if sumAltFieldRefField == 0:
            return 0
        else:
            return getattr(call.data, alt_field_name)[tier_index] / sumAltFieldRefField
    else:
        sumTIRandTAR = call.data.TIR[tier_index] + call.data.TAR[tier_index]
        # prevent division by 0
        if sumTIRandTAR == 0:
            return 0
        else:
            return call.data.TIR[tier_index] / sumTIRandTAR


def is_in_lc_region(bed_content, chromosome, position):
    for region in bed_content[chromosome]:
        if region.start <= position <= region.end:
            return False
    return True


def lc_filter(call, tier_index, qs):
    return allelic_frequency(call, tier_index) > LC_AF_THRESHOLD and qs > LC_QS_THRESHOLD


class GenomeRegion:
    def __init__(self, start, end):
        self.start = start
        self.end = end


def load_bed(bed_fh):
    bed_content = defaultdict(list)
    with open(bed_fh, 'r') as bed_file:
        for line in bed_file:
            chromosome, start, end = line.split("\t")
            bed_content[chromosome].append(GenomeRegion(start, end))
    return bed_content

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-b", "--bed_file", required=True, help="path to high confidence bed file")
    parser.add_argument("-v", "--vcf_file", nargs="?", type=argparse.FileType('r'), default=sys.stdin, help="path/to/file.vcf")
    args = parser.parse_args()
    bed_records = load_bed(args.bed_file)
    filter_strelka(args.vcf_file, bed_records)
