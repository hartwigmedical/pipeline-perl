#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import os
import sys
import vcf
import argparse

TUMOR_SAMPLE = "TUMOR"
THRESHOLD = 1.3


def filter_strelka(input_fh):
    input_vcf = vcf.Reader(input_fh, strict_whitespace=True)
    output_vcf = vcf.Writer(sys.stdout, template=input_vcf)
    passed_records = (record for record in input_vcf if check_record(record))
    for record in passed_records:
        output_vcf.write_record(record)


def check_record(record):
    if len(record.ALT) > 1:
        print("WARN: more than one alt, keeping: {} {} {}".format(record, record.INFO, record.samples), file=sys.stderr)
        return True
    try:
        qs, tier = quality_score(record)
        for call in record.samples:
            if call.sample != TUMOR_SAMPLE:
                continue
            return allelic_frequency(call, tier - 1) * qs > THRESHOLD
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
        return getattr(call.data, alt_field_name)[tier_index] \
            / (getattr(call.data, alt_field_name)[tier_index] + getattr(call.data, ref_field_name)[tier_index])
    else:
        return call.data.TIR[tier_index] / (call.data.TIR[tier_index] + call.data.TAR[tier_index])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--vcf_file", nargs="?", type=argparse.FileType('r'), default=sys.stdin, help="path/to/file.vcf")
    args = parser.parse_args()
    filter_strelka(args.vcf_file)
