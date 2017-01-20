#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

"""
filterFreebayes.py
Applies BCBIO developed filters to the freebayes somatic output.
2 main filters applied:
1. LOD Tumor, LOD Normal > 3.5
2. Freq of Tumor / Freq Normal > 2.7
"""

import argparse
import sys

CHROM_INDEX = 0
POS_INDEX = 1
ALT_PARTS_INDEX = 4
FORMAT_PARTS_INDEX = 8
NORMAL_PARTS_INDEX = 9
TUMOR_PARTS_INDEX = 10

# Thresholds are like phred scores, so 3.5 = phred35
LOD_NORMAL_THRESHOLD = 3.5
LOD_TUMOR_THRESHOLD = 3.5

FREQ_NORMAL_THRESHOLD = 0.001
FREQ_RATIO_THRESHOLD = 2.7


def customFilterFreebayes(vcf_file):
    stripped_lines = (line.strip("\n") for line in vcf_file)
    somatic_lines = (line for line in stripped_lines if check_line(line))
    print("\n".join(somatic_lines))


def check_line(line):
    parts = line.split("\t")
    return line.startswith("#") \
        or (parts[ALT_PARTS_INDEX] != "." \
            and check_lods(parts, LOD_TUMOR_THRESHOLD, LOD_NORMAL_THRESHOLD) \
            and check_freqs(parts, FREQ_NORMAL_THRESHOLD, FREQ_RATIO_THRESHOLD))


def check_lods(parts, tumor_threshold, normal_threshold):
    """
    Ensure likelihoods for tumor and normal pass thresholds.

    Skipped if no FreeBayes GL annotations available.
    """
    try:
        gl_index = parts[FORMAT_PARTS_INDEX].split(":").index("GL")
    except ValueError:
        print("skipping LOD check, no GL in '{}'".format(parts[FORMAT_PARTS_INDEX]), file=sys.stderr)
        return True
    try:
        tumor_gls = [float(x) for x in parts[TUMOR_PARTS_INDEX].split(":")[gl_index].split(",") if x != "."]
        if tumor_gls:
            tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
        else:
            tumor_lod = -1.0
    # No GL information, no tumor call (so fail it)
    except IndexError as e:
        tumor_lod = -1.0
        print("assigning {} to tumor LOD for '{}' due to {}".format(tumor_lod, parts[TUMOR_PARTS_INDEX], e), file=sys.stderr)
    except Exception as e:
        print("{} for parts:\n{}".format(e, "\n".join(parts)), file=sys.stderr)
        raise
    try:
        normal_gls = [float(x) for x in parts[NORMAL_PARTS_INDEX].split(":")[gl_index].split(",") if x != "."]
        if normal_gls:
            normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
        else:
            normal_lod = normal_threshold
    # No GL information, no normal call (so pass it)
    except IndexError as e:
        normal_lod = normal_threshold
        print("assigning {} to normal LOD for '{}' due to {}".format(normal_lod, parts[NORMAL_PARTS_INDEX], e), file=sys.stderr)
    except Exception as e:
        print("{} for parts:\n{}".format(e, "\n".join(parts)), file=sys.stderr)
        raise
    result = normal_lod >= normal_threshold and tumor_lod >= tumor_threshold
    if args.debug and not result:
        print('{} LOD filtered (normal: {}, tumor: {})'.format(location(parts), normal_lod, tumor_lod), file=sys.stderr)
    return result


def calc_freq(item, ro_index, ao_index):
    try:
        if ao_index is not None and ro_index is not None:
            ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
            ro = int(item.split(":")[ro_index])
            freq = ao / (ao + ro)
        else:
            freq = None
            print("failing to calculate frequency for {} due to missing AO/RO".format(item), file=sys.stderr)
    except (IndexError, ValueError, ZeroDivisionError) as e:
        freq = None
        print("assigning {} to frequency for '{}' due to {}".format(freq, item, e), file=sys.stderr)
    return freq


def check_freqs(parts, normal_threshold, ratio_threshold):
    """
    Ensure frequency of tumor to normal passes a reasonable threshold.

    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    try:
        ro_index = parts[FORMAT_PARTS_INDEX].split(":").index("RO")
        ao_index = parts[FORMAT_PARTS_INDEX].split(":").index("AO")
    except ValueError:
        ao_index, ro_index = None, None
    if ao_index is None:
        raise NotImplementedError("Unexpected format annotations: %s" % parts[0])
    tumor_freq = calc_freq(parts[TUMOR_PARTS_INDEX], ro_index, ao_index)
    normal_freq = calc_freq(parts[NORMAL_PARTS_INDEX], ro_index, ao_index)
    result = normal_freq is not None and tumor_freq is not None and (
        normal_freq <= normal_threshold or
        normal_freq <= tumor_freq / ratio_threshold)
    if args.debug and not result:
        print('{} FREQ filtered (normal: {}, tumor: {})'.format(location(parts), normal_freq, tumor_freq), file=sys.stderr)
    return result


def location(item):
    return "CHROM {} POS {}".format(item[CHROM_INDEX], item[POS_INDEX])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required_named = parser.add_argument_group("required named arguments")
    required_named.add_argument("-v", "--vcf_file", nargs="?", type=argparse.FileType('r'), default=sys.stdin, help="path/to/file.vcf")
    optional_named = parser.add_argument_group("optional named arguments")
    optional_named.add_argument("-d", "--debug", action="store_true", help="extra debug logging")

    args = parser.parse_args()
    customFilterFreebayes(args.vcf_file)
