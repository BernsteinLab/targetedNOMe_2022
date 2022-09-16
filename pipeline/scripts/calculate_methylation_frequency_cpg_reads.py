#! /usr/bin/env python

"""
Modified by KD, version 0.4
- Modified to work with the cpggpc methylation caller
- Modified to output BED file for CpG model
- Modified to keep strands and reads separate in table
"""
import math
import sys
import csv
import argparse
from collections import namedtuple

def make_key(c, s, e, name, strand):
    return c + ":" + str(s) + ":" + str(e) + ":" + name + ":" + strand

def split_key(k):
    f = k.split(":")
    return (f[0], int(f[1])+1, int(f[2])+2, f[3], f[4])

class SiteStats:
    def __init__(self, g_size, g_seq):
        self.num_reads = 0
        self.posterior_methylated = 0
        self.called_sites = 0
        self.called_sites_methylated = 0
        self.group_size = g_size
        self.sequence = g_seq

def update_call_stats(key, num_called_cpg_sites, is_methylated, sequence):
    if key not in sites:
        sites[key] = SiteStats(num_called_cpg_sites, sequence)

    sites[key].num_reads += 1
    sites[key].called_sites += num_called_cpg_sites
    if is_methylated > 0:
        sites[key].called_sites_methylated += num_called_cpg_sites

parser = argparse.ArgumentParser( description='Calculate methylation frequency at genomic CpG sites')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=1.5)
parser.add_argument('-i', '--input', type=str, required=False)
parser.add_argument('-s', '--split-groups', action='store_true')
args = parser.parse_args()
assert(args.call_threshold is not None)

sites = dict()

if args.input:
    in_fh = open(args.input)
else:
    in_fh = sys.stdin
csv_reader = csv.DictReader(in_fh, delimiter='\t')

for record in csv_reader:
    
    num_sites = int(record['num_motifs']) 
    llr = float(record['log_lik_ratio'])

    # Skip ambiguous call
    if abs(llr) < args.call_threshold:
        continue
    sequence = record['sequence']

    is_methylated = llr > 0
    
    # if this is a multi-cpg group and split_groups is set, break up these sites
    if args.split_groups and record['start'] != record['end']:
        c = record['chromosome'] 
        s = int(record['start'])
        e = int(record['end'])
        name=str(record['read_name'])
        strand=str(record['strand'])

        # find the position of the first CG dinucleotide and potential first GC dinucleotide
        sequence = record['sequence']
        cg_pos = sequence.find("MG")
        first_gc_pos = sequence.find("GC", 0, cg_pos)
        first_gm_pos = sequence.find("GM", 0, cg_pos+1)
        first_cg_pos = cg_pos

        if first_gc_pos > -1:
            if first_gm_pos > -1:
                s = s - min([first_gc_pos, first_gm_pos]) + first_cg_pos
            else:
                s = s - first_gc_pos + first_cg_pos
        elif first_gm_pos > -1:
            s = s - first_gm_pos + first_cg_pos

        while cg_pos != -1:
            key = make_key(c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos, name, strand)
            update_call_stats(key, 1, is_methylated, "split-group")
            cg_pos = sequence.find("MG", cg_pos + 1)
    else:
        key = make_key(record['chromosome'], record['start'], record['end'], record['read_name'], record['strand'])
        update_call_stats(key, num_sites, is_methylated, sequence)

# header
# print("\t".join(["chromosome",  "start", "end", "name", "methylation", "strand"]))

sorted_keys = sorted(sites.keys(), key = lambda x: split_key(x))

for key in sorted_keys:
    if sites[key].called_sites > 0:
        (c, s, e, name, strand) = key.split(":")
        f = float(sites[key].called_sites_methylated) / sites[key].called_sites
        print("%s\t%s\t%s\t%s\t%d\t%s" % (c, int(s), int(e)+1, name, sites[key].called_sites_methylated, strand))