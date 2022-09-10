# Code to extract matrix from 1kb cooler

# modules
import argparse
import cooler
import numpy as np

parser = argparse.ArgumentParser(description='Extract specific region from 1 kb cooler')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-r', '--region', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('-f', '--flip', action='store_true')
args = parser.parse_args()

# Load cooler and parse coordinates
clr = cooler.Cooler(args.input + '::resolutions/1000')
chrom, coord = args.region.split(':')
start, stop = coord.split('-')
region = (chrom, int(start), int(stop))

# Option to flip if region is inverted
if args.flip:
	np.savetxt(args.output, np.flip(clr.matrix(balance=True).fetch(region)))
else:
	np.savetxt(args.output, clr.matrix(balance=True).fetch(region))