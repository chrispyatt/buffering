'''
Author: Chris Pyatt

Script that takes mosdepth output regions bed files, calculates total bases per interval, and plots by chromosome.
Assumes all regions bed files are from the same sample and are named regions_[limits].bed.gz, where the limits are
the upper and lower bounds of the interval lengths, e.g. regions_600-1500.bed.gz
'''


import matplotlib.pyplot as plt
import pandas as pd
import sys
import gzip
import re

# regions files must be named "regions_[limits].bed.gz"
# e.g. "regions_600-1500.bed.gz" where 600 is the lower limit of
# interval lengths associated with that file


def parse_regions_file(regions_file):
	# read regions bed into dataframe & ensure chromosome column is string
	with gzip.open(regions_file, 'rb') as fh:
		df = pd.read_csv(fh, sep="\t", header=None, names=["chrom","start","end","coverage"])
	df['chrom'] = df['chrom'].astype(str)
	return df


def calculate_bases_per_interval(regions_df):
	# add extra columns for interval length (end - start) and bases per interval (length * coverage)
	regions_df['length'] = regions_df['end'] - regions_df['start']
	regions_df['bases'] = regions_df['length'] * regions_df['coverage']
	return regions_df


def make_plot(df):
	limits = df.limits.unique()
	chroms = df.chrom.unique()
	fig, axes = plt.subplots(nrows=8, ncols=3, sharex=True)
	for idx, ax in enumerate(axes.flatten()):
		chrom = chroms[idx]
		df_chrom = df.loc[df.chrom == chrom]
		ax.hist(
			[df_chrom.loc[df_chrom.limits == x, 'bases'] for x in limits],
			label=limits,
			bins=100,
			)
		ax.set_xlim(-5000, 200000)
		ax.title.set_text(chrom)
		ax.grid(True, axis='both')
	plt.suptitle(f"TWE Bases-per-interval when intervals constrained to various limits")
	plt.tight_layout(h_pad=-1)
	handles, labels = plt.gca().get_legend_handles_labels()
	fig.legend(handles, labels, loc='upper right')


def main():
	dfs = []
	for file in sys.argv[1:]:
		limits = re.split("_|\.", file)[1]
		df = parse_regions_file(file)
		df_bases = calculate_bases_per_interval(df)
		df_bases['limits'] = limits
		dfs.append(df_bases)


	big_df = pd.concat(dfs)

	make_plot(big_df)
	plt.show()


if __name__=="__main__":
    main()
