'''
Author: Chris Pyatt

Script that takes mosdepth output regions bed files, calculates total bases per interval, and plots by chromosome.
Assumes all regions bed files are from the same sample and are named regions_[insert size]-[upper limit].bed.gz, 
e.g. regions_250-1500.bed.gz
'''


import matplotlib.pyplot as plt
import pandas as pd
import sys
import gzip
import re
import natsort
import argparse


def parse_arguments():
	parser = argparse.ArgumentParser(description="Plot bases per interval by chromosome")
	parser.add_argument(
		"-e", "--excluded_regions",
		type=str, nargs='+', required=False,
		help="Excluded regions bed files, named to match regions files. E.g. excluded_250-1500.bed"
		)
	parser.add_argument(
		"-r", "--regions_files",
		type=str, nargs='+',
		help="Regions bed files, named to specify limits. E.g. regions_250-1500.bed.gz"
		)
	return parser.parse_args()


def parse_regions_file(regions_file):
	# read regions bed into dataframe & ensure chromosome column is string
	with gzip.open(regions_file, 'rb') as fh:
		df = pd.read_csv(fh, sep="\t", header=None, names=["chrom","start","end","coverage"], dtype={'chrom': 'str'})
	df['chrom'] = df['chrom'].astype(str)
	return df


def calculate_bases_per_interval(regions_df):
	# add extra columns for interval length (end - start) and bases per interval (length * coverage)
	regions_df['length'] = regions_df['end'] - regions_df['start']
	regions_df['bases'] = regions_df['length'] * regions_df['coverage']
	return regions_df


def remove_excluded_regions(df, excluded_regions):
	# remove excluded regions from dataframe
	excluded = pd.read_csv(excluded_regions, sep="\t", header=None, names=["chrom","start","end","name","quality","strand"], dtype={'chrom': 'str'})
	df['chrom'] = df['chrom'].astype(str)
	# merge excluded regions with regions df
	merged = pd.merge(df, excluded, how='outer', indicator=True)
	# keep only regions that are in regions df but not in excluded df
	df = merged.loc[merged['_merge'] == 'left_only'].drop(columns=['_merge', 'name', 'quality', 'strand'])
	return df


def make_plot(df):
	# get lists of limits & chromosomes to use in labels
	limits = df.limits.unique()
	chroms = df.chrom.unique()
	# create empty subplots
	fig, axes = plt.subplots(nrows=8, ncols=3, sharex=True)
	# populate subplots with histograms per chromosome
	for idx, ax in enumerate(axes.flatten()):
		chrom = chroms[idx]
		df_chrom = df.loc[df.chrom == chrom]
		ax.hist(
			[df_chrom.loc[df_chrom.limits == x, 'bases'] for x in limits],
			label=limits,
			bins=100,
			)
		ax.set_xlim(-5000, 200000)
		ax.set_xticks([0, 50000, 100000, 150000, 200000])
		ax.title.set_text(chrom)
		ax.grid(True, axis='both')
	# add legends & titles
	plt.suptitle(f"TWE Bases-per-interval when intervals constrained to various limits")
	plt.tight_layout(h_pad=-1)
	handles, labels = plt.gca().get_legend_handles_labels()
	fig.legend(handles, labels, loc='upper right')


def get_stats(df):
    # calculate stats
	min = "{:.2f}".format(df.min(axis=0)['bases'])
	max = "{:.2f}".format(df.max(axis=0)['bases'])
	mean = "{:.2f}".format(df['bases'].mean(axis=0))
	median = "{:.2f}".format(df['bases'].median(axis=0))
	std_dev = "{:.2f}".format(df['bases'].std(axis=0))
	# return stats
	return (min, max, mean, median, std_dev)



def main():
	# empty list of dataframes
	dfs = []
	# parse arguments
	excluded_files = parse_arguments().excluded_regions
	regions_files = parse_arguments().regions_files
	# turn inputs into dictionaries so we can map the regions files to the excluded files
	regions_dict = {}
	for file in regions_files:
		limits = re.split("_|\.", file)[1]
		regions_dict[limits] = file
	excluded_dict = {}
	if excluded_files:
		for file in excluded_files:
			limits = re.split("_|\.", file)[1]
			excluded_dict[limits] = file
	# make dataframe per regions file given
	for limits, regions_file in regions_dict.items():
		df = parse_regions_file(regions_file)
		df_bases = calculate_bases_per_interval(df)
		df_bases['limits'] = limits
		# remove excluded intervals if provided
		if limits in excluded_dict:
			df_included = remove_excluded_regions(df_bases, excluded_dict[limits])
		else:
			df_included = df_bases
			print(f"No excluded regions file provided for {limits}")
		dfs.append(df_included)
		# print basic stats to console
		stats = get_stats(df_included)
		print(
			f"{limits} df stats:\n"
			f"\tMin: {stats[0]}\n"
			f"\tMax: {stats[1]}\n"
			f"\tMean: {stats[2]}\n"
			f"\tMedian: {stats[3]}\n"
			f"\tStd Deviation: {stats[4]}\n"
		)
	# concatenate dataframes and sort by chromosome
	big_df = pd.concat(dfs)
	big_df_sorted = big_df.sort_values(['chrom', 'start', 'end']).reset_index(drop=True)
	big_df_sorted = big_df_sorted.iloc[natsort.index_humansorted(big_df_sorted.chrom)]
	# generate plot
	make_plot(big_df_sorted)
	plt.show()


if __name__=="__main__":
    main()
