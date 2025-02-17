import pandas as pd
import sys
import gzip
import subprocess

try:
    REGIONS = sys.argv[1]
    PANELS = sys.argv[2]
    THRESHOLD = int(sys.argv[3])
except:
	print("----------\nUsage: python3 check_panels.py REGIONS PANELS THRESHOLD\n----------\n")


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


def subset_regions(df, threshold):
    # subset regions to only those with bases < threshold
    df = df.loc[df['bases'] < threshold]
    return df


def intersect_panels():
    # merge overlapping intervals
    command = f"bedtools intersect -a {PANELS} -b annotated_regions.bed"
    subprocess.run(command, shell=True)


def write_amended_bed(df):
    # write amended data to intermediate bed file
    df.to_csv("annotated_regions.bed", sep="\t", header=False, index=False, columns=['chrom', 'start', 'end', 'length', 'coverage', 'bases'])


def clean_temp_files():
    # remove intermediate files
     subprocess.run("rm annotated_regions.bed", shell=True)


def main():
    # parse regions file
    regions_df = parse_regions_file(REGIONS)
    # calculate bases per interval
    regions_df = calculate_bases_per_interval(regions_df)
    # subset regions to those with bases < threshold
    regions_df = subset_regions(regions_df, THRESHOLD)
    # write amended data to intermediate bed file
    write_amended_bed(regions_df)
    # intersect panels with annotated regions
    intersect_panels()
    # remove intermediate files
    clean_temp_files()


if __name__=="__main__":
    main()
