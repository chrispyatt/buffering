'''
Author: Chris Pyatt

Script to iteratively split, merge, & buffer an intervals BED file to conform to given insert size and maximum length
'''

# read capture bed in

import sys
import subprocess
import pandas as pd

try:
    INFILE = sys.argv[1]
    INSERT_SIZE = int(sys.argv[2])
    UPPER_LIMIT = int(sys.argv[3])
    OUTFILE_PREFIX = sys.argv[4]
except:
	print("----------\nUsage: python3 buffer_intervals.py INFILE INSERT_SIZE UPPER_LIMIT OUTFILE_PREFIX\n----------\n")


def split_interval(chrom, start, end):
    length = end - start
    num_intervals = round(length / UPPER_LIMIT)
    if length % UPPER_LIMIT < (INSERT_SIZE * 2):
        num_intervals = num_intervals + 1
    # prevent infinite loop where num_intervals is 1 even though length > upper
    if num_intervals == 1:
        num_intervals = 2
    boundary = round(length / num_intervals)
    intervals = [(start, start+boundary)]
    while intervals[-1][1]+1 < end:
        interval_start = intervals[-1][1] + 1
        interval_end = interval_start + boundary
        intervals.append((interval_start, interval_end))
    df = pd.DataFrame(intervals, columns=['buffered_start', 'buffered_end'])
    df['chrom'] = chrom
    df['start'] = None
    df['end'] = None
    return df


def merge_overlaps():
    # merge overlapping intervals
    command="bedtools merge -d -1 -i targets.bed > merged.bed"
    subprocess.run(command, shell=True)


def write_targets_bed(content):
    # write amended data to intermediate bed file
     with open("targets.bed", "w") as fh:
        content.to_csv(fh, sep="\t", header=False, columns=['chrom', 'buffered_start', 'buffered_end'])


def clean_temp_files():
    # remove intermediate files
     subprocess.run("rm targets.bed; rm merged.bed", shell=True)
     

def split_and_buffer(infile, buffer, upper):
    # read in bed file
    with open(infile, 'r') as fh:
        df = pd.read_csv(fh, sep="\t", header=None, names=["chrom", "start", "end"])
    df['chrom'] = df['chrom'].astype(str)
    # buffer all intervals
    if buffer < 1:
        buffer = 1
    df['buffered_start'] = df['start'] - buffer
    df['buffered_end'] = df['end'] + buffer
    # split larger intervals
    df_long = df[df['buffered_end'] - df['buffered_start'] > upper]
    df_ok = df[df['buffered_end'] - df['buffered_start'] <= upper]
    split_intervals = df_long.apply(lambda x: split_interval(x.chrom, x.start, x.end), axis=1)
    # add split intervals back onto those below upper limit
    print('commencing interval loop')
    for interval_df in split_intervals:
        print(interval_df['buffered_start'][0])
        df_ok = pd.concat([df_ok, interval_df])
    # sort the dataframe & reset indices
    df_ok = df_ok.sort_values(['chrom', 'buffered_start', 'buffered_end']).reset_index().drop('index', axis=1)
    # return sorted dataframe
    return df_ok


def get_stats(bed_file):
    # make df
    with open(bed_file, "r") as fh:
        df = pd.read_csv(fh, sep="\t", header=None)
    # make length column
    df[3] = df[2] - df[1]
    # calculate stats
    min = df.min(axis=0)[3]
    max = df.max(axis=0)[3]
    mean = "{:.0f}".format(df[3].mean(axis=0))
    median = "{:.0f}".format(df[3].median(axis=0))
    std_dev = "{:.0f}".format(df[3].std(axis=0))
    # return stats
    return (min, max, mean, median, std_dev)


def main():
    # print stats of input file
    stats = get_stats(INFILE)
    print(
        "Input file stats:\n"
        f"\tMin: {stats[0]}\n"
        f"\tMax: {stats[1]}\n"
        f"\tMean: {stats[2]}\n"
        f"\tMedian: {stats[3]}\n"
        f"\tStd Deviation: {stats[4]}\n"
    )
    # round 1
    file_content = split_and_buffer(INFILE, INSERT_SIZE, UPPER_LIMIT)
    write_targets_bed(file_content)
    merge_overlaps()
    # iterate if needed
    count = 0
    while ( get_stats("merged.bed")[1] > UPPER_LIMIT ):
        file_content = split_and_buffer("merged.bed", INSERT_SIZE, UPPER_LIMIT)
        write_targets_bed(file_content)
        merge_overlaps()
        count = count + 1
    # report final output stats
    stats = get_stats("merged.bed")
    print(
        "Output file stats:\n"
        f"\tMin: {stats[0]}\n"
        f"\tMax: {stats[1]}\n"
        f"\tMean: {stats[2]}\n"
        f"\tMedian: {stats[3]}\n"
        f"\tStd Deviation: {stats[4]}\n"
        f"\nIterations required for calculation: {count}\n"
    )
    # create final output
    fname = f"{OUTFILE_PREFIX}_{INSERT_SIZE}_{UPPER_LIMIT}.bed"
    subprocess.run(f"cp merged.bed {fname}", shell=True)
    # clean up
    clean_temp_files()


if __name__=="__main__":
    main()
