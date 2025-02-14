'''
Author: Chris Pyatt

Script to iteratively split, merge, & buffer an intervals BED file to conform to given insert size and maximum length
'''

# read capture bed in

import sys
import subprocess
import pandas as pd
import time

try:
    INFILE = sys.argv[1]
    INSERT_SIZE = int(sys.argv[2])
    UPPER_LIMIT = int(sys.argv[3])
    OUTFILE_PREFIX = sys.argv[4]
except:
	print("----------\nUsage: python3 buffer_intervals.py INFILE INSERT_SIZE UPPER_LIMIT OUTFILE_PREFIX\n----------\n")


def split_interval(chrom, start, end, buffer, upper):
    length = end - start
    num_intervals = round(length / upper)
    if length % upper < (buffer * 2):
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
    df = pd.DataFrame(intervals, columns=['start', 'end'])
    df['chrom'] = chrom
    return df


def merge_overlaps():
    # merge overlapping intervals
    command = "bedtools merge -d -1 -i targets.bed > merged.bed"
    subprocess.run(command, shell=True)


def write_targets_bed(content):
    # write amended data to intermediate bed file
    content.to_csv("targets.bed", sep="\t", header=False, index=False, columns=['chrom', 'start', 'end'])


def clean_temp_files():
    # remove intermediate files
     subprocess.run("rm targets.bed; rm merged.bed", shell=True)


def buffer_intervals(infile, buffer):
    # read in bed file
    df = pd.read_csv(infile, sep="\t", header=None, names=["chrom", "start", "end"])
    df['chrom'] = df['chrom'].astype(str)
    # buffer all intervals
    if buffer < 1:
        buffer = 1
    df['start'] = df['start'] - buffer
    df['end'] = df['end'] + buffer
    # sort the dataframe & reset indices
    df = df.sort_values(['chrom', 'start', 'end']).reset_index(drop=True)
    # return sorted dataframe
    return df


def split_intervals(df, buffer, upper):
    # split larger intervals
    df_long = df[df['end'] - df['start'] > upper]
    df_ok = df[df['end'] - df['start'] <= upper]
    df_split = pd.concat(df_long.apply(lambda x: split_interval(x.chrom, x.start, x.end, buffer, upper), axis=1).tolist())
    # add split intervals back onto those below upper limit
    df_ok = pd.concat([df_ok, df_split])
    # sort the dataframe & reset indices
    df_ok = df_ok.sort_values(['chrom', 'start', 'end']).reset_index(drop=True)
    # return sorted dataframe
    return df_ok


def get_stats(bed_file):
    # make df
    df = pd.read_csv(bed_file, sep="\t", header=None)
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
        "\nInput file stats:\n"
        f"\tMin: {stats[0]}\n"
        f"\tMax: {stats[1]}\n"
        f"\tMean: {stats[2]}\n"
        f"\tMedian: {stats[3]}\n"
        f"\tStd Deviation: {stats[4]}\n"
    )
    # round 1
    start = time.time()
    buffered_intervals = buffer_intervals(INFILE, INSERT_SIZE)
    file_content = split_intervals(buffered_intervals, INSERT_SIZE, UPPER_LIMIT)
    write_targets_bed(file_content)
    merge_overlaps()
    end = time.time()
    elapsed = "{:.2f}".format(end - start)
    print(f"Time elapsed iteration 1: {elapsed} seconds")
    # iterate if needed
    count = 1
    while (get_stats("merged.bed")[1] > UPPER_LIMIT):
        start = time.time()
        df_merged = pd.read_csv("merged.bed", sep="\t", header=None, names=["chrom", "start", "end"])
        file_content = split_intervals(df_merged, INSERT_SIZE, UPPER_LIMIT)
        write_targets_bed(file_content)
        merge_overlaps()
        count = count + 1
        end = time.time()
        elapsed = "{:.2f}".format(end - start)
        print(f"Time elapsed iteration {count}: {elapsed} seconds")
    # report final output stats
    stats = get_stats("merged.bed")
    print(
        "\nOutput file stats:\n"
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
