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


def split_interval(start, end, num):
    length = end - start
    boundary = round(length / num)
    intervals = [(start, start+boundary)]
    while intervals[-1][0]+boundary+1 < end:
        interval_start = intervals[-1][0] + boundary + 1
        interval_end = interval_start + boundary
        intervals.append((interval_start, interval_end))
    return intervals


def merge_overlaps():
    command="bedtools merge -d -1 -i targets.bed > merged.bed"
    subprocess.run(command, shell=True)


def write_targets_bed(content):
     with open("targets.bed", "w") as fh:
        fh.write(content)


def clean_temp_files():
     subprocess.run("rm targets.bed; rm merged.bed", shell=True)
     

def split_and_buffer(infile, buffer, upper):
    new_file_contents = ""
    with open(infile, 'r') as fh:
        for line in fh:
            chrom = line.split('\t')[0]
            start = int(line.split('\t')[1])
            end = int(line.split('\t')[2])
            # buffer all intervals
            if buffer < 1:
                buffer = 1
            buffered_start = start - buffer
            buffered_end = end + buffer
            # split larger intervals
            length = buffered_end - buffered_start
            if length > upper:
                num_intervals = round(length / upper)
                if length % upper < (buffer * 2):
                    num_intervals = num_intervals + 1
                # prevent infinite loop where num_intervals is 1 even though length > upper
                elif num_intervals == 1:
                    num_intervals = 2
                for interval in split_interval(buffered_start, buffered_end, num_intervals):
                    new_line = f'{chrom}\t{interval[0]}\t{interval[1]}\n'
            else:
                new_line = f'{chrom}\t{buffered_start}\t{buffered_end}\n'
            new_file_contents = new_file_contents + new_line

    new_file_contents = new_file_contents.strip()
    return new_file_contents


def get_stats(bed_file):
    # make df
    with open(bed_file, "r") as fh:
        df = pd.read_csv(fh, sep="\t", header=None)
    # make length column
    df[3] = df[2] - df[1]
    # calculate stats
    min = df.min(axis=0)[3]
    max = df.max(axis=0)[3]
    mean = df[3].mean(axis=0)
    median = df[3].median(axis=0)
    std_dev = df[3].std(axis=0)
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
