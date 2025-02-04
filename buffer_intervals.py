'''
Author: Chris Pyatt

Script to iteratively split, merge, & buffer an intervals BED file to conform to given minima and maxima
'''

# read capture bed in

import sys
import subprocess
import pandas as pd

try:
    INFILE = sys.argv[1]
    LOWER_LIMIT = int(sys.argv[2])
    UPPER_LIMIT = int(sys.argv[3])
    OUTFILE_PREFIX = sys.argv[4]
except:
	print("----------\nUsage: python3 buffer_intervals.py INFILE LOWER_LIMIT UPPER_LIMIT OUTFILE_PREFIX\n----------\n")


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
     

def split_and_buffer(infile, lower, upper):
    new_file_contents = ""
    with open(infile, 'r') as fh:
        for line in fh:
            chrom = line.split('\t')[0]
            start = int(line.split('\t')[1])
            end = int(line.split('\t')[2])
            length = end - start
            if length < lower:
                buffer = int((lower - length) / 2)
                new_start = start - buffer
                new_end = end + buffer
                new_file_contents = new_file_contents + (f'{chrom}\t{new_start}\t{new_end}\n')
            elif length > upper:
                num_intervals = round(length / upper)
                if length % upper < lower:
                    num_intervals = num_intervals + 1
                for interval in split_interval(start, end, num_intervals):
                    new_file_contents = new_file_contents + (f'{chrom}\t{interval[0]}\t{interval[1]}\n')
            else:
                new_file_contents = new_file_contents + (f'{chrom}\t{start}\t{end}\n')

    new_file_contents = new_file_contents.strip()
    return new_file_contents


def check_min_max():
    # make df
    with open("merged.bed", "r") as fh:
        df = pd.read_csv(fh, sep="\t", header=None)
    # make length column
    df[3] = df[2] - df[1]
    # return min and max
    return (df.min(axis=0)[3], df.max(axis=0)[3])


def main():
    # round 1
    file_content = split_and_buffer(INFILE, LOWER_LIMIT, UPPER_LIMIT)
    write_targets_bed(file_content)
    merge_overlaps()

    # iterate if needed
    while ( check_min_max()[0] < LOWER_LIMIT or check_min_max()[1] > UPPER_LIMIT ):
        subprocess.run("cat merged.bed | wc -l", shell=True)
        print(check_min_max())
        file_content = split_and_buffer("merged.bed", LOWER_LIMIT, UPPER_LIMIT)
        write_targets_bed(file_content)
        merge_overlaps()

    # create final output
    fname = f"{OUTFILE_PREFIX}_{LOWER_LIMIT}_{UPPER_LIMIT}.bed"
    subprocess.run(f"cp merged.bed {fname}", shell=True)

    # clean up
    clean_temp_files()


if __name__=="__main__":
    main()