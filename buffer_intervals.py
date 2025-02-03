'''
Author: Chris Pyatt

Script to iteratively split, merge, & buffer an intervals BED file to conform to given minima and maxima
'''

# read capture bed in

import sys
import subprocess

try:
	INFILE = sys.argv[1]
	LOWER_LIMIT = int(sys.argv[2])
	UPPER_LIMIT = int(sys.argv[3])
except:
	print("----------\nUsage: python3 buffer_intervals.py INFILE LOWER_LIMIT UPPER_LIMIT\n----------\n")


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
    subprocess.run(command)


def write_targets_bed(content):
     with open("targets.bed", "w") as fh:
        fh.write(content)


def clean_temp_files():
     subprocess.run("rm targets.bed", "rm merged.bed")
     

def split_and_buffer(infile, lower, upper):
    new_file_contents = ""
    with open(infile, 'r') as fh:
        for line in fh:
            chrom = line.split('\t')[0]
            start = int(line.split('\t')[1])
            end = int(line.split('\t')[2])
            length = end - start
            if length <= lower:
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
     with open("merged.bed", "r") as fh:
          # make df
          # make length column
          # sort length column
          # return min and max
          pass
