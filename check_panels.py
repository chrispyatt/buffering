import pandas as pd
import sys

try:
    REGIONS = sys.argv[1]
    PANELS = sys.argv[2]
except:
	print("----------\nUsage: python3 check_panels.py REGIONS PANELS\n----------\n")


with open(REGIONS, 'r') as fh:
      df = pd.read_csv(fh, header=None, sep="\t", names=["chrom","start","end","coverage","length","bases"])
      